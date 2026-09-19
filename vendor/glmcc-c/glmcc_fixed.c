/*
 * GLMCC - Generalized Linear Model for Spike Cross-Correlation
 * C implementation with OpenMP parallelization
 *
 * Based on: Kobayashi et al., Nature Communications, 2019
 * Original Python: https://github.com/NII-Kobayashi/GLMCC
 *
 * Compile: gcc -O3 -fopenmp -o glmcc glmcc.c -lm
 * Usage: ./glmcc <data_dir> <n_neurons> <sim|exp> <GLM|LR>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define WIN       50.0
#define DELTA     1.0
#define NPAR      102
#define NPAR_M2   (NPAR-2)    /* 100 */
#define MAX_ITER  1000
/* T_SEC was a compile-time constant of 5400 s, matching the recordings in
 * Kobayashi et al. It is used only to filter spikes to [0, T_SEC*1000) ms and to
 * name the output file -- it does not enter the rate baseline, which comes from
 * the correlogram itself. On recordings longer than 5400 s, or not starting at 0,
 * the filter silently discards data: this dataset spans 10,000-19,800 s starting
 * near 6,000 s, so correct millisecond input left ZERO spikes. Now a CLI argument. */
static double T_SEC = 5400.0;

/* Exponential integral Ei(x) */
static double expi(double x) {
    if (x == 0.0) return -1e30;
    double sum = 0.0, term = 1.0;
    double absx = fabs(x);
    for (int n = 1; n < 300; n++) {
        term *= x / n;
        double add = term / n;
        sum += add;
        if (fabs(add) < 1e-15 * (1.0 + fabs(sum))) break;
    }
    return 0.5772156649015329 + log(absx) + sum;
}

static inline double func_f(double sec, int delay, double tau) {
    return (sec >= delay) ? exp(-(sec - delay) / tau) : 0.0;
}

static inline double k_delta(int i, int j) {
    return (i == j) ? 1.0 : 0.0;
}

/* ---- per-pair workspace (heap-allocated) ---- */
typedef struct {
    double par[NPAR], new_par[NPAR];
    double Gk[NPAR_M2], new_Gk[NPAR_M2];
    double grad[NPAR];
    double hessian[NPAR * NPAR];
    double delta[NPAR];
} Workspace;

static void init_par(double *par, double rate) {
    double lr = log(rate);
    for (int i = 0; i < NPAR; i++) par[i] = lr;
    par[NPAR-2] = 0.1;
    par[NPAR-1] = 0.1;
}

static void calc_Gk(const double *par, double tau0, double tau1,
                    int delay_synapse, double *Gk) {
    for (int i = 0; i < NPAR_M2; i++) {
        double x_k = (i + 1) * DELTA - WIN;
        if (x_k <= -delay_synapse) {
            double arg = par[NPAR-1] * func_f(-x_k, delay_synapse, tau1);
            if (fabs(arg) > 1e-6) {
                double arg2 = par[NPAR-1] * func_f(-x_k + DELTA, delay_synapse, tau1);
                Gk[i] = (expi(arg) - expi(arg2)) * exp(par[i]) * tau1;
            } else {
                Gk[i] = DELTA * exp(par[i]);
            }
        } else if (x_k > delay_synapse) {
            double arg = par[NPAR-2] * func_f(x_k - DELTA, delay_synapse, tau0);
            if (fabs(arg) > 1e-6) {
                double arg2 = par[NPAR-2] * func_f(x_k, delay_synapse, tau0);
                Gk[i] = (expi(arg) - expi(arg2)) * exp(par[i]) * tau0;
            } else {
                Gk[i] = DELTA * exp(par[i]);
            }
        } else {
            Gk[i] = DELTA * exp(par[i]);
        }
    }
}

static double calc_log_post(const double *par, double beta, double tau0, double tau1,
                            const double *c, const double *Gk) {
    double ll = 0.0;
    for (int i = 0; i < NPAR; i++) ll += par[i] * c[i];
    for (int i = 0; i < NPAR_M2; i++) ll -= Gk[i];
    double reg = 0.0;
    for (int i = 0; i < NPAR_M2 - 1; i++) {
        double d = par[i+1] - par[i];
        reg += d * d;
    }
    return ll - (beta / (2.0 * DELTA)) * reg;
}

static void calc_grad(const double *par, double beta, double tau0, double tau1,
                      const double *c, int n_sp, const double *t_sp,
                      int ds, const double *Gk, double *grad) {
    memset(grad, 0, NPAR * sizeof(double));
    for (int i = 0; i < NPAR_M2; i++) {
        double tmp = 0.0;
        grad[i] = -Gk[i];
        if (i == 0) tmp = -(par[0] - par[1]);
        else if (i == NPAR_M2 - 1) tmp = -(par[i] - par[i-1]);
        else tmp = -(par[i] - par[i-1]) - (par[i] - par[i+1]);
        grad[i] += (beta / DELTA) * tmp + c[i];
    }
    double tmp_ij = 0.0, tmp_ji = 0.0;
    for (int s = 0; s < n_sp; s++) {
        if (t_sp[s] > ds) tmp_ij += func_f(t_sp[s], ds, tau0);
        else if (t_sp[s] < -ds) tmp_ji += func_f(-t_sp[s], ds, tau1);
    }
    for (int i = 0; i < NPAR_M2; i++) {
        double x_k = (i + 1) * DELTA - WIN;
        if (x_k > ds) {
            if (fabs(par[NPAR-2]) < 1e-3)
                tmp_ij -= tau0 * exp(par[i]) * func_f(x_k-DELTA, ds, tau0) * (1-exp(-DELTA/tau0));
            else
                tmp_ij -= (tau0*exp(par[i])/par[NPAR-2]) *
                    (exp(par[NPAR-2]*func_f(x_k-DELTA,ds,tau0)) - exp(par[NPAR-2]*func_f(x_k,ds,tau0)));
        } else if (x_k <= -ds) {
            if (fabs(par[NPAR-1]) < 1e-3)
                tmp_ji -= tau1 * exp(par[i]) * func_f(-x_k, ds, tau1) * (1-exp(-DELTA/tau1));
            else
                tmp_ji -= (tau1*exp(par[i])/par[NPAR-1]) *
                    (exp(par[NPAR-1]*func_f(-x_k,ds,tau1)) - exp(par[NPAR-1]*func_f(-x_k+DELTA,ds,tau1)));
        }
    }
    grad[NPAR-2] = tmp_ij;
    grad[NPAR-1] = tmp_ji;
}

static void calc_hess(const double *par, double beta, double tau0, double tau1,
                      const double *c, const double *Gk, int ds, double *H) {
    memset(H, 0, NPAR * NPAR * sizeof(double));
    for (int i = 0; i < NPAR_M2; i++) {
        double x_k = (i + 1) * DELTA - WIN;
        for (int j = 0; j < NPAR; j++) {
            if (j == NPAR-2 && x_k > ds) {
                if (fabs(par[NPAR-2]) < 1e-3)
                    H[i*NPAR+j] = tau0*exp(par[i])*func_f(x_k-DELTA,ds,tau0)*(1-exp(-DELTA/tau0));
                else
                    H[i*NPAR+j] = -(tau0*exp(par[i])/par[NPAR-2]) *
                        (exp(par[NPAR-2]*func_f(x_k-DELTA,ds,tau0)) - exp(par[NPAR-2]*func_f(x_k,ds,tau0)));
            } else if (j == NPAR-1 && x_k <= -ds) {
                if (fabs(par[NPAR-1]) < 1e-3)
                    H[i*NPAR+j] = tau1*exp(par[i])*func_f(-x_k,ds,tau1)*(1-exp(-DELTA/tau1));
                else
                    H[i*NPAR+j] = -(tau1*exp(par[i])/par[NPAR-1]) *
                        (exp(par[NPAR-1]*func_f(-x_k,ds,tau1)) - exp(par[NPAR-1]*func_f(-x_k+DELTA,ds,tau1)));
            } else if (j < NPAR_M2) {
                if (i == j)
                    H[i*NPAR+j] = -Gk[i] + (beta/DELTA)*(k_delta(i,0)+k_delta(i,NPAR_M2-1)-2);
                else
                    H[i*NPAR+j] = (beta/DELTA)*(k_delta(i-1,j)+k_delta(i+1,j));
            }
        }
    }
    for (int i = NPAR-2; i < NPAR; i++) {
        int l = i - (NPAR-2);  /* 0 or 1 */
        double taul = (l==0) ? tau0 : tau1;
        double J = par[i];
        for (int k = 0; k < NPAR_M2; k++) {
            double x_k = (k + 1) * DELTA - WIN;
            int sign = (l==0) ? 1 : -1;
            double sx = sign * x_k;
            if ((l==0 && sx > ds) || (l==1 && sx > ds)) {
                double f0 = func_f(sx-DELTA, ds, taul);
                double f1 = func_f(sx, ds, taul);
                if (fabs(J) < 1e-3) {
                    double tmp = (taul/2)*f0*f0*(1-exp(-2*DELTA/taul));
                    H[i*NPAR+i] -= tmp;
                } else {
                    double t0 = (J*f0-1)*exp(J*f0) - (J*f1-1)*exp(J*f1);
                    H[i*NPAR+i] -= (taul*exp(par[k])/(J*J))*t0;
                }
            }
        }
        for (int j = 0; j < NPAR_M2; j++)
            H[i*NPAR+j] = H[j*NPAR+i];
    }
}

/* Solve Ax=b via Gaussian elimination with partial pivoting */
static int solve_system(double *A, double *b, double *x, int n) {
    double *a = (double*)malloc(n * n * sizeof(double));
    double *bb = (double*)malloc(n * sizeof(double));
    memcpy(a, A, n * n * sizeof(double));
    memcpy(bb, b, n * sizeof(double));
    for (int k = 0; k < n; k++) {
        int maxr = k; double maxv = fabs(a[k*n+k]);
        for (int i = k+1; i < n; i++) if (fabs(a[i*n+k]) > maxv) { maxv=fabs(a[i*n+k]); maxr=i; }
        if (maxv < 1e-20) { free(a); free(bb); return -1; }
        if (maxr != k) {
            for (int j=k;j<n;j++) { double t=a[k*n+j]; a[k*n+j]=a[maxr*n+j]; a[maxr*n+j]=t; }
            double t=bb[k]; bb[k]=bb[maxr]; bb[maxr]=t;
        }
        for (int i=k+1;i<n;i++) {
            double f=a[i*n+k]/a[k*n+k];
            for (int j=k;j<n;j++) a[i*n+j]-=f*a[k*n+j];
            bb[i]-=f*bb[k];
        }
    }
    for (int i=n-1;i>=0;i--) {
        x[i]=bb[i];
        for (int j=i+1;j<n;j++) x[i]-=a[i*n+j]*x[j];
        x[i]/=a[i*n+i];
    }
    free(a); free(bb); return 0;
}

/* Levenberg-Marquardt optimizer. Returns 1 on convergence, 0 on max_iter. */
static int lm_optimize(Workspace *W, double beta, double tau0, double tau1,
                       const double *c, int n_sp, const double *t_sp,
                       int ds, int cond) {
    double C_lm = 0.01, eta = 0.1;
    if (cond > 0) W->par[NPAR - 3 + cond] = 0.0;

    for (int iter = 0; iter <= MAX_ITER; iter++) {
        calc_Gk(W->par, tau0, tau1, ds, W->Gk);
        double lp = calc_log_post(W->par, beta, tau0, tau1, c, W->Gk);
        calc_grad(W->par, beta, tau0, tau1, c, n_sp, t_sp, ds, W->Gk, W->grad);
        calc_hess(W->par, beta, tau0, tau1, c, W->Gk, ds, W->hessian);

        /* Build reduced system if cond > 0 */
        int sn = NPAR;
        double *H = W->hessian, *g = W->grad;
        double Hr[NPAR*NPAR], gr[NPAR];
        int idx = -1;
        if (cond > 0) {
            idx = NPAR - 3 + cond;
            sn = NPAR - 1;
            int ri = 0;
            for (int i = 0; i < NPAR; i++) {
                if (i == idx) continue;
                int rj = 0;
                for (int j = 0; j < NPAR; j++) {
                    if (j == idx) continue;
                    Hr[ri*sn+rj] = H[i*NPAR+j];
                    rj++;
                }
                gr[ri] = g[i];
                ri++;
            }
            H = Hr; g = gr;
        }

        /* Add damping: H += C_lm * diag(diag(H)) */
        for (int i = 0; i < sn; i++)
            H[i*sn+i] += C_lm * H[i*sn+i];

        if (solve_system(H, g, W->delta, sn) < 0) {
            C_lm /= eta;
            continue;
        }

        /* Update */
        if (cond > 0) {
            int ri = 0;
            for (int i = 0; i < NPAR; i++) {
                if (i == idx) { W->new_par[i] = 0.0; continue; }
                W->new_par[i] = W->par[i] - W->delta[ri++];
            }
        } else {
            for (int i = 0; i < NPAR; i++)
                W->new_par[i] = W->par[i] - W->delta[i];
        }
        /* Clamp J */
        for (int i = NPAR-2; i < NPAR; i++) {
            if (W->new_par[i] < -3.0) W->new_par[i] = -3.0;
            if (W->new_par[i] > 5.0)  W->new_par[i] = 5.0;
        }

        calc_Gk(W->new_par, tau0, tau1, ds, W->new_Gk);
        double nlp = calc_log_post(W->new_par, beta, tau0, tau1, c, W->new_Gk);

        if (nlp >= lp) {
            memcpy(W->par, W->new_par, NPAR * sizeof(double));
            C_lm *= eta;
        } else {
            C_lm /= eta;
            continue;
        }
        if (fabs(nlp - lp) < 1e-4) return 1;
    }
    return 0;
}

/* Read spike times from file, filtering [0, T*1000) */
static int read_spikes(const char *path, double *buf, int maxn) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    int n = 0;
    double v, tmax = T_SEC * 1000.0;
    while (n < maxn && fscanf(f, "%lf", &v) == 1)
        if (v >= 0.0 && v < tmax) buf[n++] = v;
    fclose(f);
    return n;
}

/* Build cross-correlogram. Returns spike count. */
static int cross_corr(const double *c1, int n1, const double *c2, int n2,
                      double *t_sp, double *hist, int hist_len) {
    int cnt = 0, w = (int)WIN;
    memset(hist, 0, hist_len * sizeof(double));
    int min_idx = 0;
    for (int i = 0; i < n2; i++) {
        double lo = c2[i] - w, hi = c2[i] + w;
        while (min_idx < n1 && c1[min_idx] <= lo) min_idx++;
        for (int j = min_idx; j < n1 && c1[j] < hi; j++) {
            double diff = c1[j] - c2[i];
            if (diff < WIN && diff > -WIN) {
                if (cnt < 5000000) t_sp[cnt++] = diff;
                int bin = (int)((diff + w) / DELTA);
                if (bin >= 0 && bin < hist_len) hist[bin] += 1.0;
            }
        }
    }
    return cnt;
}

int main(int argc, char **argv) {
    if (argc != 5 && argc != 6) {
        fprintf(stderr, "Usage: %s <data_dir> <n> <sim|exp> <GLM|LR> [T_seconds]\n", argv[0]);
        return 1;
    }
    const char *dir = argv[1];
    int N = atoi(argv[2]);
    int is_exp = (strcmp(argv[3], "exp") == 0);
    int is_LR  = (strcmp(argv[4], "LR") == 0);
    if (argc == 6) T_SEC = atof(argv[5]);
    double beta = is_LR ? 10000.0 : 4000.0;
    double tau[2] = {4.0, 4.0};
    int hist_len = (int)(2 * WIN / DELTA);

    fprintf(stderr, "GLMCC C: dir=%s N=%d mode=%s method=%s beta=%.0f\n",
            dir, N, is_exp?"exp":"sim", is_LR?"LR":"GLM", beta);
    fprintf(stderr, "Threads: %d\n", omp_get_max_threads());

    double *W = (double*)calloc(N * N, sizeof(double));
    int total_pairs = N * (N - 1) / 2;
    int done = 0;

    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 1; i < N; i++) {
        /* Thread-local buffers */
        double *c1 = (double*)malloc(50000000 * sizeof(double));
        double *c2 = (double*)malloc(50000000 * sizeof(double));
        double *t_sp = (double*)malloc(50000000 * sizeof(double));
        double *hist = (double*)calloc(hist_len, sizeof(double));
        Workspace *ws = (Workspace*)calloc(1, sizeof(Workspace));
        char p1[1024], p2[1024];

        for (int j = 0; j < i; j++) {
            snprintf(p1, sizeof(p1), "%s/cell%d.txt", dir, i);
            snprintf(p2, sizeof(p2), "%s/cell%d.txt", dir, j);
            int n1 = read_spikes(p1, c1, 5000000);
            int n2 = read_spikes(p2, c2, 5000000);
            if (n1 <= 0 || n2 <= 0) continue;

            int n_sp = cross_corr(c1, n1, c2, n2, t_sp, hist, hist_len);
            if (n_sp == 0) continue;

            double best_par[NPAR], best_lp = 0, best_ll = 0;
            int best_ds = 1;

            if (!is_exp) {
                init_par(ws->par, (double)n_sp / (2*WIN));
                lm_optimize(ws, beta, tau[0], tau[1], hist, n_sp, t_sp, 3, 0);
                memcpy(best_par, ws->par, NPAR*sizeof(double));
                best_ds = 3;
            } else {
                for (int m = 1; m <= 4; m++) {
                    init_par(ws->par, (double)n_sp / (2*WIN));
                    lm_optimize(ws, beta, tau[0], tau[1], hist, n_sp, t_sp, m, 0);
                    double lp = calc_log_post(ws->par, beta, tau[0], tau[1], hist, ws->Gk);
                    if (m == 1 || (!is_LR && lp > best_lp) || (is_LR && lp > best_ll)) {
                        memcpy(best_par, ws->par, NPAR*sizeof(double));
                        best_lp = lp;
                        best_ds = m;
                    }
                }
            }

            /* Connection parameters */
            int nb = (int)(WIN/DELTA);
            double cc0[2]={0,0}, Jmin[2];
            for (int l = 0; l < 2; l++) {
                int mt = (int)(tau[l]+0.1);
                for (int m = 0; m < mt; m++)
                    cc0[l] += exp(best_par[nb + (l==0?best_ds:-best_ds) + (l==0?m:-m)]);
                cc0[l] /= mt;
                Jmin[l] = sqrt(16.3/tau[l]/cc0[l]);
                if (tau[l]*cc0[l] <= 10) best_par[NPAR-2+l] = 0;
            }

            double c_E=2.532, c_I=0.612, scale=1.277;
            /* Est_Data.py writes par[NPAR-1] as J_+ and par[NPAR-2] as J_-, then sets
             * W[i][j] = calc_PSP(J_+, Jmin[1]*scale) and W[j][i] = calc_PSP(J_-, Jmin[0]*scale).
             * This port had the two slots the other way round, so each direction received
             * the other's coefficient and threshold. */
            double Jp=best_par[NPAR-1], Jm=best_par[NPAR-2];
            double Wij=0, Wji=0;

            if (is_LR) {
                /* Likelihood ratio test */
                init_par(ws->par, (double)n_sp/(2*WIN));
                lm_optimize(ws, beta, tau[0], tau[1], hist, n_sp, t_sp, best_ds, 1);
                double ll_p = calc_log_post(ws->par, beta, tau[0], tau[1], hist, ws->Gk);
                init_par(ws->par, (double)n_sp/(2*WIN));
                lm_optimize(ws, beta, tau[0], tau[1], hist, n_sp, t_sp, best_ds, 2);
                double ll_n = calc_log_post(ws->par, beta, tau[0], tau[1], hist, ws->Gk);
                double D1 = best_ll - ll_p, D2 = best_ll - ll_n;
                double z_a = 15.14;
                if (2*D1 > z_a) Wij = (Jp>=0) ? c_E*Jp : c_I*Jp;
                if (2*D2 > z_a) Wji = (Jm>=0) ? c_E*Jm : c_I*Jm;
            } else {
                if (Jp > Jmin[1]*scale) Wij = Jp*c_E;
                else if (Jp < -Jmin[1]*scale) Wij = Jp*c_I;
                if (Jm > Jmin[0]*scale) Wji = Jm*c_E;
                else if (Jm < -Jmin[0]*scale) Wji = Jm*c_I;
            }

            #pragma omp critical
            {
                W[i*N+j] = Wij;
                W[j*N+i] = Wji;
                done++;
                if (done % 100 == 0 || done == total_pairs)
                    fprintf(stderr, "\r  %d/%d pairs (%.0f%%)", done, total_pairs, 100.0*done/total_pairs);
            }
        }
        free(c1); free(c2); free(t_sp); free(hist); free(ws);
    }
    fprintf(stderr, "\n");

    /* Write W CSV */
    char out[1024];
    snprintf(out, sizeof(out), "W_py_%.0f.csv", T_SEC);
    FILE *fout = fopen(out, "w");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(fout, "%.6f", W[i*N+j]);
            fprintf(fout, j < N-1 ? ", " : "\n");
        }
    }
    fclose(fout);
    fprintf(stderr, "Output: %s\n", out);
    free(W);
    return 0;
}
