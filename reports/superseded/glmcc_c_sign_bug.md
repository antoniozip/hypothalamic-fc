# C GLMCC sign bug: root cause found and fixed

**Date:** 2026-09-19 · `vendor/glmcc-c/glmcc_fixed.c`

The port is now at parity with Kobayashi's Python reference in all four modes. This records
what the defect was, since two earlier rounds of fixes had guessed wrong about it.

## Root cause: the design vector was two entries short

Python's `GLMCC()` builds a **102**-element design vector: 100 correlogram bins plus two
coupling sufficient statistics.

```python
new_c[NPAR-2] += np.exp((-1)*(t_sp[i]-delay_synapse)/tau[0])   # over lags t >  delay
new_c[NPAR-1] += np.exp((t_sp[i]+delay_synapse)/tau[1])        # over lags t < -delay
```

`calc_log_posterior` then consumes all 102:

```python
for i in range(0, NPAR):
    log_likelihood += (par[i][0]*c[i])
```

The port passed the bare **100**-bin histogram instead. So `calc_log_post`'s loop ran two
elements off the end of the buffer:

```
==866131==ERROR: AddressSanitizer: heap-buffer-overflow
READ of size 8 ... in lm_optimize
allocated by thread T0 here: ... calloc   <- double *hist = calloc(hist_len=100, ...)
```

On the diagnostic fixture the port read `c[100] = 0` and `c[101] = 4.35e-319` where the
reference has `2122.85` and `278.08`.

## Why that inverted every sign

`par[NPAR-2]*c[NPAR-2] + par[NPAR-1]*c[NPAR-1]` are the only terms in the objective that
**reward** a positive coupling. Without them, all that remains for J is `-sum(Gk)`, and Gk
grows with J because a larger coupling predicts more spikes. The objective therefore fell
monotonically in J, so LM's accept test (`nlp >= lp`) preferred a more negative J at every
step and drove both couplings to the `-3.0` clamp:

```
J_c  before fix:  J_+ = -1.743017   J_- = -3.000000   <- at the clamp
J_py             :  J_+ =  0.091939   J_- =  4.033509
```

The gradient was never affected — `calc_grad` recomputes the data term from `t_sp` directly —
which is why magnitudes still correlated +0.90 with the reference while the signed
correlation was -0.86. The fit of the baseline rate was fine; only J's sign was wrong.

The two previous reports proposed the correlogram lag convention and the synaptic basis sign
as candidates. Both were correct as written; neither was the defect.

## Also fixed, found by line-by-line comparison against `glmcc.py`

| | reference | port had |
|---|---|---|
| `calc_hess` d²P/dJ₊² | f at `(-x_k, -x_k+Δ)`, guard `x_k <= -ds` | f at `(-x_k-Δ, -x_k)`, guard `-x_k > ds` |
| `calc_Gk` branch 2 guard | `\|J₋·f(x_k)\| > 1e-6` | `\|J₋·f(x_k-Δ)\| > 1e-6` |
| design-vector binning | exact-integer lag split 0.5/0.5 across two bins | 1.0 into one bin |
| LR: `best_ll` | log-likelihood at the chosen delay | never assigned, always 0 |
| LR: D₁/D₂ | D₂ (cond=2, J₊ pinned) gates `W[i][j]` | swapped |
| LR: test statistic | log **likelihood** | log **posterior** |
| LR: coupling stats | rebuilt for the chosen delay | left at those of m = 4 |
| delay selection | reads the converged log posterior | recomputed against a stale `Gk` |

`Hr`/`gr` also moved off the stack (83 KB inside an OpenMP region) into the heap workspace.

## Verification

`tests/test_glmcc_c_parity.py` runs both implementations on one fixture in all four modes:

| mode | max \|C − Python\| | support | signs vs ground truth |
|---|---|---|---|
| exp/GLM | 1e-6 | identical | both correct |
| sim/GLM | 1e-6 | identical | both correct |
| exp/LR | 2e-6 | identical | both correct |
| sim/LR | 1e-6 | identical | both correct |

On the authors' own `vendor/GLMCC/simulation_data` (20 neurons, 190 pairs, sim/GLM):
max \|C − Python\| = 1e-6, identical support, 31 non-zero, 19 positive in both. 1e-6 is the
CSV print precision, so this is agreement to every digit either side writes.

Red-green checked: rebuilding the pre-fix source and re-running the test gives 11 failures
and exit 1.

AddressSanitizer and UndefinedBehaviorSanitizer are clean over 190 pairs in sim/LR.

## Status

The port is usable. It is ~210x faster than the Python at identical output (0.2 s vs 42.0 s
on the 12-neuron crosscheck fixture).

`vendor/glmcc-c/glmcc.c` was a stale duplicate of the same program carrying the same defect,
and `vendor/glmcc-c/glmcc` — the binary `run_jitter_control.py`, `run_gapped_ongoing*.py` and
`run_glmcc_c_missing.sh` invoke — was built from it. `vendor/glmcc-c/Makefile` now builds both
binary names from the single tracked source, so they cannot drift apart again.
