---
type: results-report
date: 2026-09-19
experiment_line: hypothalamic-functional-connectivity
round: 1
purpose: corrected-pipeline-methods-and-results
status: active
branch: analysis/glmcc-corrected
source_artifacts:
  - results/glmcc_regeneration_summary.csv
  - results/revalidation_summary_lightON_ongoing.csv
  - results/stats_first_level_glmcc.csv
  - results/stats_second_level_glmcc.csv
  - results/glmcc/{rich_club,edge_consistency,centrality,null_models,small_world}.csv
  - results/{gapped_ongoing,aligned_ongoing,jitter_control}_summary.csv
  - results/glmcc_positive_control.csv
  - tests/test_glmcc_c_parity.py
---

# Hypothalamic functional connectivity: methods and results

**Round 1 on the corrected pipeline.** Everything computed before 2026-09-19 16:00 was
discarded; see `reports/superseded/README.md` for why. This document is the only current record
of what was run and what came out.

---

## 1. Executive summary

The C implementation of GLMCC used to produce every connectivity matrix in this project
inverted the sign of every coupling. The cause was a two-element buffer error, found and fixed
on 2026-09-19; the port is now at parity with Kobayashi's reference implementation to the
printed precision, in all four of its modes, on synthetic fixtures, on the reference
`simulation_data`, and on this dataset. All 26 animal × condition matrices and all controls were
recomputed.

The corrected pipeline runs end to end and the condition effect is large and highly significant
in every graph metric.

**A second defect, in the control scripts, was found and fixed on the same day.**
`run_jitter_control.py`, `run_gapped_ongoing.py`, `run_gapped_ongoing_aligned.py` and
`run_glmcc_c_missing.sh` invoked the binary with four arguments and no explicit `T`, on spike
files written in **seconds**. GLMCC reads `WIN = 50`, `DELTA = 1` and `tau = 4` in the file's own
unit, so those runs built a **±50 second** correlogram at **1 second** resolution. Every control
arm was therefore measuring at a timescale three orders of magnitude too coarse. The main
matrices were never affected — `regenerate_glmcc.py` converts to milliseconds and passes the true
duration.

**Re-run correctly, both negative controls pass:**

1. **Jitter (§5.2).** Displacing spikes destroys the edge set at exactly the monosynaptic scale:
   16.8% survive ±5 ms, 5.6% ±10 ms, 0.7% ±25 ms, and a full within-window shuffle leaves 1.1%.
   By 25 ms the edge set is already at the shuffle floor. The aligned-ongoing edges are
   timing-dependent.
2. **Geometry (§5.3, withdrawn).** The claimed sign flip with duty cycle was an artifact of
   comparing a millisecond arm against second arms. In matched units the gapped and aligned arms
   are 70–100% negative, the same direction as continuous `ongoing` at 97–100%.

What survives as a genuine caveat is narrower: **spike count gates detection.** Pairs that
receive an edge have a far larger spike-count product than those that do not (day3 median
1,159,561 vs 1,674, p = 8.5e-297) — GLMCC's `J_min ∝ 1/sqrt(cc₀)` threshold. Firing rate decides
which pairs are testable; timing decides whether a testable pair gets an edge. Since `ongoing`
carries 10–30× more spikes than `lightON`, the condition contrast remains confounded with rate.

A third problem is unchanged: the surrogate validation layer removes 0.0% of `lightON` edges and
4.6% of `ongoing` edges, so "validated" connectivity is essentially raw GLMCC output.

**The rate confound was then tested directly (§5.4) and it accounts for the condition effect.**
Matching observation geometry and per-unit spike count between `lightON` and a light-off control
collapses the ANOVA F for condition from 378→3.9 (node strength), 1049→4.9 (clustering),
533→3.3 (local efficiency) and 409→0.47 (hub score). Nothing survives Holm correction across the
four metrics. A small residual remains — the light-off arm is ~2% denser, reproducible across 9
of 10 thinning seeds — in the opposite direction to the original claim.

**Recommendation: the connectivity is defensible as a timing-based measurement, but the
`lightON` vs `ongoing` contrast as reported is a detection-power artifact and should be retired.**
Section 9 lists what remains.

---

## 2. Experiment identity and decision context

| | |
|---|---|
| Question | Does hypothalamic functional connectivity differ between light-evoked (`lightON`) and spontaneous (`ongoing`) epochs? |
| Subjects | 13 recordings (6 "day", 7 "night"), 126–24 units each, 682 units total |
| Conditions | `lightON` (stimulus epoch), `ongoing` (spontaneous) |
| Estimator | GLMCC (Kobayashi et al., *Nat Commun* 2019) |
| Decision at stake | whether the condition contrast is publishable as a biological result |

The decision context for *this* round is narrower than the scientific question: the previous
body of results was invalid, so this round asks only **"what does the corrected pipeline
actually produce, and does it survive its own controls?"**

---

## 3. Methods

### 3.1 Spike data and preprocessing

Spike times are stored in seconds, spanning roughly 10,000–19,800 s and beginning near 6,000 s.
Both GLMCC implementations expect 0-based **milliseconds** over a recording of length *T*, and
use a ±50-unit correlogram window. Feeding the data unconverted gives a ±50 **second** window;
feeding raw milliseconds trips the hard-coded *T* = 5400 s filter and discards every spike.

`scripts/regenerate_glmcc.py` therefore converts each animal to 0-based milliseconds
(subtracting the global first spike time) and passes the true recording duration (padded 1%) so
nothing is truncated.

Per-dataset counts after preprocessing:

| animal | cond | units | dur (s) | spikes | animal | cond | units | dur (s) | spikes |
|---|---|---|---|---|---|---|---|---|---|
| day1 | lightON | 126 | 10,915 | 169,282 | day1 | ongoing | 126 | 17,783 | 2,238,378 |
| day2 | lightON | 31 | 10,805 | 57,213 | day2 | ongoing | 31 | 16,876 | 676,988 |
| day3 | lightON | 47 | 10,491 | 94,808 | day3 | ongoing | 47 | 16,124 | 964,399 |
| day4 | lightON | 61 | 11,753 | 144,841 | day4 | ongoing | 61 | 17,283 | 1,535,331 |
| day5 | lightON | 71 | 11,260 | 114,298 | day5 | ongoing | 71 | 16,727 | 1,362,752 |
| day6 | lightON | 56 | 11,116 | 170,393 | day6 | ongoing | 56 | 16,801 | 1,942,698 |
| night1 | lightON | 31 | 10,436 | 60,573 | night1 | ongoing | 31 | 16,837 | 748,756 |
| night2 | lightON | 34 | 10,943 | 47,248 | night2 | ongoing | 34 | 17,065 | 585,619 |
| night3 | lightON | 33 | 11,321 | 21,428 | night3 | ongoing | 33 | 18,041 | 276,920 |
| night4 | lightON | 35 | 10,247 | 37,497 | night4 | ongoing | 35 | 17,806 | 388,943 |
| night5 | lightON | 24 | 11,142 | 33,545 | night5 | ongoing | 24 | 18,567 | 406,608 |
| night6 | lightON | 54 | 11,167 | 70,945 | night6 | ongoing | 54 | 17,220 | 973,259 |
| night7 | lightON | 79 | 13,510 | 344,673 | night7 | ongoing | 79 | 19,765 | 3,562,230 |

17,029,625 spikes total. Unit counts match `data/processed/neurons.csv` exactly for all 13
animals in both conditions.

### 3.2 GLMCC estimation

GLMCC fits, to each pair's cross-correlogram, a GLM whose rate is a smooth baseline over 100
1-ms bins spanning ±50 ms plus two exponential coupling terms:

- `WIN` = 50 ms, `DELTA` = 1 ms, `NPAR` = 102 (100 bins + J₋ + J₊)
- τ = 4 ms both directions; β = 4000 (GLM) / 10000 (LR)
- delay scanned over 1–4 ms in `exp` mode, selected by log posterior
- J clamped to [−3, 5]; Levenberg–Marquardt, convergence at |Δ log-posterior| < 1e-4
- a coupling is reported when |J| exceeds `J_min · 1.277`, with
  `J_min = sqrt(16.3 / τ / cc₀)`; PSP = J · 2.532 if excitatory, J · 0.612 if inhibitory

Run in `exp`/`GLM` mode via `vendor/glmcc-c/glmcc_fixed`.

### 3.3 The C port and its validation

The C port is used because it is ~210× faster than the reference Python (0.2 s vs 44.1 s on a
12-neuron fixture), which makes 26 datasets and three control arms tractable.

**It was previously wrong.** `GLMCC()` in the reference builds a 102-element design vector — 100
correlogram bins plus two coupling sufficient statistics,
`Σ_{t>d} exp(−(t−d)/τ₀)` and `Σ_{t<−d} exp((t+d)/τ₁)`. The port passed only the 100-bin
histogram, so `calc_log_post` read two elements past the end of the buffer (confirmed by
AddressSanitizer). Those two terms are the only ones in the objective that reward a positive
coupling; without them the objective falls monotonically in J, and Levenberg–Marquardt drove
every coupling to the −3 clamp. Every weight came out inhibitory.

Seven further divergences from the reference were fixed at the same time (Hessian ∂²/∂J₊² off by
one bin and dropping a boundary bin; a `calc_Gk` guard on the wrong endpoint; missing 0.5/0.5
bin-edge splitting; and four defects that made `LR` mode structurally non-functional).

**Parity is now asserted by `tests/test_glmcc_c_parity.py`,** which runs both implementations on
one fixture in all four modes:

| mode | max \|C − Python\| | support | ground-truth signs |
|---|---|---|---|
| exp/GLM | 1e-6 | identical | both correct |
| sim/GLM | 1e-6 | identical | both correct |
| exp/LR | 2e-6 | identical | both correct |
| sim/LR | 1e-6 | identical | both correct |

1e-6 is the CSV print precision, so this is agreement to every digit either side writes. The
test was red-green verified: rebuilt from the pre-fix source it fails with 11 errors.

Additional checks: on the authors' own `vendor/GLMCC/simulation_data` (20 neurons, 190 pairs),
max |diff| = 1e-6 with identical support and identical sign distribution (31 non-zero, 19
positive in both). On a **real** recording from this dataset (night5 lightON, 24 units, 276
ordered pairs), max |diff| = 1e-6, identical support, 8 negative / 5 positive in both.
AddressSanitizer and UndefinedBehaviorSanitizer are clean over 190 pairs.

`vendor/glmcc-c/Makefile` builds both binary names (`glmcc`, `glmcc_fixed`) from the single
source so the stale-duplicate failure that caused this cannot recur.

> **Note for anyone re-running the reference:** `vendor/GLMCC/Est_Data.py` **crashes on this
> dataset.** A pair with no coincident spikes within ±50 ms gives `rate = 0` and `init_par` calls
> `math.log(0)`. The C port skips such pairs and leaves the weight at zero.

### 3.4 Edge validation

`scripts/revalidate_all.py` masks GLMCC weights by an independent CCG test: cross-correlogram
peaks against 100 shuffle-ISI surrogates, Benjamini–Hochberg step-up at q = 0.05. All 26 runs
completed. Because the p-values come from spike trains, they are independent of the GLMCC
weights.

### 3.5 Graph metrics and statistics

Per-animal metrics (`src/r/compute_graph_metrics.R`, igraph): node strength, clustering
coefficient, local efficiency, hub score. 26/26 computed, 0 skipped, 0 failed.

- **First level:** `lmer(metric ~ condition * region + (1|animal/neuron))`, type-III ANOVA with
  Kenward–Roger degrees of freedom (falls back to type II on failure). n = 1364 observations.
- **Second level:** `lmer(density ~ condition * region_pair * day_night + (1|animal))`,
  type-III / Kenward–Roger.
- **Tier 1:** rich club (degree-preserving nulls), edge consistency across conditions,
  within/between-region decomposition, centrality convergence.
- **Tier 2:** surrogate null models (modularity, clustering, path length), small-world σ and ω.
  Seeded with `set.seed(42)` before any stochastic step.

### 3.6 Controls

| control | construction | purpose |
|---|---|---|
| Positive control | Poisson trains with injected couplings at 3 ms latency | does GLMCC recover known connections at this operating point? |
| Gapped ongoing | `ongoing` gapped to each animal's own lightON duty cycle | separates condition from recording geometry |
| Aligned ongoing | `ongoing` gapped to the real pulse geometry, phase-shifted 45 s into light-off | same, with identical window count/duration/period to lightON |
| Jitter | every spike displaced ±25 ms | destroys 1–5 ms synaptic structure, preserves slow rate structure |

---

## 4. Main findings

### 4.1 Connectivity

26/26 datasets estimated (82 s total). 23,172 edges: 593 positive (2.6%), 22,579 negative.
Density ranges 0.28% (night3 lightON, 3 edges) to 62.67% (night7 ongoing).

Density tracks spike count, as expected from GLMCC's threshold: `J_min ∝ 1/sqrt(cc₀)`, so more
spikes admit smaller couplings. `ongoing` (0.3–3.6 M spikes) runs 16.8–62.7% dense; `lightON`
(21 k–345 k spikes) runs 0.3–22.7%.

### 4.2 Validation removes almost nothing

| condition | raw edges | validated | removed |
|---|---|---|---|
| lightON | 3,179 | 3,179 | **0 (0.0%)** |
| ongoing | 19,993 | 19,067 | 926 (4.6%) |

Median edge-survival across the 26 runs is 0.948 (range 0.647–1.000). The shuffle-ISI null is
close to saturated: it passes nearly every pair, so the masking step is very nearly a no-op.

### 4.3 First-level LME

| metric | term | F | p |
|---|---|---|---|
| node_strength | condition | 378.48 | 3.2e-67 |
| node_strength | region | 2.99 | 2.6e-03 |
| node_strength | condition:region | 14.47 | 1.4e-19 |
| clustering_coefficient | condition | 1049.02 | 2.9e-170 |
| clustering_coefficient | region | 1.13 | 0.34 |
| clustering_coefficient | condition:region | 11.78 | 2.6e-16 |
| local_efficiency | condition | 533.36 | 1.5e-99 |
| local_efficiency | region | 2.60 | 8.0e-03 |
| local_efficiency | condition:region | 3.56 | 4.3e-04 |
| hub_score | condition | 408.51 | 2.4e-71 |
| hub_score | region | 5.12 | 3.3e-06 |
| hub_score | condition:region | 2.02 | 0.042 |

Twelve fixed-effect terms are significant at p < 0.05:

| metric | term | β | p |
|---|---|---|---|
| node_strength | region VM-thalamus | −11.652 | 5.3e-03 |
| node_strength | lightON:VM-thalamus | +17.547 | 7.7e-04 |
| node_strength | lightON:ZI | +22.971 | 2.3e-03 |
| clustering_coefficient | lightON | −0.551 | 3.0e-08 |
| clustering_coefficient | region ARH | −0.259 | 7.9e-03 |
| clustering_coefficient | region ZI | −0.215 | 0.048 |
| clustering_coefficient | lightON:ARH | +0.302 | 0.021 |
| clustering_coefficient | lightON:DMH | +0.270 | 9.0e-03 |
| clustering_coefficient | lightON:PH | +0.216 | 0.039 |
| clustering_coefficient | lightON:ZI | +0.397 | 6.8e-03 |
| local_efficiency | lightON | −0.205 | 9.8e-03 |
| hub_score | lightON | −0.290 | 2.8e-03 |

All four models converged. Residuals are non-normal in every case (Shapiro–Wilk p < 0.001);
clustering_coefficient and node_strength fits are singular in the random effects.

### 4.4 Second-level LME

condition F = 314.94, p < 2e-16; day_night F = 0.0016, p = 0.97;
condition:day_night F = 0.073, p = 0.79.

**The second level supports no pairwise claim.** The model is rank deficient — 76 columns are
dropped for missing `region_pair × day_night` cells — and every emmeans contrast returns
`nonEst`. This is a design problem in the model specification, not a consequence of the
corrected data.

### 4.5 Tier 1

- **Edge consistency: mean 0.075** (max 0.183). Most consistent region pairs: PH↔PVH 0.183,
  VMH↔PH 0.167, PH↔VMH 0.151, DMH↔PH 0.144, MC↔DMH 0.136.
- **Rich club: absent in 9 of 13 animals** (zero significant k-levels). Present in night7 (57
  levels), day4 (11), day5 (11), day6 (2) — the four densest graphs.
- **Region decomposition:** group mean 43.4% within / 56.6% between, with enormous spread
  (day1 97.9% within; day2 21.1%).
- **Centrality convergence:** mean Spearman ρ across strength/betweenness/eigenvector/pagerank
  ranges 0.298 (night3) to 0.939 (day6), tracking density.

### 4.6 Tier 2

**Null models.** Modularity exceeds its degree-preserving null in **3 of 13** animals
(night7 z = 5.08, day5 z = 2.51, day1 z = 1.67); clustering in **2 of 13** (night6 z = 2.40,
night7 z = 2.26). Mean z_modularity 1.28, mean z_clustering 0.43. Path-length z is undefined
everywhere.

**Small-world propensity cannot be computed.** σ and ω require a single connected component, and
every `lightON` graph is disconnected:

| animal | nodes | edges | density | components | largest |
|---|---|---|---|---|---|
| day1 | 126 | 565 | 0.036 | 10 | 117 |
| day2 | 31 | 38 | 0.041 | 16 | 16 |
| day3 | 47 | 110 | 0.051 | 18 | 30 |
| day4 | 61 | 326 | 0.089 | 9 | 53 |
| day5 | 71 | 185 | 0.037 | 16 | 56 |
| day6 | 56 | 386 | 0.125 | 11 | 46 |
| night1 | 31 | 51 | 0.055 | 7 | 25 |
| night2 | 34 | 22 | 0.020 | 18 | 16 |
| night3 | 33 | 3 | 0.003 | 30 | 3 |
| night4 | 35 | 20 | 0.017 | 18 | 18 |
| night5 | 24 | 13 | 0.024 | 14 | 7 |
| night6 | 54 | 59 | 0.021 | 26 | 29 |
| night7 | 79 | 1,401 | 0.227 | 7 | 72 |

---

## 5. Control results

### 5.1 Positive control — the estimator works, but not at `lightON` spike counts

| scenario | spikes/neuron | recall | precision | sign correct |
|---|---|---|---|---|
| generous | 36,448 | 92.5% | 52.1% | **100.0%** |
| ongoing-like | 25,553 | 80.0% | 47.8% | **96.9%** |
| lightON-like | 788 | **0.0%** | 0.0% | — |

GLMCC recovers 80–92% of injected couplings with near-perfect sign accuracy when given
25–36 k spikes/neuron. At 788 spikes/neuron it recovers nothing. The real `lightON` datasets
carry 21 k–345 k spikes across 24–126 units, i.e. roughly 0.6–4.4 k per unit — at or below that
floor for the smaller animals.

### 5.2 Jitter — **this control passes** (corrected)

Superseded numbers: the original ±25 ms control reported 41,827 of 41,842 edges surviving. That
run used a ±50 second correlogram (see §1). Re-run in milliseconds, the same control produces
446 edges against the 2,758-edge aligned-ongoing baseline — ±25 ms removes most of the edge set.

The full sweep is in §5.2a.

### 5.2a Jitter sensitivity sweep — edges break at the monosynaptic scale

Full analysis in `reports/jitter_sensitivity/`; figure `figures/main/figJ1_jitter_sensitivity.png`.
13 animals, 91/91 runs, aligned-ongoing arm, geometry fixed across levels.

| displacement | edges surviving (mean ± SD) | median | range |
|---|---|---|---|
| 0 ms (identity check) | 100.00 ± 0.00 % | 100.00 | — |
| 5 ms | 16.77 ± 17.81 % | 9.09 | 2.44 – 66.67 |
| 10 ms | 5.60 ± 6.81 % | 2.44 | 0.00 – 19.64 |
| 25 ms | 0.73 ± 1.05 % | 0.00 | 0.00 – 3.14 |
| 50 ms | 2.18 ± 2.97 % | 1.22 | 0.00 – 11.11 |
| 100 ms | 0.42 ± 0.51 % | 0.00 | 0.00 – 1.22 |
| **full within-window shuffle (floor)** | **1.11 ± 1.28 %** | 0.99 | 0.00 – 4.55 |

The 0 ms level reproduces the baseline exactly in all 13 animals. **By 25 ms the edge set has
already reached the shuffle floor** — further displacement removes nothing more, because nothing
finer than 25 ms is left. Wilcoxon rejects a 100% null at every level from 5 ms (p = 2.4e-04).

Absolute counts: 2,758 baseline edges → 457 at 5 ms, 225 at 25 ms, 134 after full shuffle; of
those only 242 / 57 / 30 coincide with a baseline edge.

**What still holds from the earlier analysis.** Spike count gates *detection*: pairs with an edge
have a much larger spike-count product than pairs without (night5 4.0×, day3 693×, night6 15×;
p ≤ 2.1e-08). That is GLMCC's `J_min = sqrt(16.3/τ/cc₀)` threshold. Firing rate decides which
pairs are testable; timing decides the outcome for a testable pair.

### 5.3 Geometry — **withdrawn**

This section reported that the same spikes re-windowed to the lightON duty cycle flipped from
97–100% negative to ~100% positive. That compared a millisecond arm against second arms (§1).

With both in milliseconds there is no flip:

| arm | edges | negative |
|---|---|---|
| ongoing, continuous | 19,993 | 97–100% |
| ongoing, gapped to lightON duty cycle | 3,463 | 70.4–100% |
| ongoing, gapped to real pulse geometry (45 s shifted) | 2,758 | 76.2–100% |

All three agree in direction. The duty-cycle geometry reduces the edge count roughly six-fold,
which is expected from the ~10× reduction in observation time, but does not reverse the sign.

No claim in this document now rests on a geometry-driven sign flip.

### 5.4 Rate matching — the condition effect does not survive

Full analysis in `reports/rate_matched/`; figure `figures/main/figR1_rate_matched.png`;
matrices in `results/glmcc_rm/`.

Across the 26 main matrices, Spearman ρ between spikes-per-unit and edge density is **+0.955
(p = 3.7e-14)**, and `ongoing` carries 10.2–13.7× more spikes per unit than `lightON`. GLMCC's
threshold `J_min ∝ 1/sqrt(cc₀)` falls as spikes accumulate, so density follows spike count by
construction.

**Design.** Two arms differing only in whether the light was on: `lightON` inside its true
stimulus windows, and a light-off control built from `ongoing` gapped to those same windows
shifted +45 s. Window count, duration and period are identical, so `T_observed` matches; each
unit is then thinned to the smaller of its two counts, so `n_i` matches per unit. Thinning is
uniform sampling without replacement, which preserves correlogram shape while lowering counts.

**Result.** 26/26 arms estimated. Median density 3.09% (lightON) vs 3.10% (light-off).

| metric | F raw | p raw | F matched | p matched | p Holm |
|---|---|---|---|---|---|
| node_strength | 378.48 | 3.2e-67 | 3.93 | 0.048 | 0.143 |
| clustering_coefficient | 1049.02 | 2.9e-170 | 4.87 | 0.028 | 0.111 |
| local_efficiency | 533.36 | 1.5e-99 | 3.33 | 0.068 | 0.143 |
| hub_score | 408.51 | 2.4e-71 | 0.47 | 0.492 | 0.492 |

F collapses by two to three orders of magnitude; 8 of 12 terms lose significance; **nothing
survives multiple-comparison correction.** The `conditionlightON` fixed effects are all
non-significant (p = 0.49–0.99).

**The residual.** Repeating over 10 thinning seeds, the light-off arm is consistently ~0.07 pp
(~2% relative) denser, significant in 9 of 10. This is not explained by the lightON arm being
thinned harder: the excess is *smaller* where thinning was more asymmetric (ρ = −0.311,
p < 0.001, n = 130), the opposite of that prediction. It is ~2% of the effect size originally
reported and runs in the opposite direction, so it is a different phenomenon and needs the
symmetric-thinning control before interpretation.

### 5.5 Putative E/I class explains nothing

Full analysis in `reports/ei_differential/`; figure `figures/main/figE1_ei_differential.png`.

`{animal}_put_inh.txt` marks narrow-spiking units (half-width < 0.6 ms and trough-to-peak
< 0.65 ms, `tim_stuff.m:714`), a waveform criterion independent of GLMCC. **Only 7 of 13 animals
were ever classified** — those label 38.7–51.1% of units inhibitory, the other six label
0–3.2%, with nothing in between. Analysable set: 399 units, 178 I / 221 E.

| analysis | result |
|---|---|
| firing rate by class | I/E = 1.01×, paired p = 0.33 — **no difference** |
| outgoing coupling sign | E 0.988 vs I 0.984 fraction negative, p = 0.93 — uninformative: 94.9% of units emit only negative couplings |
| node metrics (4) | all p ≥ 0.37 |
| edge types E→E / E→I / I→E / I→E | all Holm p = 1.00; direction holds in only 5–6 of 14 |

The pooled edge-type numbers look orderly (E→E 9.38% > I→I 7.29%) but that is Simpson's paradox
— the per-animal medians reverse the sign, and the paired test is null.

**Two methodological findings outweigh the nulls.** First, `src/py/ei_classification.py` assigns
`"E"` to every unit not listed as inhibitory, so **279 units (41% of the dataset) in the six
unclassified animals are silently labelled excitatory**; `data/processed/neurons_ei.csv` carries
those labels. Use `data/processed/neurons_ei_strict.csv` instead. Second, a narrow-spiking
population that does not fire faster than the broad-spiking one, and that is 44.6% of all units,
is not behaving like a fast-spiking interneuron population — the criterion may not be separating
cell types here at all, which caps what any of these nulls can mean.

### 5.6 Waveform re-classification — the data will not support an E/I split

`reports/ei_differential/waveform-reclassification.md`. Re-derived from `Units(i).wf` (81
samples, 20 kHz) using the spikeMAP method (Mahallati et al., *eLife* 2026): k-means on full
width at half-maximum and peak-to-peak duration, k by Calinski-Harabasz, no fixed thresholds.
That paper reports 5.25% inhibitory in mouse PFC.

| variant | inhibitory |
|---|---|
| pooled k-means | 39.0% |
| pooled, recordings centred | 81.4% |
| per recording, as published | 27.1% (range 3.2–59.6%) |

**None is credible.** The cause is that median peak-to-peak duration varies **0.400–1.050 ms
between recordings**, a 2.6× range exceeding the cell-type separation it should detect. Under
pooled clustering a recording's inhibitory fraction is almost entirely predicted by where it
sits: Spearman **−0.978, p = 8e-09**. Clustering per recording removes that, but then the
"inhibitory" centroid itself ranges 0.222–0.891 ms — 4× — so the classes are not comparable
across animals.

There is some real structure: peak-to-peak is bimodal pooled and after centring (dip p < 1e-4),
and in 6 of 13 recordings individually; full width is unimodal (p = 0.87). Not enough to assign
classes. Agreement with the shipped labels is 70.2%, κ = +0.385.

**A genuine bug in the original code was found:** `tim_stuff.m:702` takes its half-width from
`findpeaks`, which returns *positive* peaks, so it measures the last positive bump rather than
the spike trough — the likely source of the 44.6% inhibitory fraction.

**Diagnosed (§5.7): an acquisition change mid-study.** Ordering by date, median PP splits
cleanly at 2018-01-11 → 2018-01-31 (early n=6: 0.775 ms; late n=7: 0.500 ms; Mann-Whitney
p = 0.0031, Spearman with date order −0.715, p = 0.0060), together with a drop in trough depth
and a rise in positive overshoot — the signature of a stronger high-pass filter. Windowing,
alignment and sampling are identical across all 13 (81 samples, trough at sample 40, no
truncation), and no filter settings are stored in the files.

That also **corrects the §5.5 claim that six animals were never classified**: they were, but the
fixed 0.65 ms threshold sits above the early block's entire distribution. Applying it to
recomputed PP reproduces the 0–3.2% / 41–51% split (Spearman +0.722, p = 0.005).

**Cause confirmed (§5.7).** The nominal acquisition filter was 600–6000 Hz throughout. Sweeping a
compensating high-pass over 50–1525 Hz to equalise median PP, the six early recordings
independently require **550–775 Hz (clustered on ~600)** while the late seven require **0–450 Hz,
mostly ~0**. The 600–6000 Hz filter reached the late block and effectively did not reach the early
block. Compensation cuts the between-recording spread from 2.62× to 1.11× and **removes the block
effect** from the classification (early 21.7% vs late 31.0% inhibitory, p = 0.366, was p = 0.0031).

It does not rescue the classification: still 28.8% inhibitory against a 5.25% reference, still
0–42.6% per animal, full width still unimodal. Compensating a mean, already-filtered waveform is
not equivalent to filtering the raw trace, and **no pre-spike-sorting data survives**, so that
ceiling cannot be lifted. A further untestable consequence: spike amplitude halves across the
boundary, which can affect sorting yield and unit isolation in ways that cannot now be checked.

**Consequence: the §5.5 E/I nulls are untestable rather than negative.** Neither the original
labels nor a literature-standard re-derivation yields a classification worth testing. Resolve
the between-recording waveform variation first, or drop the E/I dichotomy for this dataset.

## 6. What changed our belief

| before this round | after |
|---|---|
| GLMCC does not recover connectivity on data like this (0–5% sign accuracy) | It recovers 80–92% with 97–100% sign accuracy — the earlier failure was the port |
| The `lightON`/`ongoing` sign asymmetry is a duty-cycle artifact of GLMCC's baseline | Not established. The evidence for it compared a millisecond arm against second arms; in matched units the asymmetry is not there |
| Connectivity is small-world with rich-club structure | Small-world is not computable (all graphs disconnected); rich club survives in 4 of 13 |
| Edges are consistent across conditions (0.636) | Edge consistency is 0.075 |
| Surrogate validation constrains the edge set | It removes 0.0% of lightON and 4.6% of ongoing edges |
| *(mid-round)* The edge set carries no timing information and is not synaptic | **Withdrawn.** That rested on control runs at 1 s resolution. Corrected, edges break between 5 and 25 ms |
| The condition effect is the finding | **Tested: the condition effect is the confound.** Matching geometry and spike count collapses F by 100–1000× and nothing survives correction (§5.4) |

### 5.8 File concatenation does not affect the correlograms (checked, negative)

Each recording is a concatenation of 34–35 files with real gaps between them — 2,361 s across a
17,783 s span in day1, **13.3% of the timeline**. That raised the question of whether spike pairs
straddling a file join produce false coincidences, since GLMCC reads one continuous time axis.

**They cannot.** Checked across all 13 recordings:

| check | result |
|---|---|
| spikes falling inside an inter-file gap | **0** in 12 of 13; **1** in night6 (of ~1.5 M) |
| smallest gap, any recording | **128 ms** (day4), vs the ±50 ms correlogram window |

The spike timeline is **absolute**, not re-based: the gaps are genuine periods with no recording
and no spike sits inside one. Every gap exceeds the correlogram window by at least 2.6× (up to
17.6×), so no pair of spikes can straddle a boundary and still fall within ±50 ms. The
concatenation is invisible to GLMCC, which is correct. No correction is needed and none was made.

The 13.3% gap figure is still the source of the "~89% continuous" characterisation of `ongoing`
used in the duty-cycle reasoning; only the coincidence concern is dismissed.

## 7. Limitations and negative results

1. **The headline condition effect does not exist once detection power is equalised** (§5.4).
   This is the principal negative result of the round.
2. **The validation null is inert** — 0.0% of lightON and 4.6% of ongoing edges removed (§4.2).
3. **`lightON` is at or below the detection floor** for several animals (§5.1). night3 lightON
   yields 3 edges from 33 units.
4. **The second-level model is rank deficient** and supports no contrast (§4.4).
5. **LME residuals are non-normal** in all four first-level models; two fits are singular.
6. **Small-world and path-length statistics are undefined** on these graphs (§4.6).
7. **Two implementation defects were found in one day** — a sign inversion in the port and a
   time-base error in every control script. Both produced plausible-looking output for months.
   Anything not re-derived since 2026-09-19 16:00 should be treated as unverified.
8. **Not re-run this round:** community detection, participation coefficient, hub stability,
   circos plots, main manuscript figures, tier-3 analyses.

## 8. Figure index

Only five figures survive this round; the rest were built on discarded data.

| figure | what it shows | status |
|---|---|---|
| `figures/main/figT1_rich_club.png` | φ_norm vs k with degree-preserving nulls | current; supports §4.5 — flat in 9/13 animals |
| `figures/main/figT1_edge_consistency.png` | edge overlap across conditions by region pair | current; supports §4.5 — mean 0.075 |
| `figures/main/figT1_region_decomposition.png` | within- vs between-region edge share | current; supports §4.5 — high variance across animals |
| `figures/main/figT1_centrality_convergence.png` | rank agreement among four centrality measures | current; supports §4.5 |
| `figures/main/figT2_null_models.png` | z-scores for modularity / clustering vs nulls | current; supports §4.6 — 3/13 and 2/13 significant |
| `figures/main/figJ1_jitter_sensitivity.png` | edge survival vs displacement, with shuffle floor and the spike-count mechanism | current; supports §5.2a — edges are timing-dependent |
| `figures/main/figR1_rate_matched.png` | the rate confound, F before/after matching, paired matched densities | current; supports §5.4 — the key negative result |
| `figures/main/figE1_ei_differential.png` | classifiability, absent rate confound, degenerate sign distribution, pooled-vs-paired edge types | current; supports §5.5 |

`figT2_small_world.png` was **not** produced: there is nothing to plot (§4.6).

---

## 9. Next actions

1. ~~Find the scale at which jitter breaks the edge set.~~ **Done (§5.2a).** It breaks between
   5 and 25 ms. The edges are timing-dependent and the synaptic interpretation is available.
2. ~~Break the rate confound.~~ **Done (§5.4).** The condition effect is the confound. Retire the
   raw `lightON` vs `ongoing` contrast; `results/glmcc_rm/` holds the defensible comparison.
3. **Replace the shuffle-ISI null with a within-window shuffle null**, per pair at q = 0.05.
   `scripts/jitter_sensitivity.py` already implements the surrogate. The ~1% shuffle survival in
   §5.2a predicts the real edges will largely pass, so this is a strengthening step rather than a
   demolition — but it is needed before any edge count is quotable.
4. **Regenerate everything derived from the control arms.** The duty-cycle and geometry analyses,
   and any figure or table quoting them, were computed at the wrong time base.
5. **Respecify the second-level model** so it is not rank deficient.
6. **Audit for further unit errors.** `scripts/glmcc_units.py` centralises the conversion;
   `run_glmcc_c_missing.sh` has been disabled. Any new call site must use it.
7. Only then regenerate community/participation/hub-stability/circos and the manuscript figures.

Items 1 and 2 are settled. The condition contrast as framed does not survive. Item 3 decides
whether the edge counts are quotable; the symmetric-thinning control decides whether the small
light-off excess is worth pursuing as a question in its own right.

## 10. Reproducibility index

**Branch:** `analysis/glmcc-corrected`. Discarded artifacts are archived outside the repo (see
the purge commit message); superseded reports are in `reports/superseded/`.

| stage | command | runtime |
|---|---|---|
| Build the port | `make -C vendor/glmcc-c` | <1 s |
| Verify the port | `python3 tests/test_glmcc_c_parity.py` | ~39 s, expect 20 PASS / 0 FAIL |
| Connectivity | `python3 scripts/regenerate_glmcc.py` | 82 s |
| Gapped control | `python3 scripts/run_gapped_ongoing.py` | ~5 min |
| Aligned control | `python3 scripts/run_gapped_ongoing_aligned.py` | ~5 min |
| Jitter control | `python3 scripts/run_jitter_control.py` | ~12 min |
| Positive control | `python3 scripts/glmcc_positive_control.py` | 7 s |
| Rate-matched contrast | `python3 scripts/rate_matched_contrast.py` | ~1 min |
| Matched metrics + LME | `ESTIMATOR=glmcc_rm bash scripts/rebuild_downstream.sh` | ~3 min |
| Matched figure | `python3 scripts/plot_rate_matched.py` | ~10 s |
| E/I differential | `python3 scripts/ei_differential.py && python3 scripts/plot_ei_differential.py` | ~20 s |
| Waveform re-classification | `python3 scripts/ei_classify_waveforms.py` (needs scikit-learn, diptest) | ~40 s |
| Filter harmonisation | `python3 scripts/waveform_harmonise.py` | ~3 min |
| Jitter sensitivity sweep | `python3 scripts/jitter_sensitivity.py` | ~45 min |
| Sweep figure | `python3 scripts/plot_jitter_sensitivity.py` | ~30 s |
| Edge validation | `python3 scripts/revalidate_all.py` | ~35 min |
| Metrics + first level | `bash scripts/rebuild_downstream.sh` | ~3 min |
| Second level | `Rscript src/r/stats_second_level.R --estimator glmcc` | ~1 min |
| Tier 1 | `Rscript src/r/tier1_analyses.R` | ~2 min |
| Tier 2 | `Rscript src/r/tier2_analyses.R` | ~2 min |

**Environment.** Python 3.12 (numpy 2.4.6, scipy 1.15.3); R 4.3.3 (igraph, lme4, lmerTest,
emmeans, pbkrtest); gcc with OpenMP. Seeds: `set.seed(42)` in tier1/tier2;
`np.random.default_rng(42)` in the controls; `default_rng(20260919)` in the parity test.

**`pbkrtest` installation note.** It is required by `stats_second_level.R` for Kenward-Roger
degrees of freedom, and was missing — the script died before writing, leaving a stale CSV on
disk. `Deriv` is not available for R 4.3.3 from the current CRAN index, which blocks `doBy` and
then `pbkrtest`. Install from archive first:

```r
install.packages("https://cloud.r-project.org/src/contrib/Archive/Deriv/Deriv_4.1.3.tar.gz",
                 repos = NULL, type = "source")
install.packages(c("doBy", "pbkrtest"))
```

**Code changed this round.** `vendor/glmcc-c/glmcc_fixed.c` (the port fix), new
`vendor/glmcc-c/Makefile`, new `tests/test_glmcc_c_parity.py`, and a robustness guard in
`src/r/tier2_analyses.R` so it reports "small-world undefined" instead of halting on an empty
result set. No statistic was altered.
