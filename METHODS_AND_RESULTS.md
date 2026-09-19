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
in every graph metric. **But two negative controls fail**, and they fail in a way that undercuts
a synaptic interpretation of the result:

1. **Jitter.** Displacing every spike by ±25 ms destroys all 1–5 ms structure. **41,827 of 41,842
   edges survive it (100.0%).**
2. **Geometry.** The sign of the weights is set by the recording's duty cycle, not by the spike
   trains: continuous `ongoing` yields 97–100% negative weights, the *same spikes* gapped to the
   `lightON` duty cycle yield ~100% positive.

A third, softer problem: the surrogate validation layer removes 0.0% of `lightON` edges and 4.6%
of `ongoing` edges, so "validated" connectivity is essentially raw GLMCC output.

**Recommendation: do not build manuscript claims on these matrices yet.** The estimator is now
correct; what it is measuring on this data is not established. Section 9 lists what would settle it.

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

### 5.2 Jitter — **this control fails**

Each aligned-ongoing matrix compared with its jittered twin (identical window geometry,
identical 45 s offset; only `jitter_ms` differs):

| animal | real | jittered | overlap | surviving |
|---|---|---|---|---|
| day1 | 15,721 | 15,722 | 15,721 | 100.0% |
| day2 | 870 | 870 | 870 | 100.0% |
| day3 | 1,508 | 1,507 | 1,507 | 99.9% |
| day4 | 3,356 | 3,355 | 3,355 | 100.0% |
| day5 | 4,774 | 4,773 | 4,773 | 100.0% |
| day6 | 2,824 | 2,823 | 2,823 | 100.0% |
| night1 | 924 | 924 | 924 | 100.0% |
| night2 | 978 | 977 | 977 | 99.9% |
| night3 | 1,029 | 1,028 | 1,028 | 99.9% |
| night4 | 1,077 | 1,076 | 1,076 | 99.9% |
| night5 | 530 | 531 | 530 | 100.0% |
| night6 | 2,452 | 2,449 | 2,449 | 99.9% |
| night7 | 5,799 | 5,794 | 5,794 | 99.9% |
| **TOTAL** | **41,842** | **41,829** | **41,827** | **100.0%** |

A ±25 ms displacement destroys every 1–5 ms feature a synapse could produce. **Essentially no
edge notices.** Whatever these matrices index, it is not millisecond-scale coupling.

(`retention_pct` in `results/jitter_control_summary.csv` is the fraction of *spikes* retained,
~99.87% — not edge survival. The table above is the edge statistic.)

### 5.3 Geometry — sign is set by duty cycle

| arm | edges | negative |
|---|---|---|
| ongoing, continuous | 19,993 | 97–100% |
| ongoing, gapped to lightON duty cycle | 43,285 | **0.0–0.2%** |
| ongoing, gapped to real pulse geometry (45 s shifted) | 41,842 | **0.0–0.1%** |

The same spikes, re-windowed, flip from almost entirely inhibitory to almost entirely
excitatory, and roughly double in edge count. Both extremes are implausible as biology. The
`lightON`/`ongoing` contrast the study rests on differs in exactly this way.

---

## 6. What changed our belief

| before this round | after |
|---|---|
| GLMCC does not recover connectivity on data like this (0–5% sign accuracy) | It recovers 80–92% with 97–100% sign accuracy — the earlier failure was the port |
| The `lightON`/`ongoing` sign asymmetry is a duty-cycle artifact of GLMCC's baseline | Partly; but the headline contrast compared a Python-derived arm against a defective-C arm, so the effect size was never measured cleanly |
| Connectivity is small-world with rich-club structure | Small-world is not computable (all graphs disconnected); rich club survives in 4 of 13 |
| Edges are consistent across conditions (0.636) | Edge consistency is 0.075 |
| Surrogate validation constrains the edge set | It removes 0.0% of lightON and 4.6% of ongoing edges |
| The condition effect is the finding | The condition effect is large and significant — and is also exactly what the geometry artifact predicts |

---

## 7. Limitations and negative results

1. **The jitter control fails outright** (§5.2). This is the single most important result in
   this document and it is negative.
2. **The sign is a function of recording geometry** (§5.3), and the two conditions differ in
   that geometry by design.
3. **`lightON` is at or below the detection floor** for several animals (§5.1). night3 lightON
   yields 3 edges from 33 units.
4. **Density confounds with spike count**, and `ongoing` has 10–30× more spikes than `lightON`.
   Any unnormalised comparison between conditions inherits this.
5. **The second-level model is rank deficient** and currently supports no contrast (§4.4).
6. **LME residuals are non-normal** in all four first-level models; two fits are singular.
7. **Small-world and path-length statistics are undefined** on these graphs.
8. **Not re-run this round:** community detection, participation coefficient, hub stability,
   circos plots, main manuscript figures, tier-3 analyses. Their previous outputs were deleted
   rather than refreshed, because regenerating them before §5.2 is resolved would waste the effort.

---

## 8. Figure index

Only five figures survive this round; the rest were built on discarded data.

| figure | what it shows | status |
|---|---|---|
| `figures/main/figT1_rich_club.png` | φ_norm vs k with degree-preserving nulls | current; supports §4.5 — flat in 9/13 animals |
| `figures/main/figT1_edge_consistency.png` | edge overlap across conditions by region pair | current; supports §4.5 — mean 0.075 |
| `figures/main/figT1_region_decomposition.png` | within- vs between-region edge share | current; supports §4.5 — high variance across animals |
| `figures/main/figT1_centrality_convergence.png` | rank agreement among four centrality measures | current; supports §4.5 |
| `figures/main/figT2_null_models.png` | z-scores for modularity / clustering vs nulls | current; supports §4.6 — 3/13 and 2/13 significant |

`figT2_small_world.png` was **not** produced: there is nothing to plot (§4.6).

---

## 9. Next actions

**Before any further analysis, settle what the estimator is measuring.** In priority order:

1. **Explain or accept the jitter result.** A ±25 ms jitter that removes no edges means the edge
   set is determined by structure coarser than 25 ms. Either identify the mechanism (slow rate
   co-modulation entering through the baseline fit is the obvious candidate) or treat the
   connectivity as non-synaptic. Concrete test: jitter at 5, 10, 25, 50, 100 ms and find the
   scale at which edges start to disappear.
2. **Make the conditions geometrically comparable.** The aligned-ongoing control already builds
   a light-off arm with identical window count, duration and period. The `lightON` vs
   aligned-ongoing contrast is the only one currently free of the geometry confound; the raw
   `lightON` vs `ongoing` contrast is not, and should be retired.
3. **Replace the shuffle-ISI null.** At 0.0–4.6% rejection it does no work. A jitter-based null
   at the timescale identified in (1) would be a null the data can actually fail.
4. **Match spike counts, or normalise density explicitly**, before comparing conditions.
5. **Respecify the second-level model** so it is not rank deficient — likely by collapsing
   `region_pair` or dropping the three-way interaction.
6. Only then regenerate community/participation/hub-stability/circos and the manuscript figures.

Items 1 and 2 decide whether there is a paper here. Items 3–5 are prerequisites for trusting
any number in §4.

---

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
