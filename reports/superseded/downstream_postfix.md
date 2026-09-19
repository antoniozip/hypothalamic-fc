# Graph metrics and tier 1/2 analyses on corrected connectivity

**Date:** 2026-09-19 · inputs: `results/glmcc/` rebuilt from the fixed C port

## What was run

1. Corrected matrices installed into `results/glmcc/adj_*.csv` (previous contents backed up).
2. `scripts/revalidate_all.py` — 26/26 ok, fresh CCG + 100 shuffle-ISI surrogates per dataset.
3. `scripts/rebuild_downstream.sh` — graph metrics **26 ok, 0 skipped, 0 failed**; first-level LME.
4. `src/r/stats_second_level.R` — reran after installing `pbkrtest`.
5. `src/r/tier1_analyses.R`, `src/r/tier2_analyses.R`.

Two prerequisites had to be fixed to get here; both are recorded at the end.

## First-level LME — omnibus ANOVA

| metric | term | F before | p before | F after | p after |
|---|---|---|---|---|---|
| node_strength | condition | 1253.66 | 5.0e-193 | 378.48 | 3.2e-67 |
| node_strength | region | 0.20 | 0.99 | 2.99 | 2.6e-03 |
| node_strength | condition:region | 37.59 | 5.6e-54 | 14.47 | 1.4e-19 |
| clustering_coefficient | condition | 3.43 | 0.065 | **1049.02** | **2.9e-170** |
| clustering_coefficient | region | 2.28 | 0.021 | 1.13 | 0.34 |
| clustering_coefficient | condition:region | 6.99 | 8.6e-09 | 11.78 | 2.6e-16 |
| local_efficiency | condition | 3511.22 | <1e-300 | 533.36 | 1.5e-99 |
| local_efficiency | region | 0.16 | 1.00 | 2.60 | 8.0e-03 |
| local_efficiency | condition:region | 50.33 | 4.2e-71 | 3.56 | 4.3e-04 |
| hub_score | condition | 53.45 | 4.6e-13 | 408.51 | 2.4e-71 |
| hub_score | region | 0.76 | 0.64 | 5.12 | 3.3e-06 |
| hub_score | condition:region | 2.92 | 3.1e-03 | 2.02 | 0.042 |

Every omnibus term stays significant, but the magnitudes move by one to two orders of
magnitude in both directions. `clustering_coefficient ~ condition` goes from p = 0.065 to
p = 2.9e-170. `region` was non-significant for three of four metrics before and is now
significant for three of four.

### A sign reversal in a headline effect

`local_efficiency ~ conditionlightON` is significant before and after, in **opposite
directions**:

| | estimate | p |
|---|---|---|
| before | **+1.341** | 9.5e-10 |
| after | **−0.205** | 9.8e-03 |

Eight of the local_efficiency terms that were significant before are not any more; eight
clustering_coefficient and node_strength terms that were not significant have become so.
`hub_score ~ conditionlightON` goes from p = 0.80 to p = 2.8e-03.

## Second-level LME

Omnibus, Kenward-Roger: condition F = 314.94, p < 2e-16; day_night F = 0.0016, p = 0.97;
condition:day_night F = 0.073, p = 0.79. All 156 coefficient estimates changed. The four
main terms are non-significant before and after (|t| < 2 throughout).

The model remains **rank deficient — 76 columns dropped** for missing `region_pair × day_night`
cells, so all emmeans contrasts come back `nonEst`. That is a pre-existing design problem,
unchanged by this work, and it means the second level currently supports no pairwise claim.

## The validation step removes almost nothing

Comparing raw against validated matrices for the 13 canonical animals:

| condition | raw edges | validated edges | removed |
|---|---|---|---|
| lightON | 3,179 | 3,179 | 0 (0.0%) |
| ongoing | 19,993 | 19,067 | 926 (4.6%) |

Median edge-survival rate across the 26 runs is 0.948. The CCG + shuffle-ISI + BH-FDR layer
passes essentially every pair, so "validated" connectivity is very nearly the raw GLMCC output.
This restates, on corrected matrices, the concern already recorded in section 4 of
`reports/glmcc_sign_asymmetry_diagnosis.md`: the surrogate null is close to saturated.

## Tier 1

| | before | after |
|---|---|---|
| **Edge consistency** (mean) | 0.636 | **0.075** |
| Rich club, rows | 1,145 | 435 |
| Rich club, φ_norm (mean) | 0.989 | 1.008 |
| Centrality, mean strength | 326.8 | 5.85 |
| Centrality, mean eigenvector | 0.472 | 0.215 |

**Edge consistency falls 8.5-fold**, from 0.636 to 0.075. Rich-club structure is now absent in
9 of 13 animals (0 significant k-levels); it survives in night7 (57 levels), day4 (11), day5
(11) and day6 (2).

Top consistent region pairs are now PH↔PVH 0.183, VMH↔PH 0.167, PH↔VMH 0.151, DMH↔PH 0.144,
MC↔DMH 0.136 — all far below the old mean.

Within- vs between-region decomposition: group mean 43.4% within / 56.6% between, but the
spread across animals is enormous (day1 97.9% within; day2 21.1% within).

Centrality convergence (mean Spearman ρ across strength/betweenness/eigenvector/pagerank)
ranges 0.298 (night3) to 0.939 (day6).

## Tier 2

**Surrogate null models.** Mean Z_modularity falls from 3.51 to 1.28 and Z_clustering rises
from −2.78 to +0.43. Modularity now exceeds its degree-preserving null in **3 of 13** animals,
clustering in **2 of 13**. Z_pathlength is undefined everywhere (disconnected graphs).

**Small-world propensity is no longer computable.** σ and ω require a single connected
component. Seven animals had them before; **none does now** — every lightON graph is
disconnected:

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

The corrected lightON connectivity is roughly ten times sparser than before (3,179 edges
across all animals, against 32,633). night3 has three edges.

## How to read this

The condition effect does not go away — it gets *stronger* for clustering and hub_score. But
`reports/glmcc_regeneration_postfix.md` shows the sign and density of these matrices are set
largely by recording geometry (continuous ongoing ≈ 98% negative; the same spikes gapped to the
lightON duty cycle ≈ 100% positive), and that essentially every aligned-ongoing edge survives
±25 ms jitter. A large, highly significant lightON-vs-ongoing contrast is exactly what that
artifact predicts. These numbers are now computed correctly; whether the contrast they describe
is biological is a separate question, and the jitter control argues against it.

Edge consistency at 0.075, rich club gone in 9/13, small-world undefined and modularity nulls
passing in 3/13 all point the same way: after correcting the estimator there is much less
network structure than the previous results reported.

## Two things that had to be fixed to run this

**1. `pbkrtest` was missing**, so `src/r/stats_second_level.R` died at
`anova(model, type="III", ddf="Kenward-Roger")` before writing anything — the second-level CSV
on disk was stale. Installed via CRAN archive (`Deriv 4.1.3` is not available for R 4.3.3 from
the current CRAN index, which blocks `doBy` and then `pbkrtest`). The script itself was not
changed, so it still uses its intended Kenward-Roger method. Note `stats_first_level.R` already
had a KR → type-II fallback; the second level has none.

**2. `src/r/tier2_analyses.R` crashed on an empty result set.** It skips animals whose graph is
disconnected, then computes `mean(sw_results$sigma)` and pivots `c(sigma, omega)` unguarded. With
no animal connected, `sw_results` has zero rows and no such columns. Added a guard so it reports
"small-world indices are undefined" and skips `figT2_small_world.png` instead of halting. No
change to any statistic — the only reason this surfaced now is that the corrected connectivity is
sparse enough to disconnect every graph.

## Not re-run

`community_detection.R`, `compute_participation.R`, `hub_stability.R`, `make_figures.R`,
`plot_circos.R`, `plot_second_level.R`, and the tier3 Python analyses. `figures/main/figT2_small_world.png`
is stale (10:18) because there was nothing to plot.
