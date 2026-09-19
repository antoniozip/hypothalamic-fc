# Statistics appendix — jitter sensitivity sweep

**Date:** 2026-09-19 · source: `results/jitter_sensitivity.csv` (91 rows, 13 animals × 7 levels)

## Design

| | |
|---|---|
| Unit of analysis | animal (n = 13) |
| Arm | aligned-ongoing: `ongoing` gapped to the real pulse geometry, phase-shifted 45 s into light-off |
| Baseline | the same arm with no displacement, re-estimated per animal |
| Primary metric | edge survival = \|E_base ∩ E_jitter\| / \|E_base\| (higher = less effect) |
| Secondary | Jaccard(E_base, E_jitter); Pearson r of weights over surviving edges; sign agreement |
| Levels | 0, 5, 10, 25, 50, 100 ms uniform ±, plus full within-window shuffle |
| Seed | 42 (one jitter draw per animal × level) |
| Estimator | `vendor/glmcc-c/glmcc exp GLM`, parity-verified (`tests/test_glmcc_c_parity.py`) |

Geometry is identical across every level within an animal, so displacement is the only variable.

## Validity checks

- **Identity check.** Level 0 gives survival 100.00% (SD 0.00), Jaccard 1.0000, weight r 1.0000,
  41,842 → 41,842 edges in all 13 animals. The comparison machinery introduces no drift.
- **Spike retention.** Displaced spikes leaving their window are dropped, not clipped. Retention
  falls monotonically 100.000 → 99.495% from 0 to 100 ms. The shuffle level has 100% retention by
  construction and shows the *largest* effect, so spike loss is not the mechanism.
- **Metric direction.** Higher survival = smaller effect of jitter. A synaptic edge set should
  show survival falling toward chance as displacement exceeds the synaptic lag.

## Descriptive statistics

| level | survival mean ± SD | median [IQR] | min | Jaccard | weight r | sign agree | retention |
|---|---|---|---|---|---|---|---|
| 0 | 100.00 ± 0.00 | 100.00 [100.00, 100.00] | 100.00 | 1.0000 | 1.0000 | 100.00% | 100.000% |
| 5 | 99.98 ± 0.03 | 100.00 [100.00, 100.00] | 99.90 | 0.9998 | 1.0000 | 100.00% | 99.971% |
| 10 | 99.97 ± 0.04 | 100.00 [99.96, 100.00] | 99.90 | 0.9997 | 0.9999 | 100.00% | 99.948% |
| 25 | 99.95 ± 0.05 | 99.96 [99.91, 100.00] | 99.88 | 0.9993 | 0.9998 | 100.00% | 99.868% |
| 50 | 99.93 ± 0.09 | 99.98 [99.90, 100.00] | 99.68 | 0.9993 | 0.9996 | 100.00% | 99.746% |
| 100 | 99.92 ± 0.09 | 99.98 [99.85, 100.00] | 99.73 | 0.9988 | 0.9993 | 100.00% | 99.495% |
| shuffle | 99.30 ± 1.14 | 99.58 [99.40, 99.80] | 95.63 | 0.9919 | 0.9602 | 99.98% | 100.000% |

## Per-animal edge survival (%)

| animal | 0 | 5 | 10 | 25 | 50 | 100 | shuffle |
|---|---|---|---|---|---|---|---|
| day1 | 100.00 | 100.00 | 100.00 | 100.00 | 100.00 | 99.99 | 99.98 |
| day2 | 100.00 | 100.00 | 100.00 | 100.00 | 100.00 | 100.00 | 99.54 |
| day3 | 100.00 | 100.00 | 100.00 | 99.93 | 99.87 | 99.73 | 99.40 |
| day4 | 100.00 | 99.97 | 99.97 | 99.97 | 100.00 | 99.85 | 99.58 |
| day5 | 100.00 | 100.00 | 100.00 | 99.98 | 99.98 | 99.98 | 99.20 |
| day6 | 100.00 | 99.93 | 99.93 | 99.96 | 99.68 | 99.79 | 99.89 |
| night1 | 100.00 | 100.00 | 100.00 | 100.00 | 100.00 | 100.00 | 100.00 |
| night2 | 100.00 | 100.00 | 100.00 | 99.90 | 100.00 | 100.00 | 99.80 |
| night3 | 100.00 | 99.90 | 99.90 | 99.90 | 99.90 | 100.00 | 95.63 |
| night4 | 100.00 | 100.00 | 100.00 | 99.91 | 99.91 | 99.91 | 99.16 |
| night5 | 100.00 | 100.00 | 100.00 | 100.00 | 100.00 | 100.00 | 99.62 |
| night6 | 100.00 | 100.00 | 99.96 | 99.88 | 99.84 | 99.88 | 99.43 |
| night7 | 100.00 | 100.00 | 99.91 | 99.91 | 99.91 | 99.84 | 99.69 |
## Inferential tests

One-sample Wilcoxon signed-rank against a 100%-survival null (the null that jitter has no
effect). n = 13 paired observations per level. Holm correction across the six non-zero levels.
Non-parametric because survival is bounded at 100 and strongly left-skewed; Shapiro–Wilk rejects
normality at every level.

| level | median survival | W | p raw | p Holm | median loss (pts) |
|---|---|---|---|---|---|
| 5 ms | 100.00 | 0.0 | 0.250 | 0.250 | 0.00 |
| 10 ms | 100.00 | 0.0 | 0.063 | 0.125 | 0.00 |
| 25 ms | 99.96 | 0.0 | 0.0039 | **0.020** | 0.04 |
| 50 ms | 99.98 | 0.0 | 0.016 | **0.047** | 0.02 |
| 100 ms | 99.98 | 0.0 | 0.0078 | **0.031** | 0.02 |
| shuffle | 99.58 | 0.0 | 0.00049 | **0.0029** | 0.42 |

### How to read these rejections

They are real but they are not the finding. Between-animal variance is on the order of 0.1
percentage points, so the test has power to detect a median loss of 0.02 points. **Effect size is
the relevant quantity, and the largest across all levels is 0.42 percentage points.** For scale, an
edge set resting on 1–5 ms structure would lose on the order of 94 points by 25 ms.

Reporting p-values alone here would invert the conclusion. The correct statement is: the effect of
destroying all sub-window timing structure is statistically detectable and practically nil.

## Mechanism test

Mann–Whitney U on the spike-count product $n_i n_j$, pairs with an edge vs pairs without, in the
unjittered baseline. Three animals spanning the size range.

| animal | n units | edges / ordered pairs | median $n_i n_j$ present | median absent | ratio | U p |
|---|---|---|---|---|---|---|
| night5 | 24 | 530 / 552 | 1,193,248 | 302,022 | 4.0× | 2.07e-08 |
| day3 | 47 | 1,508 / 2,162 | 1,159,561 | 1,674 | 693× | 8.51e-297 |
| night6 | 54 | 2,452 / 2,862 | 234,234 | 15,351 | 15.3× | 4.16e-227 |

Spearman correlation between $n_i n_j$ and |W| across all ordered pairs is −0.706 (night5),
+0.301 (day3), −0.222 (night6). The sign is inconsistent, so **no claim is made about weight
magnitude**; the consistent and extreme result is on edge *presence*.

## Blockers and limitations

- **Single arm.** Only aligned-ongoing was swept. The mechanism is a property of GLMCC's
  `J_min ∝ 1/sqrt(cc₀)` threshold rather than of this arm, but `lightON` was not tested.
- **Single jitter draw per cell.** Seed 42 throughout. Reported variability is between animals,
  not between jitter realisations. Repeating with several seeds would tighten the SD but cannot
  change a result this far from the alternative.
- **Mechanism panel is 3 of 13 animals.** Direction identical, p-values extreme, but not exhaustive.
- **No test of the alternative.** This analysis shows the edges do not depend on timing and do
  depend on spike counts. It does not establish that spike count is a *sufficient* generative
  model — that would need a rate-matched surrogate that reproduces the observed edge set.
- The schematic synaptic-decay curve in figure panel A is illustrative, not fitted or measured.
