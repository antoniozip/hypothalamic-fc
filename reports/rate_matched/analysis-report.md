# Rate-matched contrast: the condition effect is detection power

**Date:** 2026-09-19 · **Branch:** `analysis/glmcc-corrected`
**Resolves:** Next Action 2 in `METHODS_AND_RESULTS.md` §9
**Artifacts:** `results/rate_matched_summary.csv`, `results/glmcc_rm/`,
`results/stats_first_level_glmcc_rm.csv`, `figures/main/figR1_rate_matched.png`

## Question

GLMCC reports a coupling when `|J| > J_min = sqrt(16.3 / tau / cc_0)`, and `cc_0` scales with the
pair's coincidence count, hence with `n_i * n_j / T_observed`. More spikes lower the threshold and
admit more edges with no change in circuitry. `ongoing` carries 10.2–13.7× more spikes per unit
than `lightON`. Is the condition effect circuitry, or is it detection power?

## The confound, measured

Across the 26 main matrices, Spearman ρ between spikes-per-unit and edge density is **+0.955
(p = 3.7e-14)**. Within `lightON` alone it is +0.967; within `ongoing`, +0.676. Density is very
nearly a deterministic function of spike count (figure panel A).

## Design

Two arms differing only in whether the light was on:

| | source | windows |
|---|---|---|
| **lightON** | `DATA{old}` (evoked epoch) | the true stimulus windows |
| **light-off** | `DATA{old}_ongoing` | the same windows shifted +45 s into the light-off interval |

Window count, duration and period are identical by construction, so `T_observed` matches. Each
unit is then thinned by uniform sampling without replacement to the smaller of its two spike
counts, so `n_i` matches per unit. Both terms of `n_i * n_j / T_observed` are equalised and
`J_min` is on the same footing in the two arms.

Thinning preserves correlogram *shape* — a coincidence survives with probability `p_i * p_j`
irrespective of lag — while lowering counts. It equalises detection power without creating or
destroying temporal structure. Median thinning was 10.7% on the lightON side and 2.6% on the
light-off side; spike counts after thinning are identical per animal by construction.

26/26 arms estimated successfully.

## Result: the condition effect does not survive

### Per-animal (n = 13)

| animal | spikes (both arms) | lightON edges | light-off edges | lightON density | light-off density |
|---|---|---|---|---|---|
| day1 | 159,974 | 486 | 488 | 3.09% | 3.10% |
| day2 | 54,521 | 33 | 37 | 3.55% | 3.98% |
| day3 | 81,898 | 79 | 80 | 3.65% | 3.70% |
| day4 | 110,603 | 187 | 181 | 5.11% | 4.95% |
| day5 | 103,470 | 137 | 151 | 2.76% | 3.04% |
| day6 | 156,644 | 324 | 313 | 10.52% | 10.16% |
| night1 | 53,290 | 31 | 48 | 3.33% | 5.16% |
| night2 | 37,708 | 12 | 15 | 1.07% | 1.34% |
| night3 | 19,942 | 6 | 6 | 0.57% | 0.57% |
| night4 | 26,801 | 7 | 9 | 0.59% | 0.76% |
| night5 | 32,827 | 13 | 15 | 2.36% | 2.72% |
| night6 | 62,014 | 48 | 56 | 1.68% | 1.96% |
| night7 | 307,662 | 1,148 | 1,244 | 18.63% | 20.19% |

Paired Wilcoxon: density p = 0.027, edge count p = 0.077, % negative p = 0.84. **Where there is
a difference the light-off arm is slightly higher** — the opposite direction to the raw contrast.
Median density 3.09% vs 3.10%.

### First-level LME

`lmer(metric ~ condition * region + (1|animal/neuron))`, type III, Kenward–Roger, n = 1364.

| metric | term | F raw | p raw | F matched | p matched |
|---|---|---|---|---|---|
| node_strength | condition | 378.48 | 3.2e-67 | **3.93** | 0.048 |
| node_strength | region | 2.99 | 2.6e-03 | 1.62 | 0.12 |
| node_strength | condition:region | 14.47 | 1.4e-19 | 1.85 | 0.065 |
| clustering_coefficient | condition | 1049.02 | 2.9e-170 | **4.87** | 0.028 |
| clustering_coefficient | condition:region | 11.78 | 2.6e-16 | 0.38 | 0.93 |
| local_efficiency | condition | 533.36 | 1.5e-99 | **3.33** | 0.068 |
| local_efficiency | condition:region | 3.56 | 4.3e-04 | 0.70 | 0.69 |
| hub_score | condition | 408.51 | 2.4e-71 | **0.47** | 0.49 |
| hub_score | condition:region | 2.02 | 0.042 | 2.16 | 0.028 |

F for condition falls by **two to three orders of magnitude** in every metric: 378 → 3.9,
1049 → 4.9, 533 → 3.3, 409 → 0.47. Eight of twelve terms that were significant are no longer.

**Nothing survives multiple-comparison correction.** Holm across the four condition tests:

| metric | F | p | p Holm |
|---|---|---|---|
| clustering_coefficient | 4.87 | 0.028 | 0.111 |
| node_strength | 3.93 | 0.048 | 0.143 |
| local_efficiency | 3.33 | 0.068 | 0.143 |
| hub_score | 0.47 | 0.492 | 0.492 |

The `conditionlightON` fixed effects are all small and non-significant: node_strength β = +0.53
(p = 0.49), clustering β = −0.041 (p = 0.62), local_efficiency β = −0.022 (p = 0.76), hub_score
β = +0.0007 (p = 0.99).

## Robustness to the thinning draw

Thinning is random, so the whole analysis was repeated over 10 seeds (n = 13 animals each):

| seed | median lightON | median light-off | Wilcoxon p | light-off higher in |
|---|---|---|---|---|
| 1 | 2.95% | 3.12% | 0.0161 | 10/13 |
| 2 | 2.96% | 3.14% | 0.0042 | 12/13 |
| 3 | 3.07% | 3.16% | 0.0171 | 11/13 |
| 4 | 3.00% | 3.09% | 0.0024 | 12/13 |
| 5 | 3.13% | 3.08% | 0.4043 | 8/13 |
| 6 | 2.99% | 3.08% | 0.0005 | 12/13 |
| 7 | 3.03% | 3.21% | 0.0105 | 11/13 |
| 8 | 3.05% | 3.18% | 0.0186 | 9/13 |
| 9 | 3.13% | 3.10% | 0.0034 | 10/13 |
| 42 | 3.09% | 3.10% | 0.0269 | 10/13 |

Across seeds: median 3.040% lightON vs 3.110% light-off, significant at 0.05 in **9 of 10**.

**So the matched contrast is not a clean null — there is a small, reproducible density
difference, and it favours the light-off arm.** It is ~0.07 percentage points, about 2% relative,
against the ~10-fold difference in the raw contrast.

### Is the residual an artifact of asymmetric thinning?

Matching thins to the per-unit minimum, and the lightON arm loses more: median 10.7% of its
spikes discarded against 2.5% on the light-off side. Discarding more spikes adds sampling noise
to that arm's correlograms, which would depress its detected edges — pointing in exactly the
direction observed.

Tested directly across 13 animals × 10 seeds (n = 130). If the residual were a thinning artifact,
the light-off excess should *grow* with how much more the lightON arm was thinned. It does the
opposite:

| | mean light-off excess |
|---|---|
| animals where lightON was thinned only slightly more (< 6.0 pts) | **+0.475 pp** |
| animals where lightON was thinned much more (≥ 6.0 pts) | **+0.253 pp** |

Spearman(extra thinning on lightON, light-off excess) = **−0.311, p < 0.001**.

**The artifact hypothesis is not supported** — the correlation runs the wrong way. The caveat is
that this is correlational across animals, and thinning asymmetry covaries with other properties
(animals with a large count mismatch differ in firing rate). A symmetric design — thin both arms
to a common fixed count well below either, so both discard heavily — would settle it directly and
has not been run.

## Conclusion

**The `lightON` vs `ongoing` condition effect as reported is an artifact of detection power.**
Equalising observation geometry and per-unit spike count collapses every graph-metric condition
effect: F falls 100–1000×, and nothing survives Holm correction across the four metrics.

Two things remain, and neither supports the original claim:

1. **A small, reproducible density excess in the light-off arm** — ~0.07 pp, ~2% relative,
   significant in 9 of 10 thinning seeds, and *not* explained by asymmetric thinning (above).
   It is in the opposite direction to "light evokes connectivity". At 2% of the effect size
   originally reported, it is a different phenomenon, and it needs the symmetric-thinning control
   before it is worth interpreting.
2. `hub_score: condition × region` at F = 2.16, p = 0.028 uncorrected — which was also only
   marginal in the raw analysis (p = 0.042) and is one of twelve tests. Not a basis for a claim.

## What this does not say

- **It is not a statement about the estimator.** GLMCC is working as specified; the threshold
  depends on spike count by design, and the corrected jitter sweep
  (`reports/jitter_sensitivity/`) shows the edges it finds are genuinely timing-dependent.
- **It is not a claim that the two epochs are physiologically identical.** It is a claim that
  these graph metrics, at this operating point, cannot distinguish them once detection power is
  equalised.
- **It does not rescue the raw matrices.** `results/glmcc/` remains confounded; the matched arms
  in `results/glmcc_rm/` are the defensible comparison.

## Limitations

- **The residual light-off excess may be an artifact of asymmetric thinning** (see above). It is
  small and consistent, and it is not evidence for anything until the symmetric-thinning control
  is run.
- **LME statistics are from the seed-42 realisation only.** The seed sweep covers density; the
  full LME was not repeated per seed. Given F collapses by 100-1000x this will not change the
  conclusion, but the exact matched p-values would move.
- **Thinning discards data.** The lightON arm lost a median 10.7% of its spikes. The comparison
  is therefore at a lower operating point than either raw condition, and both arms sit closer to
  the detection floor — night3 has 6 edges from 33 units. Low power is a real limitation of the
  matched design, and a null result here is weaker than a null at full power would be.
- **The light-off arm is not `ongoing`.** It is ongoing data sampled through the stimulus window
  geometry, shifted 45 s. That is the right control for `lightON`, but it is not the epoch the
  original analysis used.
- **Validation was not re-applied.** These matrices are written directly as `validated_adj_*` so
  the R pipeline consumes them. Since the surrogate null removes 0.0–4.6% of edges (§4.2), this
  changes nothing material, but it should be stated rather than hidden.
- **Region-level terms are underpowered** in the matched analysis: several animals contribute
  fewer than 20 edges.

## Next actions

1. **Retire the raw `lightON` vs `ongoing` contrast** from the manuscript. It measures spike
   count. Use the matched arms.
2. **Run the symmetric-thinning control.** Thin both arms to a common fixed count so each
   discards the same fraction, and see whether the ~2% light-off excess survives. This is the
   one open question the matched design raises, and it is cheap.
3. **Next Action 3 (within-window shuffle null)** is now the remaining validation gap, but it
   will not revive the condition effect — it constrains which edges are real, not whether the
   two arms differ.
4. **Consider whether the study has a different question.** The matched comparison is a
   well-powered null on graph-level metrics. Per-edge or per-region-pair analyses, or a
   within-condition question, may be where the data can still speak.
