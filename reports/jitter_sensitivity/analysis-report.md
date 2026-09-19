# Jitter sensitivity: at what timescale does the GLMCC edge set break?

**Date:** 2026-09-19 · **Branch:** `analysis/glmcc-corrected`
**Resolves:** Next Action 1 in `METHODS_AND_RESULTS.md` §9
**Artifacts:** `results/jitter_sensitivity.csv`, `figures/main/figJ1_jitter_sensitivity.png`,
`results/glmcc/sensitivity/`

## Analysis question

The ±25 ms jitter control left 100.0% of aligned-ongoing edges standing. Either the edge set is
not synaptic, or the control was not doing what it claimed. This sweep varies the displacement
and asks: **at what scale do edges start to disappear?**

## Answer

**They never do.** There is no timescale between 5 ms and the full 10-second observation window
at which the edge set meaningfully changes.

| displacement | edges surviving (mean ± SD, n=13) | weight r vs baseline | spike retention |
|---|---|---|---|
| 0 ms (identity) | 100.00 ± 0.00 % | 1.0000 | 100.000% |
| 5 ms | 99.98 ± 0.03 % | 1.0000 | 99.971% |
| 10 ms | 99.97 ± 0.04 % | 0.9999 | 99.948% |
| 25 ms | 99.95 ± 0.05 % | 0.9998 | 99.868% |
| 50 ms | 99.93 ± 0.09 % | 0.9996 | 99.746% |
| 100 ms | 99.92 ± 0.09 % | 0.9993 | 99.495% |
| **full within-window shuffle** | **99.30 ± 1.14 %** | **0.9602** | 100.000% |

91/91 runs completed. The 0 ms level reproduces the baseline exactly (survival 100.00%,
r = 1.0000, 41,842 → 41,842 edges), so the comparison machinery is sound and every other number
here is measured against a correct reference.

**The decisive row is the last one.** A full within-window shuffle redraws every spike time
uniformly inside the 10 s window it came from. It preserves each cell's spike count per window
and the window geometry, and destroys every temporal relationship finer than a window — there is
no correlogram structure left to find. **99.30% of edges survive it, and the weights are still
correlated with the originals at r = 0.96.**

An edge set built on 1–5 ms synaptic structure cannot behave this way. By 25 ms a real
monosynaptic peak is smeared across half the ±50 ms correlogram window; essentially every such
edge should be gone.

## What the edge set does depend on

If not timing, then what? Edge presence is predicted by the pair's spike-count product:

| animal | edges | median $n_i n_j$, edge present | median $n_i n_j$, absent | ratio | Mann–Whitney p |
|---|---|---|---|---|---|
| night5 | 530 / 552 | 1,193,248 | 302,022 | 4.0× | 2.1e-08 |
| day3 | 1,508 / 2,162 | 1,159,561 | 1,674 | **693×** | 8.5e-297 |
| night6 | 2,452 / 2,862 | 234,234 | 15,351 | 15× | 4.2e-227 |

This is GLMCC's detection threshold behaving exactly as specified rather than malfunctioning.
The threshold is `J_min = sqrt(16.3 / τ / cc₀)`, where `cc₀` is the fitted correlogram height near
the synaptic delay. `cc₀` grows with the pair's coincidence count, which for independent trains
is proportional to $n_i n_j$. More spikes ⇒ smaller threshold ⇒ edge. When the correlogram carries
no peak — as after a full shuffle — the threshold still falls with spike count, and `J` still
absorbs whatever mismatch remains between the stiff smoothed baseline (β = 4000 over 100 bins) and
the observed level. The result is an edge call driven by counts, not coupling.

This also explains the density pattern already noted in `METHODS_AND_RESULTS.md` §4.1: density
ranges 0.28% to 62.67% and tracks spike count across datasets.

## Statistical note

Wilcoxon signed-rank tests against a 100%-survival null (n = 13, Holm-corrected across the six
non-zero levels) reject at 25 ms (p = 0.020), 50 ms (p = 0.047), 100 ms (p = 0.031) and shuffle
(p = 0.0029); 5 ms (p = 0.25) and 10 ms (p = 0.13) do not reject.

**These rejections are not the finding.** The between-animal variance is so small that a median
loss of 0.02–0.42 percentage points is detectable. The effect size is what matters: the worst
case across every level is a **0.42 point** median loss, against the ~94 point loss a synaptic
edge set would show by 25 ms. Statistically detectable, materially nil.

## What this changes

1. **Next Action 1 is resolved, negatively.** The question was "find the scale at which edges
   start to disappear". There is none. The prior ±25 ms result was not an artifact of a badly
   chosen displacement; the edge set simply does not encode spike timing.

2. **The synaptic interpretation is not available for these matrices.** Whatever
   `adj_*_ongoing_aligned.csv` and, by the same mechanism, the main `adj_*_{lightON,ongoing}.csv`
   measure, it is not millisecond-scale coupling. Network statistics computed on them — the
   condition effects in §4.3, rich club, modularity, edge consistency — describe a graph whose
   edges are set by firing rates and recording geometry.

3. **The condition effect now has a sufficient non-biological explanation.** `ongoing` carries
   10–30× more spikes than `lightON`. If edge presence follows spike count, a large and highly
   significant `lightON` vs `ongoing` contrast follows from that alone, with no circuitry
   involved. This sits alongside the geometry result (§5.3) — sign flips with duty cycle — and
   the two together account for the contrast without appeal to connectivity.

4. **The shuffle test should replace the shuffle-ISI surrogate null.** The existing validation
   rejects 0.0% of lightON and 4.6% of ongoing edges. A within-window-shuffle null applied
   per-pair would reject ~99% of them, and is a null the data can actually fail.

## Limitations

- One arm only (aligned-ongoing). It was chosen because its geometry is fixed across levels, so
  displacement is the only variable. The mechanism is a property of GLMCC's threshold, not of
  this arm, so it should generalise — but `lightON` was not swept and that has not been shown.
- The mechanism panel covers 3 of 13 animals (night5, day3, night6), chosen for size. The
  direction is identical in all three and the p-values are extreme, but it is not the full set.
- One seed per level (42). Between-animal spread (n = 13) is the variability reported; jitter-draw
  variability within an animal is not characterised.
- Spike retention falls to 99.5% at 100 ms as displaced spikes leave their windows. This is far
  too small to account for the survival numbers, and the shuffle level has 100% retention with the
  largest effect, so loss of spikes is not driving the result.
- The red dashed curve in figure panel A is a **schematic** expectation for a synaptic edge set,
  not a measurement. It is drawn to give the flat blue curve a scale.

## Next actions

1. **Retire the synaptic framing** for the current matrices in `METHODS_AND_RESULTS.md`. Done in
   this commit.
2. **Replace the surrogate null** with a per-pair within-window shuffle at q = 0.05, and re-run
   validation. This is the single highest-value next step: it is the test that the present edge
   set fails, so it will produce an edge set that at least carries timing information.
3. **Match spike counts between conditions** before any condition contrast, or abandon the raw
   contrast. Sub-sampling `ongoing` to `lightON` counts is the cleanest version.
4. If (2) leaves too few edges to analyse — likely, given 99% would be rejected — then GLMCC at
   this operating point does not support a connectivity claim on this dataset, and that is the
   result.
