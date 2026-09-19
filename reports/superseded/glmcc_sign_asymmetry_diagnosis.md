# Why GLMCC weights are ~all-positive in lightON and ~all-negative in ongoing

**Date:** 2026-09-18
**Status:** Root cause identified and **confirmed experimentally** (section 6).
**Bottom line:** The `lightON` epoch is a concatenation of ~108 disjoint stimulus windows with a
**~10% duty cycle**. GLMCC fits its cross-correlogram baseline assuming a continuous recording, so
for `lightON` the expected coincidence rate is underestimated roughly tenfold and essentially every
pair scores as strongly excitatory. The `ongoing` epoch is ~89% continuous and behaves normally.
The condition contrast in the manuscript is therefore confounded at the input layer, before any
analysis code runs.

---

## 1. The observation

Raw GLMCC adjacency matrices (`results/glmcc/adj_{animal}_{condition}.csv`), sign of non-zero weights:

| | negative (inhibitory) weights |
|---|---|
| lightON | **0.9%** (per animal 0.0–10.6%) |
| ongoing | **92.9%** (per animal 50.7–99.1%) |

The magnitudes differ too, which rules out a simple sign-convention flip:

| condition | min | max | median |
|---|---|---|---|
| lightON (day1) | 0.724 | 11.424 | 4.772 |
| ongoing (day1) | −1.836 | 1.106 | −0.072 |

## 2. What it is not

- **Not a file mix-up.** The pattern reproduces across two independent GLMCC run sets — the
  in-use files and the `_pycharm` variants in `GLMCC/lightON_pycharm/` and `GLMCC/ongoing_pycharm/`.
- **Not a sign convention.** `calc_PSP` in `vendor/GLMCC/glmcc.py:684` is signed by construction:
  excitatory PSP = J × 2.532, inhibitory = J × 0.612. Both signs are reachable in either condition.
  (Note this 4.1× asymmetry in the scaling constants — inhibitory weights are inherently smaller
  in magnitude, which matters for any weighted metric.)
- **Not CoNNect output mistaken for GLMCC.** The `CONNECT/` matrices occupy a different range
  (0.016–1.123, median 0.263).
- **Not spike volume.** Within each condition, the negative fraction does not track spike count or
  firing rate (lightON rho=−0.34 p=0.28; ongoing rho=+0.39 p=0.21). The split is categorical.

## 3. What it is

`evoked_ts` (the `lightON` source) contains only spikes inside stimulus windows. Pooling all units
and measuring gaps >1 s in population activity:

| condition | span | observed (gap-free) time | duty cycle |
|---|---|---|---|
| lightON | 10,247–13,510 s | 947–1,089 s | **8.0–10.4%** (median 9.7%) |
| ongoing | 16,124–19,765 s | 14,208–15,874 s | **80.3–92.4%** (median 88.8%) |

Consistent across all 12 animals with no overlap between the two ranges.

**This also corrects a second misreading.** Computing firing rate as spikes ÷ span makes `lightON`
look 8–10× quieter than `ongoing` (0.12 vs 1.0 Hz). Dividing by *observed* time instead:

| | effective rate (Hz/neuron) |
|---|---|
| lightON | 0.65 – 4.04 |
| ongoing | 0.55 – 2.84 |

The rates are comparable, and `lightON` is slightly **higher** than `ongoing` in every single animal —
which is what a stimulus-driven epoch should look like. The apparent rate suppression was an artifact
of the denominator.

### Mechanism (inferred)

GLMCC estimates J from the excess of short-lag coincidences over a slowly-varying fitted baseline.
The baseline scales with observed recording time. Feeding it a train whose file span is ~10× its
observed time under-predicts the expected coincidence count by roughly that factor, so the measured
CCG sits far above the fitted baseline for nearly every pair, and J comes out positive and inflated.
With `ongoing` at ~89% duty cycle the baseline is approximately right and J distributes around zero.

This predicts exactly what is observed: `lightON` uniformly positive with large magnitudes,
`ongoing` centred slightly below zero.

## 4. Knock-on effect on edge validation

The same structure degrades the shuffle-ISI surrogate test. `shuffle_isi` permutes inter-spike
intervals, which redistributes the long inter-window gaps throughout the train and spreads spikes
more uniformly, lowering surrogate CCG peaks. Real peaks then exceed the surrogate almost always.

Fraction of pairs whose real CCG peak exceeded **all 100** surrogates (p at the 1/100 floor):

| condition | median | range |
|---|---|---|
| lightON | 95.3% | 79.8–99.3% |
| ongoing | 84.8% | 60.8–99.4% |

Directionally as predicted, though both are near-saturated — the duty cycle aggravates a test that is
weak in both conditions. This is consistent with the separate finding that BH-FDR validation removes
only 0.5% of GLMCC edges overall.

## 5. Consequences

1. The `condition` effect on any weighted metric (node strength, local efficiency, hub score)
   is confounded: `lightON` weights are systematically inflated relative to `ongoing`.
   Today's encoding fix made both sides the same *kind* of quantity; it does not remove this bias.
2. E/I composition cannot be compared across conditions as the data stand.
3. `NEXT_STEPS.md` and the manuscript should not report a condition effect on connectivity
   strength until this is resolved.

## 6. Confirmatory experiment (RUN — mechanism confirmed)

night3 (171213, 33 units). The `ongoing` spike train was gapped with night3's own measured
`lightON` window structure (9.59 s on, 89.9 s period -> 201 windows, 10.7% duty cycle, 30,225 of
276,920 spikes retained). Within-window firing rate is preserved (0.475 vs 0.465 Hz/neuron), so the
**only** difference from the control is the observation structure. Both were run through
`vendor/glmcc-c/glmcc . 33 exp GLM`.

| run | negative weights | median | nonzero |
|---|---|---|---|
| night3 ongoing, on disk (reference) | 95.6% | -0.530 | 990 |
| **CONTROL** - ongoing, unmodified | **98.6%** | -0.557 | 1036 |
| **TEST** - ongoing, gapped to 10.7% | **14.8%** | **+0.793** | 196 |
| night3 lightON, on disk (reference) | 1.7% | +5.363 | 817 |

The control reproduces the on-disk `ongoing` matrix (98.6% vs 95.6% negative; r=+0.854 across 982
shared non-zero weights), validating the invocation. Gapping the same data then inverts the sign
distribution from 98.6% negative to 85.2% positive.

**Conclusion.** The sign inversion is caused by the gapped observation structure, not by physiology.

**Qualification.** The gapped surrogate does not fully reproduce `lightON`: median +0.79 vs +5.36,
and 196 detected connections vs 817. Real `lightON` is therefore not purely an artifact — genuine
light-evoked co-activation plausibly contributes additional, real positive coupling on top of the
baseline bias. What is established is that the sign inversion and a substantial part of the weight
inflation are artifactual, which is sufficient to invalidate the cross-condition comparison of
connection weights.

## 7. All-animal duty-cycle-matched re-run

Each animal's `ongoing` recording was gapped with **its own** measured `lightON` window structure
(`scripts/run_gapped_ongoing.py`, using `scripts/gap_spike_trains.py`) and re-run through the C
GLMCC. Gapped duty cycles landed at 10.7-11.2%, matching the native `lightON` range.

| condition | median negative | range | median weight |
|---|---|---|---|
| lightON (native ~10% duty) | **0.3%** | 0.0-31.8% | +1.772 |
| ongoing (native ~89% duty) | **95.2%** | 50.7-99.1% | -0.148 |
| ongoing_gapped (~11% duty) | **31.9%** | 14.1-42.0% | +0.720 |

All 13 animals move the same way: median drop of 63.5 points in negative fraction
(range 15.8-81.5), and the median weight changes sign from -0.148 to +0.720.

Gapping closes roughly **two-thirds of the sign gap** and **~45% of the weight gap** between the
conditions. The duty cycle is therefore the dominant driver of the asymmetry, but not the entire
explanation.

**Reading the residual.** Two contributions cannot be separated by this experiment:
genuine light-evoked coupling, and the fact that the imposed mask is a fixed period fitted to each
animal's median window rather than aligned to true light onsets — real stimulus windows deliver
synchronous drive that a phase-arbitrary mask does not reproduce. Aligning the mask to actual
stimulus times (available via `mask_stimulus` / `evoked_st` in the .mat files) would tighten this.

Outputs: `results/glmcc/adj_{animal}_ongoing_gapped.csv` (+ `.meta.json`),
`results/gapped_ongoing_summary.csv`. Comparison: `scripts/compare_sign_distributions.py`.
Caveat: night2's `lightON` reference cell (31.8% negative) is the unusable 50x50 matrix.

## 7b. Does the mask need to be stimulus-aligned? (No)

Section 7 imposed a regular period fitted to each animal's median window. A stricter control reads
the true pulse edges from `mask_stimulus` (exactly 10.00 s pulses every 90.0 s) and shifts them
+45 s into the light-off interval, giving identical window count, duration and period to `lightON`
while sampling genuine light-off data (`scripts/run_gapped_ongoing_aligned.py`).

Verified first: `ongoing` is the **complement** of `lightON` — only 3 of night3's 276,920 ongoing
spikes, and 18 of day3's 964,399, fall inside a stimulus pulse. The two conditions partition the
same session, so a phase-shifted mask selects real ongoing data with no overlap.

| condition | median negative | median weight |
|---|---|---|
| lightON | 0.3% | +1.772 |
| ongoing | 95.2% | -0.148 |
| ongoing_gapped (fitted period) | 31.9% | +0.720 |
| ongoing_aligned (real geometry) | **28.1%** | +0.736 |

Per animal, the two masks differ by a mean of -0.50 points (median |diff| 2.17, max 6.15;
Wilcoxon p=0.54) against the 94.8-point gap under explanation — about **2% of the effect**.
The effect size, not the p-value, is what carries this: n=13 is too small for a null result alone
to establish equivalence.

**Conclusion.** Mask geometry does not matter; duty cycle does. The ~28-point residual between the
aligned control and `lightON` is not an artifact of mask misalignment.

**But the residual is not yet evidence of connectivity.** The control windows sample darkness while
the `lightON` windows sample a 10 s square-wave stimulus driving many neurons at once. Stimulus-locked
rate co-modulation produces CCG peaks with no synaptic connection, and GLMCC's slowly-varying baseline
may not absorb a sharp 10 s rate transient. Distinguishing common drive from coupling needs a further
control — for example jittering spikes within each window to destroy fine timing while preserving the
stimulus-locked rate envelope, and checking how much of the residual survives.

Outputs: `results/glmcc/adj_{animal}_ongoing_aligned.csv` (+ `.meta.json`),
`results/aligned_ongoing_summary.csv`.

## 7c. Jitter control: is the residual fine-timing or slow common drive?

Each spike displaced uniformly within +/-25 ms, inside its own window, destroying structure below
that scale while leaving the 10 s stimulus envelope intact (`scripts/run_jitter_control.py`,
seed 42, spike retention >99.8%). Both lightON and the light-off control were jittered.

| condition | median negative | median weight |
|---|---|---|
| lightON | 0.3% | +1.772 |
| lightON jittered | **34.4%** | +0.658 |
| light-off control | 28.1% | +0.736 |
| light-off control jittered | **27.8%** | +0.731 |

- lightON moves +30.7 points (Wilcoxon p=0.0002).
- The light-off control moves -0.3 points (p=0.79) — unchanged.
- The control-minus-lightON gap goes from +28.1 points to -2.9: **jitter closes 90% of the
  separation**, and the small remainder reverses sign.

**Reading.** lightON's distinctive positivity depends almost entirely on structure finer than
25 ms. The light-off control has no such structure to lose. So the post-duty-cycle residual is
genuinely fine-timescale, not slow rate co-modulation.

**What this does not establish.** A +/-25 ms jitter destroys *any* sub-25 ms structure, including
precisely-locked common drive: a light pulse driving many neurons with <25 ms latency spread is
erased just as a synapse is. Fine-timescale is therefore not the same claim as synaptic.

**Tension with section 7d, stated plainly.** The lag analysis finds prominent short-latency CCG
peaks equally often in both conditions (6.0% vs 7.0%), which sits awkwardly beside a jitter result
saying lightON has fine structure the control lacks. Both can hold if the structure is weak
per-pair but systematic across pairs — below the z>=5 prominence bar yet still shifting GLMCC's
GLM fit. The two analyses use different criteria and should not be collapsed into one claim.

## 7d. CCG peak lag structure (all 13 animals)

If the lightON residual were common drive, cross-correlogram peaks should cluster at lag ~= 0 and be
symmetric. If synaptic, they should sit at short non-zero lags (1-5 ms) and be asymmetric.
`scripts/ccg_lag_structure.py` measures both, over +/-50 ms at 1 ms bins.

| | lightON | ongoing_aligned (light-off) |
|---|---|---|
| median \|peak lag\|, all pairs | 25.5 ms | 26.0 ms |
| pairs with a prominent peak (z>=5 vs flanks) | **6.0%** | **7.0%** |
| of those, \|lag\| <= 5 ms | 36.6% | 33.3% |
| of those, \|lag\| <= 1 ms | 6.2% | 5.7% |

**Most CCG "peaks" are noise.** Uniform noise over a +/-50 ms window produces a median \|peak lag\|
of exactly 25.0 ms. Both conditions measure 25.5-26.0 ms, so for the great majority of pairs the
argmax is simply the noisiest bin, not a feature.

**Neither hypothesis is supported for the residual.** Prominent peaks are no more common, no more
short-latency, and no more zero-lag in lightON than in the light-off control. Whatever separates
the two conditions after duty-cycle matching does not show up as distinctive correlogram structure.

**The sharper implication concerns edge validation.** Only 6-7% of pairs have a prominent
short-latency CCG feature, yet the surrogate + BH-FDR step passes 65-100% of pairs (section 4) and
removes 0.5% of GLMCC edges. The validation is built on the CCG peak statistic, so this is a
like-for-like comparison: it is admitting roughly an order of magnitude more edges than have any
detectable correlogram feature. That is consistent with, and independent of, the earlier finding
that shuffle-ISI surrogates are a near-saturated null.

Caveat: GLMCC itself fits a GLM rather than thresholding a CCG peak, so this does not directly
indict GLMCC's own edge calls — it indicts the CCG-peak-based validation layer applied on top.

Output: `results/ccg_lag_structure.csv` (26 rows, 13 animals x 2 conditions).

## 8. Remedies

Take one small animal (e.g. night3: 33 units, 21k lightON spikes), apply the `lightON` gap mask to
its `ongoing` spike train to impose the same ~10% duty cycle, and re-run GLMCC via the C
implementation (`vendor/glmcc-c/glmcc`). If the artificially gapped `ongoing` data flips from
~96% negative to ~all-positive, the mechanism is confirmed outright.

In order of preference:
- Give GLMCC a baseline computed from observed time rather than file span, so a gapped train is
  modelled correctly.
- Or match the two conditions on duty cycle — compare `lightON` against an `ongoing` train gapped
  the same way (exactly the TEST condition above), so the bias is common to both sides.
- Or drop the cross-condition weight comparison and report only within-condition structure.

Note the second option is directly actionable now: the gapping code is in this repo's scratch work
and the C implementation runs a 33-unit animal in minutes.

## 9. Evidence trail

- Adjacency files are byte-identical copies of legacy root `result_GLMCC_*` files; neither condition
  was produced by the current pipeline.
- Raw GLMCC `J_py_*.txt` files were not preserved, so W could not be recomputed from J directly.
- `result_GLMCC_lightON_171208.csv` is 50×50 for a 34-unit recording and is unusable
  (this is the `night2_lightON` matrix already excluded from the refit).
