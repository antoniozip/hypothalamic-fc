# GLMCC positive control: the estimator does not recover known couplings

**Date:** 2026-09-19 · `scripts/glmcc_positive_control.py` · `results/glmcc_positive_control.csv`

## Setup

30 neurons, independent Poisson trains, with **40 known connections injected** (20 excitatory,
20 inhibitory — 4.6% of ordered pairs). Excitatory connections copy a presynaptic spike to the
postsynaptic train at 2 +/- 0.5 ms latency with p=0.30; inhibitory ones delete postsynaptic spikes
in a 3 ms window with the same probability. That is a strong, clearly detectable synapse.

Scenarios match the real recording regimes, so the answer is specific to this dataset.

## Result

| scenario | rate | duration | spk/neuron | density | **recall** | **precision** | **sign correct** | detected negative |
|---|---|---|---|---|---|---|---|---|
| generous | 2.00 Hz | 15,000 s | 36,448 | 9.20% | **15.0%** | **7.5%** | **16.7%** | 100.0% |
| ongoing-like | 1.40 Hz | 15,000 s | 25,553 | 7.82% | 2.5% | 1.5% | 0.0% | 98.5% |
| lightON-like | 0.65 Hz | 1,000 s | 788 | 4.37% | 2.5% | 2.6% | 0.0% | 100.0% |

Both binary modes fail on the generous scenario; `sim` is not better than `exp`:

| mode | density | recall | precision | sign correct | detected negative |
|---|---|---|---|---|---|
| exp | 9.20% | 15.0% | 7.5% | 16.7% | 100.0% |
| sim | 6.09% | 7.5% | 5.7% | 33.3% | 100.0% |

## Reading

Even in the most favourable regime — 36,000 spikes per neuron, far more than any real recording
here — the pipeline recovers **15% of true connections**, and **92.5% of what it reports is false**.

The decisive observation is the last column. Ground truth contains 20 excitatory connections, and
the estimator returns **100% negative weights**. That is not an underpowered detector; an
underpowered detector misses connections, it does not invert their sign. Something in the fit or in
the sign assignment is systematically driving J negative.

This reproduces exactly the signature seen on the real data (99.8-100% negative across all 26
datasets), and explains the other anomaly too: if reported edges are mostly false positives, their
count should track statistical power rather than circuitry — which is what the real data shows
(density vs spike count, rho=+0.938).

## What this establishes, and what it does not

**Establishes:** the GLMCC pipeline as vendored and invoked in this repository does not measure
connectivity. Network metrics computed on its output — density, hub scores, rich-club,
participation, community structure — are not interpretable, in either condition.

**Does not establish** that GLMCC the *method* is at fault. Three alternatives remain open:

1. **A porting bug in the C implementation.** The immediate next test is running the vendored
   Python `Est_Data.py` on this same fixture. If Python recovers the couplings and C does not, the
   defect is in the port and is fixable. This is the highest-value remaining check.
2. **Parameter mismatch.** Kobayashi et al. tuned for cortical data; the `Jmin`/`scale`/`z_a`
   thresholds may be wrong for hypothalamic rates.
3. **Injection-model mismatch.** The additive-copy model is standard, but if GLMCC assumes a
   different generative form the fixture could understate its performance.

## Consequence

Until a positive control passes, no connectivity result from this pipeline should be reported.
That includes the corrected LME in `reports/corrected_lme_results.md`: its inputs come from this
estimator, so the condition effects it measures cannot currently be attributed to connectivity.
