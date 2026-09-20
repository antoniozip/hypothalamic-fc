# Re-deriving E/I from waveforms using a published method

**Date:** 2026-09-20 · **Branch:** `analysis/glmcc-corrected`
**Source method:** Mahallati et al., *eLife* 2026 (spikeMAP), https://elifesciences.org/articles/106557
**Script:** `scripts/ei_classify_waveforms.py` · **Output:** `data/processed/waveform_features.csv`,
`results/ei/classification_compare.csv`

## The method, as published

Two features of the mean spike waveform:

> "the full width (FW) of spikes at half-maximum amplitude and (2) the peak-to-peak (PP)
> duration"

classified without fixed thresholds:

> "K-means clustering was performed on values of FW and PP obtained across all simultaneously
> recorded somatic electrodes, where the optimal number of clusters was chosen with a
> Calinski-Harabasz criterion"

with inhibitory cells being "a single cluster with rapid time constants for both FW and PP".
In mouse prefrontal cortex at 18 kHz they classified **5.25%** of 8,189 units as inhibitory.

Waveforms here are `Units(i).wf`, 81 samples at 20 kHz, in `DATA{old}.mat`. 661 of 682 units
yield usable features.

## A real bug in the original classification

`tim_stuff.m:702-703` computes its half-width as

```matlab
[pks,locs,ws{i},proms] = findpeaks(Units(i).wf);
half_w(i) = ws{i}(end)/20000;
```

`findpeaks` returns **positive** peaks. An extracellular spike's dominant deflection is the
negative trough, so `ws{i}(end)` is the half-prominence width of the *last positive bump* in the
waveform — a property of the trailing repolarisation, not the spike's half-width. Combined with
hard thresholds (FW < 0.6 ms AND PP < 0.65 ms), that is the likely origin of the implausible
44.6% inhibitory fraction the shipped labels give.

Here FW is measured on the trough, at half its depth below a pre-spike baseline, linearly
interpolated on both flanks.

## Result: the method does not produce a usable classification here

Three variants were run. **None gives a credible interneuron fraction.**

| variant | inhibitory | note |
|---|---|---|
| pooled k-means over all 661 units | **39.0%** | k=2 by Calinski-Harabasz |
| pooled, each recording centred on its own median | **81.4%** | splits off the slow tail; boundary in the wrong place |
| **per recording, as the paper specifies** | **27.1%** | k varies 2–6; range 3.2%–59.6% across recordings |

Against a reference of 5.25% in the source paper, and ~10–20% in most cortical datasets.

## Why it fails: recording-level variation swamps the cell-type signal

Median PP per recording spans **0.400 to 1.050 ms — a 2.6× range between animals.**

That is larger than the E/I separation it is supposed to detect, and it dominates any pooled
clustering. With pooled k-means the per-animal inhibitory fraction is almost perfectly predicted
by where that recording sits:

**Spearman(recording median PP, its % called inhibitory) = −0.978, p = 8.0e-09, n = 13.**

That is a rig or preprocessing effect being read out as biology.

Per-recording clustering removes the offset by construction, but then the classes are no longer
comparable: the "inhibitory" centroid ranges from PP 0.222 ms (night7) to 0.891 ms (night3), a
**4-fold** spread. A unit called inhibitory in night3 would be among the slowest cells in
night7. There is no common definition of "narrow-spiking" across these recordings.

## Is there any real two-class structure?

Partly, and it does not rescue the classification.

| test | result |
|---|---|
| Hartigan dip, FW pooled | D = 0.0114, p = 0.87 — **unimodal** |
| Hartigan dip, PP pooled | D = 0.0408, p < 1e-4 — bimodal |
| Hartigan dip, PP after centring each recording | D = 0.0590, p < 1e-4 — still bimodal |
| Hartigan dip, PP within each recording | bimodal in **6 of 13** |

So PP carries genuine sub-structure that survives removing the recording offset, while FW does
not. But bimodality in one feature, present in under half the recordings, with a 4× shift in
where the boundary falls, is not enough to assign cell classes.

## Agreement with the shipped labels

On the 379 units the original classification covers:

| | new E | new I |
|---|---|---|
| **old E** | 172 | 36 |
| **old I** | 77 | 94 |

Raw agreement **70.2%**, Cohen's κ = **+0.385** (fair). Better than chance (χ² p < 0.001), so the
two are picking up a common signal — but they disagree on 113 of 379 units, and neither produces
a plausible interneuron fraction.

## Conclusion

**The published method was applied faithfully and does not rescue the classification.** The
limitation is in the data, not the criterion: between-recording waveform variation exceeds the
within-recording cell-type separation, so no single threshold or clustering rule is valid across
these 13 recordings.

This supersedes the recommendation in `reports/ei_differential/analysis-report.md` to
"re-derive the waveform classification" — it has been done, and the answer is that the data will
not support it as it stands.

The E/I nulls in that report therefore stand as *untestable* rather than negative: neither the
original labels nor a literature-standard re-derivation gives a classification worth testing.

## What would be needed

1. **Find out why recordings differ 2.6× in median PP.** Different probe, filter settings,
   sampling alignment, or spike-sorting template width would all do it. If it is a filtering or
   alignment difference it can be corrected and the classification retried.
2. **Check whether waveforms were re-aligned or resampled** between sessions — the 81-sample
   window at 20 kHz should be identical across recordings if they came off one pipeline.
3. **Use an independent label** if one exists: optotagging, or the GLMCC couplings themselves
   (circular for connectivity questions, but valid for validating the waveform criterion).
4. Failing all that, drop the E/I dichotomy from this dataset rather than reporting nulls from a
   classification that does not classify.

## Limitations

- One published method, one parameterisation of the features. Other definitions of FW (e.g. at
  half-*prominence*, or on the peak rather than the trough) would shift absolute values, but
  cannot remove a 2.6× between-recording offset.
- Per-recording k-means on 24–71 units is under-powered; k ranged 2–6 and is unstable at that
  sample size. This is a real limitation of the faithful implementation on these recordings.
- 21 of 682 units gave no usable features and are excluded.
- The dip test is sensitive to sample size; the within-recording tests (n = 24–125) are weaker
  than the pooled one.
