# Why the recordings differ: an acquisition change mid-study

**Date:** 2026-09-20 · **Branch:** `analysis/glmcc-corrected`
**Follows:** `reports/ei_differential/waveform-reclassification.md`

The re-classification failed because median peak-to-peak (PP) duration varies 2.6× between
recordings. This is the diagnosis of that variation.

## Ruled out: windowing, alignment, sampling

| check | result |
|---|---|
| waveform length | **81 samples in all 13 recordings** |
| trough position | **sample 40 in all 13** — identically aligned, centred |
| peak at the window edge | **0.0% of units, every recording** — no truncation |

The window is adequate, the alignment is identical, and the sampling is common. The PP
difference is real shape, not measurement geometry.

## The cause: a step change between 2018-01-11 and 2018-01-31

Ordering the recordings by date:

| date | animal | n | median PP (ms) | peak/trough ratio | trough depth |
|---|---|---|---|---|---|
| 171019 | day1 | 126 | 0.750 | 0.228 | 35.1 |
| 171207 | night1 | 31 | 0.700 | 0.226 | 28.5 |
| 171208 | night2 | 34 | 0.800 | 0.242 | 27.6 |
| 171213 | night3 | 33 | 1.050 | 0.215 | 24.6 |
| 180110 | night4 | 35 | 0.900 | 0.250 | 24.6 |
| 180111 | night5 | 24 | 0.700 | 0.225 | 29.2 |
| — | — | — | — | — | — |
| 180131 | day2 | 31 | 0.450 | 0.400 | 11.4 |
| 180221 | night6 | 54 | 0.550 | 0.312 | 22.2 |
| 180228 | night7 | 79 | 0.300 | 0.258 | 10.2 |
| 180302 | day3 | 47 | 0.500 | 0.273 | 20.4 |
| 180419 | day4 | 61 | 0.550 | 0.240 | 20.4 |
| 180420 | day5 | 71 | 0.550 | 0.273 | 19.2 |
| 180423 | day6 | 56 | 0.450 | 0.250 | 21.0 |

- Spearman(date order, median PP) = **−0.715, p = 0.0060**
- early (≤180111, n=6) median PP **0.775 ms** vs late (≥180131, n=7) **0.500 ms**,
  Mann–Whitney **p = 0.0031**

Three things move together across that boundary:

1. **PP shortens** — 0.775 → 0.500 ms
2. **Trough depth drops** — 24.6–35.1 → 10.2–22.2
3. **Positive overshoot grows relative to the trough** — 0.215–0.250 → 0.240–0.400

That combination is the signature of a **more aggressive high-pass filter**. Stronger
high-passing differentiates the signal: it attenuates the slow trough, speeds repolarisation and
lifts the positive overshoot. A pure biological or probe-depth difference would not move all
three in that coordinated way.

Something in the acquisition or preprocessing chain changed between 2018-01-11 and 2018-01-31.

## This corrects an earlier claim of mine

`reports/ei_differential/analysis-report.md` states that six animals were "never classified",
inferring that from their 0–3.2% inhibitory fraction. **That was wrong.** The classification ran
on all 13; it returned almost nothing on the early block because the fixed threshold sits above
that block's distribution.

Applying the original PP threshold (< 0.65 ms) to recomputed PP reproduces the pattern:

| block | % of units with PP < 0.65 ms | original % labelled inhibitory |
|---|---|---|
| early (6 recordings) | 0.0 – 29.0% | 0.0 – 3.2% |
| late (7 recordings) | 59.6 – 86.2% | 41.4 – 51.1% |

Spearman between the two columns = **+0.722, p = 5.3e-03, n = 13**.

So the "classifiable / not classified" split I described was really **one fixed threshold meeting
two different waveform distributions**. The 44.6% inhibitory fraction in the late block and the
~0% in the early block are the same artifact seen from both sides.

The practical consequence is unchanged — the shipped labels are not usable — but the reason is
different and more tractable: it is a correctable acquisition difference, not missing work.

## Confirmed: the early block is missing the 600 Hz high-pass

**Information from the experimenter (2026-09-20):** the acquisition filter was nominally
**600–6000 Hz** for every session, and no pre-spike-sorting data survives.

That is decisive for the interpretation, because 600 Hz is already an aggressive high-pass — a
spike filtered there should show PP near 0.3–0.5 ms. The late block's 0.500 ms matches. The early
block's 0.775–1.05 ms does not; it is what a substantially lower corner produces.

Tested directly (`scripts/waveform_harmonise.py`): for each recording, sweep a compensating
zero-phase 2nd-order Butterworth high-pass over 50–1525 Hz and find the corner that brings its
median PP to the late-block value of 0.500 ms.

| block | compensating corner required (Hz) |
|---|---|
| **early** | 550, 600, 625, 625, 675, **775** |
| **late** | 0, 0, 50, 50, 75, 350, 450 |

The six early recordings independently converge on **~600 Hz — the nominal setting** — while the
late seven need essentially none. The grid was 50–1525 Hz and nothing constrained the answer
toward 600; that the early block lands there is strong corroboration.

**Conclusion: the 600–6000 Hz filter reached the late seven recordings and effectively did not
reach the early six.** Applying the compensation collapses the between-recording spread in median
PP from **2.62× to 1.11×** (0.400–1.050 ms → 0.450–0.500 ms).

### Harmonising helps, but does not rescue the classification

Re-running the spikeMAP clustering on harmonised waveforms:

- **The block effect is gone.** Inhibitory fraction early 21.7% vs late 31.0%, Mann–Whitney
  **p = 0.366** (it was p = 0.0031 on the raw features). The systematic artifact is removed.
- **The classification is still not usable.** 28.8% inhibitory overall against a 5.25% reference,
  per-animal still 0–42.6%, full width still unimodal (dip p = 0.99) with only peak-to-peak
  bimodal, and Calinski–Harabasz now selects k = 6 with clusters as small as n = 3.

So the confirmed artifact was real and is correctable, but correcting it is not sufficient. The
likely ceiling is the correction itself: compensating a *mean, already-filtered* waveform is not
equivalent to filtering the raw trace. The original filter's phase response is baked in, spike
averaging has smoothed it, and the two recordings needing the largest compensation (night3 at
1.050 ms, night4 at 0.900 ms) are also the two that still return 0% and 5.7%.

**Because the pre-sorting data is gone, that ceiling cannot be lifted.** The correction available
is the best one obtainable from what survives.

## What this makes possible

The variation has a named, bounded cause, which is better than the previous position.

1. **Identify the change — not recoverable from the files.** The `Exp` struct was inspected and
   carries `Fs`, `Fstart`, `Fstop`, `Files`, `LightON/OFF`, `StimON/OFF`, `units`. `Fs` is 20000
   in all 13 (confirming the sampling rate) and `Fstart` is 0.0 in all 13; `Fstop` is a
   duration-like quantity, 326–471, uncorrelated with median PP (Spearman −0.127, p = 0.68).
   **No filter settings are stored.** The change has to be identified from lab notes or the
   original acquisition software configuration for that 20-day window.
2. ~~Re-filter both blocks from raw traces.~~ **Not possible** — no pre-spike-sorting data
   survives. The waveform-level compensation above is the best available substitute and it is
   not sufficient.
3. **Use the harmonised labels only with the block caveat attached.**
   `data/processed/ei_harmonised_labels.csv` removes the block artifact but still yields an
   implausible 28.8% inhibitory fraction. It is better than the shipped labels and still not good
   enough to support an E/I claim.
4. **Re-check other waveform-derived quantities** for the same block structure before trusting
   them across the boundary.

## Does it affect the connectivity results?

Not directly. GLMCC uses spike **times**, not waveforms, and a filter difference changes waveform
shape rather than whether a spike is detected.

But spike amplitude roughly halves across the boundary, and a missing 600 Hz high-pass in the
early block would also leave more low-frequency energy in the traces that were sorted. Both can
change sorting yield and unit isolation. **This cannot be tested: no pre-spike-sorting data
survives.** It therefore stands as an untestable limitation on pooling the two blocks — worth
stating explicitly in any write-up rather than leaving implicit.

The rate-matched condition result (`reports/rate_matched/`) is between conditions *within* each
animal, so it is unaffected by a between-block difference.

## Limitations

- The filter interpretation is inferred from the direction of three coordinated changes. It could
  **not** be confirmed against acquisition settings, because none are stored in the files (see
  above). An amplifier, probe or sorting-template change
  could produce a similar signature.
- n = 6 and 7 recordings per block; the Mann–Whitney rests on 13 points.
- The boundary is placed between 180111 and 180131 because that is where the gap falls. With no
  recordings between those dates, the change could have happened anywhere in that 20-day window.
- Only mean waveforms were available (`Units(i).wf`); per-spike variability was not examined.
