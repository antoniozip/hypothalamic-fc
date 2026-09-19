# GLMCC regenerated with correct units — and what it exposes

**Date:** 2026-09-19
**Binary:** `vendor/glmcc-c/glmcc_fixed` (T as a CLI argument; original left intact)
**Driver:** `scripts/regenerate_glmcc.py` -> `results/glmcc_fixed/`, 26/26 datasets

## 1. The C implementation could not process this dataset at all

`WIN 50.0` and `DELTA 1.0` define a +/-50 unit correlogram at 1 unit bins, in whatever unit
the input carries. `read_spikes` filters to `[0, T_SEC*1000)` with `T_SEC` hardcoded at 5400,
so the code expects **0-based milliseconds over a 5400 s recording** — the recordings in
Kobayashi et al. This dataset stores **seconds**, spans 10,000-19,800 s, and starts near 6,000 s:

| input | spikes surviving | effective window |
|---|---|---|
| as-is (seconds) | **100%** | **+/-50 SECONDS** |
| converted to ms | **0%** | n/a — all truncated |

Measured on night3 lightON. Neither mode is usable. `Est_Data.py:37` hardcodes the same
`T = 5400`, so the Python path shares the assumption; `batch_surrogate_glmcc.py`'s docstring
("GLMCC expects 0-based ms times") shows this was known at some point and lost.

`T_SEC` appears only in that filter and the output filename — it never enters the rate baseline,
which comes from the correlogram itself (`n_sp/(2*WIN)`). So making it a CLI argument and feeding
0-based ms yields a genuine +/-50 ms window with no other change in behaviour.

## 2. Densities collapse

| condition | before (+/-50 s window) | after (+/-50 ms) |
|---|---|---|
| lightON | 35-94% | **median 3.44%** (0.19-22.3%) |
| ongoing | 35-94% | **median 33.6%** (16.0-62.6%) |

night3 lightON goes from 77.4% to 0.19% — 2 edges among 1,056 pairs. This confirms that the
implausible network densities flagged in the project review were an artifact of a coincidence
window 1000x too wide.

## 3. But the corrected output is still not trustworthy

**Density is a readout of sample size, not connectivity.**

| | Spearman(spike count, density) | p |
|---|---|---|
| all 26 datasets | **+0.938** | 1.4e-12 |
| within lightON (n=13) | +0.830 | 0.0005 |
| within ongoing (n=13) | +0.714 | 0.0061 |

The relationship holds *within* each condition, so it is not a condition confound. More spikes
collected means more pairs cross significance — which is what a test's power does, not what a
circuit does. ongoing carries 5-15x the spikes of lightON and lands at ~10x the density.

**Every dataset comes out essentially all-inhibitory**: 99.8-100% negative weights across all 26.
Hypothalamus has both excitatory and inhibitory connectivity; an all-inhibitory result in every
animal and both conditions indicates a systematic problem in the fit or in the sign assignment,
not biology.

**The correlograms are data-starved in lightON.** At ~0.65 Hz over ~1,000 s of observed time, the
expected count in a 1 ms bin is ~0.4 for a typical pair; for ongoing (~1.4 Hz over ~15,000 s) it is
~29. Fitting a 102-parameter GLM to a near-empty correlogram is not well posed, and the two
conditions differ ~70x in counts per bin.

## 4. Correction to an earlier claim

`reports/glmcc_sign_asymmetry_diagnosis.md` section 3 proposed that the lightON/ongoing sign
asymmetry arose because GLMCC's baseline scales with observed time and a ~10% duty cycle
under-predicts it. **That mechanism is not supported.** `T_SEC` never enters the baseline. The
empirical result stands — gapping ongoing data did flip its sign distribution, reproducibly, in all
13 animals — but the explanation was wrong. Under a +/-50 s window, 10 s stimulus windows separated
by 80 s gaps put every within-window spike pair inside a single correlogram, which is a different
pathology than the one described. The diagnosis's *observations* remain valid; its *mechanism*
should be disregarded.

## 5. What this means

The question raised in the project review — do these estimates measure connectivity? — now has a
provisional answer: **not yet.** Edge count tracks spike count (rho=+0.94), and sign is
degenerate. Fixing the window was necessary but not sufficient.

Recommended next step, before any further downstream analysis: a positive control. Inject known
synthetic couplings into surrogate trains at realistic rates and recording lengths, then check
whether this pipeline recovers them, at what rate it invents them, and whether it assigns sign
correctly. Until that passes, no network metric computed on these matrices can be interpreted.
