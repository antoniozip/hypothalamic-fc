# Jitter sensitivity: at what timescale does the GLMCC edge set break?

**Date:** 2026-09-19 (corrected same day — see §0) · **Branch:** `analysis/glmcc-corrected`
**Resolves:** Next Action 1 in `METHODS_AND_RESULTS.md` §9
**Artifacts:** `results/jitter_sensitivity.csv`, `figures/main/figJ1_jitter_sensitivity.png`,
`results/glmcc/sensitivity/`

## 0. Correction notice

The first version of this report concluded that the edge set carried **no** spike-timing
information at any scale. **That was wrong**, and so was the ±25 ms jitter result it was built on.

The control scripts — `run_jitter_control.py`, `run_gapped_ongoing.py`,
`run_gapped_ongoing_aligned.py` and the first version of `jitter_sensitivity.py` — invoked the
GLMCC binary with four arguments and no explicit `T`, on spike files written in **seconds**.
GLMCC reads `WIN = 50`, `DELTA = 1` and `tau = 4` in whatever unit the file uses, so those runs
built a **±50 second** cross-correlogram at **1 second** resolution with a 4 s synaptic time
constant and a 1–4 **second** delay scan.

At 1-second resolution a ±25 ms displacement is 1/40 of a single bin. The invariance was
guaranteed by arithmetic, not discovered in the data.

The main connectivity matrices were never affected: `regenerate_glmcc.py` converts to 0-based
milliseconds and passes the true duration. Only the control arms were wrong. Everything below is
re-run through `scripts/glmcc_units.py`, which does the same conversion.

The discrepancy that exposed it: night7 took 157 s per run in the sweep but 8 s in the main
pipeline — a ±50 s window has ~1000× more coincidences to process.

## 1. Answer

**Edges break between 5 and 25 ms** — the monosynaptic range.

| displacement | edges surviving (mean ± SD, n=13) | median | range |
|---|---|---|---|
| 0 ms (identity check) | 100.00 ± 0.00 % | 100.00 | — |
| 5 ms | 16.77 ± 17.81 % | 9.09 | 2.44 – 66.67 |
| 10 ms | 5.60 ± 6.81 % | 2.44 | 0.00 – 19.64 |
| 25 ms | 0.73 ± 1.05 % | 0.00 | 0.00 – 3.14 |
| 50 ms | 2.18 ± 2.97 % | 1.22 | 0.00 – 11.11 |
| 100 ms | 0.42 ± 0.51 % | 0.00 | 0.00 – 1.22 |
| **full within-window shuffle (floor)** | **1.11 ± 1.28 %** | 0.99 | 0.00 – 4.55 |

91/91 runs. The 0 ms level reproduces the baseline exactly in all 13 animals, so the comparison
is against a correct reference.

**By 25 ms the edge set is already at the shuffle floor.** Displacing spikes further buys
nothing, because there is nothing left to destroy. That is the signature of an edge set resting
on structure finer than 25 ms.

Absolute counts tell the same story: the aligned-ongoing baseline has 2,758 edges across the 13
animals; jittered data produces 457 at 5 ms, 225 at 25 ms, 134 after a full shuffle — and of
those, only 242 / 57 / 30 respectively coincide with a baseline edge.

Wilcoxon signed-rank against a 100% null rejects at every level from 5 ms up (p = 2.4e-04, the
floor for n = 13). Here the effect size is not in tension with the p-value: the median loss at
25 ms is 100 percentage points.

## 2. Per-animal

| animal | baseline edges | 5 ms | 25 ms | shuffle |
|---|---|---|---|---|
| day1 | 509 | 17.68% | 2.36% | 1.38% |
| day2 | 35 | 8.57% | 0.00% | 0.00% |
| day3 | 82 | 2.44% | 0.00% | 0.00% |
| day4 | 202 | 15.35% | 0.99% | 0.99% |
| day5 | 167 | 4.79% | 1.20% | 1.80% |
| day6 | 337 | 2.67% | 0.30% | 1.78% |
| night1 | 69 | 4.35% | 1.45% | 1.45% |
| night2 | 22 | 9.09% | 0.00% | 4.55% |
| night3 | 6 | 66.67% | 0.00% | 0.00% |
| night4 | 9 | 22.22% | 0.00% | 0.00% |
| night5 | 21 | 28.57% | 0.00% | 0.00% |
| night6 | 56 | 30.36% | 0.00% | 1.79% |
| night7 | 1,243 | 5.23% | 3.14% | 0.72% |

The 5 ms column is noisy because several animals have very few baseline edges — night3's
"66.67%" is 4 of 6. The 25 ms and shuffle columns are stable across the whole range.

## 3. What this does and does not license

**It does** establish that the aligned-ongoing edge set depends on millisecond spike timing, and
that ±25 ms destroys it. Both jitter controls now pass.

**It does not** rehabilitate everything. Spike count still gates *detection*: pairs receiving an
edge have a far larger spike-count product than those that do not (day3 median 1,159,561 vs
1,674). That is GLMCC's `J_min ∝ 1/sqrt(cc₀)` threshold, unchanged by this result. So which pairs
are testable is set by firing rate, even though whether a testable pair gets an edge is set by
timing. The density-vs-spike-count confound in `METHODS_AND_RESULTS.md` §4.1 stands.

## 4. Knock-on corrections

Re-running the three control scripts with the correct time base changes them all:

| control | seconds (invalid) | milliseconds (correct) |
|---|---|---|
| gapped ongoing | 43,285 edges, 0.0–0.2% negative | 3,463 edges, 70.4–100% negative |
| aligned ongoing | 41,842 edges, 0.0–0.1% negative | 2,758 edges, 76.2–100% negative |
| jitter control | 84,005 edges, 0.0–2.2% negative | 446 edges, 50–100% negative |

**The "sign flips with recording geometry" result is dead.** It compared continuous `ongoing`
(correct ms) against gapped arms (seconds). With both in milliseconds the gapped and aligned arms
are 70–100% negative, the same direction as continuous `ongoing` at 97–100%. There is no flip.
`METHODS_AND_RESULTS.md` §5.3 is withdrawn.

## 5. Limitations

- One arm (aligned-ongoing). Chosen because its geometry is fixed across levels. `lightON` was
  not swept.
- One jitter draw per level, seed 42. Reported spread is between animals.
- Several animals have <25 baseline edges, so their per-level percentages are coarse. The
  aggregate is carried by day1, day6 and night7.
- The 50 ms level sits slightly above 25 ms and 100 ms (2.18% vs 0.73% and 0.42%). With a mean of
  2–3 edges per animal at these levels this is sampling noise around the floor, not structure.
- Mechanism panel (spike-count gating) covers 3 of 13 animals.

## 6. Next actions

1. **Next Action 1 is closed, positively.** The edge set is timing-dependent. The synaptic
   interpretation is available again.
2. **Re-run anything derived from the control arms.** The duty-cycle argument, the geometry
   argument, and any figure or table quoting them need regenerating from the corrected outputs.
3. **Next Action 3 is now worth doing on its merits** rather than as a coup de grâce. A
   within-window-shuffle null replaces the near-saturated shuffle-ISI null, and the ~1% shuffle
   survival predicts the real edges will largely pass it.
4. **Audit every other call site of the GLMCC binary** for the same four-argument pattern.
   `scripts/glmcc_units.py` now centralises the conversion; `scripts/run_glmcc_c_missing.sh`
   still calls the binary directly and should be checked.
