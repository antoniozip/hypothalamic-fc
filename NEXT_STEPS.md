# Next steps

The previous contents of this file described a pipeline whose connectivity was produced by a
C GLMCC port that inverted the sign of every coupling. That work was discarded on 2026-09-19.

**Current record: `METHODS_AND_RESULTS.md`.**
**Branch: `analysis/glmcc-corrected`.**
**Why the old results went away: `reports/superseded/README.md`.**

## The state in one paragraph

The port is fixed and verified against Kobayashi's reference to the printed precision, in all
four modes, on fixtures, on the reference `simulation_data`, and on this dataset
(`tests/test_glmcc_c_parity.py`, 20 PASS / 0 FAIL). All 26 matrices, three control arms, edge
validation, graph metrics and the tier 1/2 analyses were recomputed. The condition effect is
large and significant everywhere. Two negative controls fail.

## Blocking questions, in order

1. ~~Why does ±25 ms jitter remove no edges?~~ **ANSWERED — `reports/jitter_sensitivity/`.**
   There is no timescale at which it breaks. Survival is 99.98/99.97/99.95/99.93/99.92% at
   5/10/25/50/100 ms, and **99.30% under a full within-window shuffle** that randomises every
   spike time. Weights stay correlated at r = 0.96 through the shuffle. Edge presence is set by
   the pair's spike-count product (day3: median n_i*n_j 1,159,561 with an edge vs 1,674 without,
   p = 8.5e-297) — GLMCC's J_min threshold, not coupling. **The connectivity is not synaptic.**
2. **Retire the raw `lightON` vs `ongoing` contrast.** Sign follows recording geometry —
   continuous ongoing is 97–100% negative, the same spikes gapped to the lightON duty cycle are
   ~100% positive. Use `lightON` vs aligned-ongoing, which matches window count, duration and
   period.
3. **Replace the shuffle-ISI null with the within-window shuffle** (already implemented in
   `scripts/jitter_sensitivity.py`). The current null rejects 0.0% of lightON and 4.6% of ongoing
   edges, so it constrains nothing; the within-window shuffle is a null this edge set
   demonstrably fails. **This is now the step that decides whether there is a paper.**
4. **Handle the spike-count confound.** `ongoing` has 10–30x more spikes than `lightON`, and
   GLMCC's threshold shrinks as 1/sqrt(cc0), so density tracks spike count directly.
5. **Respecify the second-level model.** It is rank deficient (76 columns dropped) and every
   contrast returns `nonEst`.

Item 1 is settled and the answer is negative. Item 3 now decides the outcome: if a
within-window-shuffle null leaves too few edges to analyse -- which 99.30% survival predicts --
then GLMCC at this operating point does not support a connectivity claim on this dataset, and
that is the result. Items 2, 4 and 5 matter only if item 3 leaves something standing.

## Deliberately not re-run

Community detection, participation coefficient, hub stability, circos plots, main manuscript
figures and tier-3 analyses. Regenerating them before (1) is settled would be wasted effort.

## Reproducing

See section 10 of `METHODS_AND_RESULTS.md` for the full command sequence, runtimes, seeds and
the `pbkrtest` installation workaround.
