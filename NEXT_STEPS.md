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

1. **Why does ±25 ms jitter remove no edges?** 41,827 of 41,842 survive (100.0%). Jitter at
   5/10/25/50/100 ms and find the timescale where the edge set starts to break. Until this has
   an answer, no connectivity claim here is safe.
2. **Retire the raw `lightON` vs `ongoing` contrast.** Sign follows recording geometry —
   continuous ongoing is 97–100% negative, the same spikes gapped to the lightON duty cycle are
   ~100% positive. Use `lightON` vs aligned-ongoing, which matches window count, duration and
   period.
3. **Replace the shuffle-ISI null.** It rejects 0.0% of lightON and 4.6% of ongoing edges, so it
   constrains nothing. Build a jitter null at the timescale from (1).
4. **Handle the spike-count confound.** `ongoing` has 10–30x more spikes than `lightON`, and
   GLMCC's threshold shrinks as 1/sqrt(cc0), so density tracks spike count directly.
5. **Respecify the second-level model.** It is rank deficient (76 columns dropped) and every
   contrast returns `nonEst`.

Items 1 and 2 decide whether there is a paper. Items 3-5 are prerequisites for trusting any
number in the results.

## Deliberately not re-run

Community detection, participation coefficient, hub stability, circos plots, main manuscript
figures and tier-3 analyses. Regenerating them before (1) is settled would be wasted effort.

## Reproducing

See section 10 of `METHODS_AND_RESULTS.md` for the full command sequence, runtimes, seeds and
the `pbkrtest` installation workaround.
