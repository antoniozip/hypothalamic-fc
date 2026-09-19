# Connectivity regenerated with the fixed C port

**Date:** 2026-09-19 · scope: `results/glmcc_fixed/` + the three C-generated controls

The C GLMCC port's sign defect is fixed (`reports/glmcc_c_sign_bug.md`). Everything the old
binary produced has been recomputed. Downstream analysis — graph metrics, the R tier1/tier2
work, LME tables, figures — has **not** been re-run.

## The port was verified on this data, not only on fixtures

Before trusting the output, the port was checked against Kobayashi's Python on a real
recording (night5_lightON, 24 units, 276 ordered pairs, prepared exactly as
`regenerate_glmcc.py` prepares it):

| | value |
|---|---|
| max \|C − Python\| | 1e-6 (the CSV print precision) |
| identical support | yes |
| negative / positive | 8 / 5 in both |

One incidental finding: **`Est_Data.py` cannot process this dataset as shipped.** Pairs with
no coincident spikes within ±50 ms give `rate = 0`, and `init_par` calls `math.log(0)`:

```
File "glmcc.py", line 285, in init_par
    par = math.log(rate)*par
ValueError: math domain error
```

The C port skips those pairs and leaves the weight at zero. The comparison above required
patching the Python to do the same. Anyone re-running the reference on this data will hit this.

## What changed

### Connectivity, `results/glmcc_fixed/` (26 matrices, 82 s)

| | edges | negative |
|---|---|---|
| before | 22,447 | 22,446 (100.0%) |
| after | 23,172 | 22,579 (97.4%) |

Positive edges went from 1 to 593. Per-dataset the negative fraction now ranges 61.5–100%
rather than a flat 100%.

The remaining negative dominance is not the old clamp artifact. Dumping J for night5_ongoing
(552 directed couplings) shows nothing pinned at either bound — 0 at −3, 0 at +5 — but the
whole distribution sits below zero: median J = −0.461, 79% negative, 5th–95th percentile
−0.816 to +0.309, against a median detection threshold of 0.571. Across all 26 matrices 0.6%
of edges land on the clamp, which is the normal rate for a converged fit rather than the
pre-fix pile-up. So this is a genuine negative offset in the estimates, not a failure to
converge — and the Python reference reproduces it exactly on the same data.

### Gapped-ongoing control, `results/glmcc/adj_*_ongoing_gapped.csv`

| | edges | negative |
|---|---|---|
| before | 25,978 | 14–42% per animal |
| after | 43,285 | **0.0–0.2%** |

### Aligned-ongoing control, `results/glmcc/adj_*_ongoing_aligned.csv`

| | edges | negative |
|---|---|---|
| before | 19,696 | 15–41% per animal |
| after | 41,842 | **0.0–0.1%** |

### Jitter control, `results/glmcc/adj_*_{lightON,ongoing_aligned}_jittered.csv`

26 runs, all completed.

## Two results that need your judgement

**1. The sign is set by recording geometry, not by the data.** Continuous ongoing comes out
97–100% negative; the same spikes gapped to the lightON duty cycle come out ~100% positive.
`run_gapped_ongoing.py`'s docstring already recorded this ("gapping ongoing flips it from
98.6% negative to 85.2% positive"); with the corrected estimator the flip is close to total.
Both extremes are implausible as biology, and they bracket the condition contrast the
manuscript rests on.

**2. The jitter negative control now fails outright.** Jittering every spike by ±25 ms
destroys 1–5 ms synaptic structure while leaving slow rate co-modulation intact, so genuine
synaptic edges should disappear. Comparing each aligned-ongoing matrix with its jittered twin
— identical window geometry, identical 45 s offset, only `jitter_ms` differs:

| animal | real | jittered | overlap | surviving |
|---|---|---|---|---|
| day1 | 15,721 | 15,722 | 15,721 | 100.0% |
| day2 | 870 | 870 | 870 | 100.0% |
| day3 | 1,508 | 1,507 | 1,507 | 99.9% |
| day4 | 3,356 | 3,355 | 3,355 | 100.0% |
| day5 | 4,774 | 4,773 | 4,773 | 100.0% |
| day6 | 2,824 | 2,823 | 2,823 | 100.0% |
| night1 | 924 | 924 | 924 | 100.0% |
| night2 | 978 | 977 | 977 | 99.9% |
| night3 | 1,029 | 1,028 | 1,028 | 99.9% |
| night4 | 1,077 | 1,076 | 1,076 | 99.9% |
| night5 | 530 | 531 | 530 | 100.0% |
| night6 | 2,452 | 2,449 | 2,449 | 99.9% |
| night7 | 5,799 | 5,794 | 5,794 | 99.9% |
| **TOTAL** | **41,842** | **41,829** | **41,827** | **100.0%** |

Essentially every edge survives. Whatever the gapped/aligned matrices are measuring, it is not
millisecond-scale coupling. Note `retention_pct` in `results/jitter_control_summary.csv` is the
fraction of *spikes* kept (~99.87%), not edge survival — the table above is the edge statistic.

This is not a porting question: the estimator now agrees with the reference to every printed
digit. It is a question about what GLMCC can measure on recordings with this geometry and these
spike counts.

## Also worth noting

The corrected positive control puts the **lightON regime below the detection floor** — at 788
spikes/neuron, recall is 0%. The previous 50% "recall" at 0% sign accuracy was the inverted
port detecting excitatory couplings and labelling them inhibitory.

Pre-fix outputs were backed up before being overwritten; ask if you need them restored.

## Not done

Graph metrics, `src/r/tier1_analyses.R` / `tier2_analyses.R`, the LME tables, community and
hub results, and Figures 1–4 all still derive from the pre-fix matrices.

---

## A provenance question worth resolving before trusting the sign-asymmetry diagnosis

`reports/glmcc_sign_asymmetry_diagnosis.md` builds its case on the raw matrices in
`results/glmcc/adj_{animal}_{condition}.csv`: lightON 0.9% negative vs ongoing 92.9% negative,
attributed to the ~10% lightON duty cycle. Those files predate the port fix.

Counting weights pinned exactly at the J = −3 clamp (i.e. −3 × c_I = −1.836):

| set | edges | at the clamp |
|---|---|---|
| `results/glmcc/adj_*_lightON.csv` | 32,633 | 3 (**0.0%**) |
| `results/glmcc/adj_*_ongoing.csv` | 33,418 | 6,280 (**18.8%**) |
| `results/glmcc_fixed/` regenerated today | 23,172 | 132 (0.6%) |
| `results/glmcc_fixed/` pre-fix | 22,447 | 133 (0.6%) |

The clamp is reachable in correct fits — 0.6% of edges land there — so it is not by itself a
broken/fixed discriminator, and the pre-fix regeneration does **not** show enrichment. But
18.8% in the original ongoing matrices against 0.0% in the original lightON matrices is a 30×
asymmetry *within the very files the diagnosis compares*. The two conditions were evidently not
produced under the same conditions.

The `.meta.json` files that survive for those paths are placeholders from failed runs
(`GLMCC exited with code 2`, `GLMCC timed out (86400s limit)`), so provenance cannot be
recovered from the repository. Some were run with `vendor/GLMCC/Est_Data.py`, some with the C
port, on dates spanning 2026-06-11/12 onward.

**This does not disprove the duty-cycle mechanism.** It does mean the headline contrast in
section 1 of that report conflates a condition difference with a provenance difference, and the
confirmatory experiments in sections 6–7c were themselves run through `vendor/glmcc-c/glmcc`,
which carried the sign defect. Both conditions need recomputing under one estimator before that
diagnosis can be relied on — and before it is used to justify anything in the manuscript.
