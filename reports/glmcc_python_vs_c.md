# The C port inverts connection sign; the Python reference is correct

**Date:** 2026-09-19 · `scripts/glmcc_python_crosscheck.py`

## First: a correction to the previous report

The positive control in `reports/glmcc_positive_control.md` wrote its fixture with **1-indexed**
cell files (`cell1..cellN`). Both GLMCC implementations loop `0..n-1`, and the real `DATA*/`
directories are 0-indexed (`cell0.txt..cell{n-1}.txt`). So the harness left `cell0` missing and
shifted every ground-truth index by one. **Those numbers were invalid and are superseded here.**

The error is mine, in the test, not in the pipeline — and it is the same class of off-by-one
indexing fault that this project has already been bitten by twice.

## Corrected positive control (0-indexed, C implementation)

| scenario | spk/neuron | density | recall | precision | sign correct | detected negative |
|---|---|---|---|---|---|---|
| generous | 36,448 | 10.00% | 50.0% | 23.0% | 0.0% | 100.0% |
| ongoing-like | 25,553 | 8.05% | 50.0% | 28.6% | 5.0% | 98.6% |
| lightON-like | 788 | 4.83% | 50.0% | 47.6% | 0.0% | 100.0% |

Recall is exactly 50% in every scenario, with sign accuracy ~0%. Breaking that down by polarity
on the generous fixture explains why:

```
true excitatory: 20  ->  detected 20 (100%), ALL reported negative
true inhibitory: 20  ->  detected  0 (0%)
```

The C implementation recovers **every** excitatory connection and labels each one inhibitory,
while detecting **no** inhibitory connections at all. That is a sign inversion, not weak
sensitivity.

## Python vs C on one identical fixture

12 neurons, 12,179 spikes/neuron, 16 true connections (8 E / 8 I), same spike trains fed to both:

| implementation | density | recall | precision | **sign correct** | detected negative | runtime |
|---|---|---|---|---|---|---|
| C (`glmcc_fixed`) | 16.67% | 56.2% | 40.9% | **11.1%** | 100.0% | 0.3 s |
| Python (vendored) | 12.88% | 62.5% | 58.8% | **100.0%** | 35.3% | 46.3 s |

Python recovers more connections, with better precision, and gets **every sign right**, returning
a plausible E/I mix. The 150x speed advantage of the C port comes with a broken sign.

## Conclusion

**GLMCC the method works on data like this.** Kobayashi's Python implementation achieves 62.5%
recall at 58.8% precision with perfect sign accuracy at this operating point. The earlier
conclusion that "the estimator does not measure connectivity" was wrong; it was the port.

**Every C-generated matrix in this repository is affected**, which includes:

- `results/glmcc_fixed/` — all 26 datasets regenerated today
- the duty-cycle experiment matrices (`*_ongoing_gapped`, `*_ongoing_aligned`, `*_jittered`)
- in the original pipeline, the six `ongoing` matrices recomputed on 2026-06-11/12, which
  `NEXT_STEPS.md` attributes to the C implementation

The 99.8-100% negative weights seen across all 26 regenerated datasets are now explained: the port
reports excitatory couplings as inhibitory and misses genuine inhibitory ones.

## Next steps

1. **Find the sign bug in `vendor/glmcc-c/glmcc.c`.** Candidate: the C code pairs `Jp` with
   `Jmin[0]` (line ~419) whereas `Est_Data.py` pairs `J_+` with `Jmin[1]` and `J_-` with `Jmin[0]`
   — the thresholds may be swapped. The cross-correlogram convention (`diff = c1[j] - c2[i]`)
   should also be checked against the Python, since a reversed convention swaps `Jp`/`Jm`.
2. **Until then, regenerate with Python.** It is ~150x slower, so a 126-unit animal is a long run,
   but correctness dominates here.
3. Re-run the positive control against any fixed C build before trusting it again.
