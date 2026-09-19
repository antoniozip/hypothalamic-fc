# C GLMCC sign bug: partially fixed, root cause localised

**Date:** 2026-09-19 · `vendor/glmcc-c/glmcc_fixed.c`

## Fixed: swapped coefficient/threshold pairing

`Est_Data.py` writes `par[NPAR-1]` as J_+ and `par[NPAR-2]` as J_-, then assigns

```
W[i][j] = calc_PSP(J_+, Jmin[1]*scale)      W[j][i] = calc_PSP(J_-, Jmin[0]*scale)
```

The port had both slots reversed, so each direction received the other's coefficient and
threshold. Correcting it improves the positive control measurably:

| | recall | precision | sign correct |
|---|---|---|---|
| before | 56.2% | 40.9% | 11.1% |
| after | **68.8%** | **50.0%** | 27.3% |

## Not fixed: J carries the wrong sign

The port still returns 100% negative weights where Python returns 35.3%. Comparing the two
matrices element-wise on an identical fixture (14 pairs non-zero in both):

| | value |
|---|---|
| corr(\|C\|, \|Python\|) | **+0.900** |
| corr(C, Python) | **-0.857** |
| median \|C\| / \|Python\| | **0.176** |
| C signs | 14 negative, 0 positive |
| Python signs | 6 negative, 8 positive |

Magnitudes track the reference closely, so the GLM fit is working. The signed correlation is
inverted, and the magnitude ratio identifies the mechanism:

```
c_I / c_E = 0.612 / 2.532 = 0.242
```

If J reaches the PSP assignment with its sign flipped, a genuinely positive J falls through
`if (J > Jth)` into `else if (J < -Jth)` and is multiplied by `c_I` rather than `c_E`, yielding a
negative weight at ~0.24x the correct magnitude. Observed ratio 0.176 and inverted sign match
that prediction.

**The defect is therefore upstream of the PSP branch — J itself is negated relative to the
reference.** The assignment code is correct as written.

## Where to look next

1. **Correlogram lag convention.** `cross_corr` computes `diff = c1[j] - c2[i]` with `c1 = cell_i`,
   `c2 = cell_j`. If Python's `cross_correlogram` uses the opposite order, the lag axis is mirrored.
2. **Synaptic basis sign in the GLM design.** The exponential kernel `exp(-(t-delay)/tau)` entering
   `calc_log_post` / `lm_optimize` is the other place a sign could invert.

A useful next diagnostic is to have both implementations dump J for the same pair — the C would
need a debug print, since it writes only W while Python also writes `J_py_*.txt`.

## Status

The C port remains unusable for any sign-dependent result. `results/glmcc_fixed/` is marked with a
README to that effect. Use the vendored Python until the port passes
`scripts/glmcc_positive_control.py`.
