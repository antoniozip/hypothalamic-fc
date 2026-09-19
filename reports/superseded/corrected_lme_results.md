# First-level LME on duty-cycle-matched inputs

**Date:** 2026-09-19
**Inputs:** `results/glmcc_dc/` (matched) and `results/glmcc_raw/` (control arm)
**Model:** `lmer(metric ~ condition * region + (1 | animal/neuron))`, n = 1330 neurons, 12 animals
paired (night2 lightON excluded: its only GLMCC matrix is 50x50 for a 34-unit recording)

## Design

Three arms isolate each change from the original analysis:

| arm | weights | ongoing side |
|---|---|---|
| OLD | surrogate-validated | native ~89% duty cycle |
| RAW | raw GLMCC | native ~89% duty cycle |
| DC | raw GLMCC | gapped to the true stimulus windows, ~10% duty cycle |

OLD vs RAW isolates the effect of skipping validation; RAW vs DC isolates duty-cycle matching.

Validation was skipped in the corrected arm for two reasons: shuffle-ISI surrogates are invalid on
gapped trains (shuffling intervals scatters the ~80 s inter-window gaps and destroys the structure
being matched), and the validation removed only 0.5% of edges in the first place.

## Result 1: surrogate validation changes nothing

| metric (condition coefficient) | OLD (validated) | RAW (unvalidated) |
|---|---|---|
| node_strength | 122.6 (p=0.24) | 123.9 (p=0.24) |
| clustering_coefficient | -0.0217 (p=0.76) | 0.0091 (p=0.90) |
| local_efficiency* | 1.341 (p=9.5e-10) | 1.354 (p=2.3e-11) |
| hub_score | 0.0391 (p=0.80) | -0.0167 (p=0.91) |

This confirms at the model level what the edge counts already showed: the surrogate + BH-FDR layer
is not doing meaningful work.

## Result 2: duty-cycle matching changes the conclusions

Marginal condition effect (lightON - ongoing), `emmeans` over regions:

| metric | duty MISmatched | duty MATCHED | change |
|---|---|---|---|
| node_strength | +287 (p=2.3e-46) | +226.4 (p=1.4e-33) | attenuated 21% |
| local_efficiency* | +2.504 (p=2.6e-184) | +2.145 (p=6.8e-134) | attenuated 14% |
| clustering_coefficient | -0.0125 (p=0.35) | **+0.119 (p=8.1e-10)** | null -> significant, sign flips |
| hub_score | +0.050 (p=0.077) | **+0.207 (p=5.2e-12)** | null -> significant |

Matching does not simply shrink effects. The unmatched `ongoing` weights were noise-dominated
(92.9% negative, median -0.148), which washed out topological contrasts; on matched inputs,
clustering and hub structure separate the conditions where they previously did not.

Mean per-neuron node strength: `ongoing` moves from -46.10 (native, sign-inverted by the baseline
artifact) to +12.59 (matched), against lightON's +330.

## Result 3: the region-specific finding is robust

`condition x region` for node_strength at **VM-thalamus** survives matching (est=449.4, p=7.3e-06),
consistent with the encoding-fix refit which also isolated VM-thalamus (p=2.7e-06). The DMH result
originally reported in `NEXT_STEPS.md` does not survive either correction.

## Caveats, in order of severity

1. **`local_efficiency` has since been fixed** (marked * above; values are post-fix). It had two
   self-cancelling errors: weights were used directly as path lengths instead of 1/|w|, and the
   metric returned `1/mean(1/d)` (harmonic mean of distances) rather than Latora & Marchiori's
   `mean(1/d)`. The two forms coincide for a uniform complete neighbour subgraph, which is why the
   bug was not visible. Correcting it roughly tripled the effect size without changing its sign, so
   the old metric was partly cancelling the real signal. Guarded by
   `tests/test_local_efficiency.R`. The OLD-arm row is pre-fix and kept only to show that
   validation changes nothing.
2. **Omnibus vs simple effects diverge sharply.** Type-II ANOVA condition terms are p=1e-163 to
   1e-239 while reference-region coefficients are null for three of four metrics. With a strong
   `condition x region` interaction the omnibus test is not an interpretable summary; the emmeans
   marginals above are. Reporting the ANOVA alone would badly overstate the effect.
3. **Weights are raw, not validated** (justified above), so these numbers are not directly
   comparable to any previously published figure from this pipeline.
4. **The duty-cycle match is close but not exact.** The control samples light-off periods with the
   same window geometry, but lightON additionally carries real stimulus drive. The jitter control
   shows the residual is fine-timescale; +/-25 ms jitter cannot separate synaptic coupling from
   precisely-locked common drive, so the residual's origin is still open.
5. **No multiplicity correction** across the 4 metrics x 20 terms.

## Reproduce

```
python scripts/assemble_dutycycle_matched.py
Rscript src/r/compute_graph_metrics.R --animal <a> --condition <c> --estimator glmcc_dc
Rscript src/r/stats_first_level.R --estimator glmcc_dc
```

Outputs: `results/stats_first_level_glmcc_dc.csv`, `results/stats_first_level_glmcc_raw.csv`,
`results/dutycycle_matched_inputs.csv`.
