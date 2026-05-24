# Hypothalamic FC Pipeline — Complete

## ✅ Completed — All Sprints

### Sprint 0 — Cleanup & ground-truth
- Pearson retired to `legacy/`, region table canonical (`data/processed/neurons.csv`)
- Filename hygiene: all 13 lightON + 8/13 ongoing GLMCC results at top level
- Repository scaffold created (`data/`, `results/`, `src/`, `figures/`, `legacy/`)
- Environments locked (`renv` + `uv`)

### Sprint 1–2 — GLMCC + TE pipelines
- 7 Python scripts: surrogates, GLMCC wrapper, edge validation, TE, CCG validation
- 5 R scripts: graph metrics, density, first-level stats, second-level stats, hub stability
- Snakemake DAG defining the full workflow

### Sprint 3 — Graph metrics & density
- Node strength, clustering coefficient, local efficiency, hub score for 21 pairs
- Region-pair density with count-bias correction

### Sprint 4 — Mixed-effects statistics
- First-level: `lmer(metric ~ condition + (1|animal/neuron))` — 4 metrics
- Second-level: `lmer(density ~ condition × day_night + (1|animal))`
- Hub stability: Spearman ρ across conditions with Wilcoxon test

### Sprint 5 — Temporal structure
- Adaptation: 7/13 animals show significant firing rate changes early→late
- Early/late evoked window: 0-500ms vs 500ms-10s within each 10s pulse
- Reported as finding; main analysis kept as-is

### Sprint 6 — Workflow & figures
- Snakemake DAG (52 leaf jobs)
- Figure cleanup: 21 `*_v2.png`, `*_fixed.png` variants archived to `legacy/figures/`
- Hub stability figure at `figures/main/fig4_hub_stability.png`

### Additional — TE pipeline unblocked
- Pure numpy TE implementation (`src/py/transfer_entropy_np.py`) replaces pyinform
- All 26 existing TE results linked to `results/te/`
- Note: numpy TE gives 0 for very sparse 5ms-binned data; existing pyinform results used
- INVESTIGATION NEEDED: TE returning 0 may indicate a fundamental issue with
  5ms binning at hypothalamic firing rates (0.15 Hz mean). Test with:
  * Broader bins (10, 20, 50 ms) to find threshold where signal emerges
  * Synthetic Poisson spike trains at matched rates as positive control
  * Alternative TE estimators (JIDT via subprocess, pyinform 0.5.x)
- See IMPLEMENTATION_PLAN.md Phase 2.3 for detailed investigation steps
- NOTE: validate_edges.py heuristic p-values have been REMOVED (see review-improvements branch).
  All edge validation now uses empirical CCG+shuffle-ISI surrogates.

## Key Results Location

| Output | Path |
|--------|------|
| Edge-validated adjacencies (CCG + shuffle-ISI) | `results/glmcc/validated_adj_*.csv` |
| P-values alongside validated edges | `results/glmcc/p_values_*.csv` |
| Graph metrics | `results/glmcc/metrics_*.csv` |
| Region-pair densities | `results/glmcc/density_*.csv` |
| First-level stats | `results/stats_first_level_glmcc.csv` |
| Second-level stats | `results/stats_second_level_glmcc.csv` |
| Hub stability | `results/hub_stability.csv` |
| Sprint 5 temporal analysis | `results/sprint5/` |
| TE matrices | `results/te/` (symlinked) |
| Pipeline scripts | `src/py/` (12), `src/r/` (5) |
| Legacy (Pearson, old scripts) | `legacy/` |

## Remaining (Collaborator-dependent)
- Manuscript (Sprint 7) — needs Storchi/Brown input
- 5 missing GLMCC ongoing runs (needs GLMCC executable)
- Final figure styling for publication
