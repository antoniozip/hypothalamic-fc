# TODO — Hypothalamic FC project

Plan derived from the revision proposal (`hypothalamic_FC_revision_proposal.docx`),
with all suggestions accepted and the FC method set restricted to **GLMCC** and
**Transfer Entropy** only. Pearson is dropped from the analysis pipeline and the
manuscript; existing `Pearson/` artefacts are archived for the record.

Animals (13): `171019, 171207, 171208, 171213, 180110, 180111, 180131, 180221,
180228, 180302, 180419, 180420, 180423`
Day/Night vector: `c(D,N,N,N,N,N,D,N,N,D,D,D,D)` → 6 day, 7 night.

Legend: `[ ]` not started · `[~]` in progress · `[x]` done.

---

## Sprint 0 — Cleanup and ground-truth (week 1)

### 0.1 Drop Pearson from the active pipeline
- [ ] Move `Pearson/` to `legacy/Pearson/` (keep, do not delete).
- [ ] Move `analysis_network_pearson.pzfx` and `R_analyses*.R` Pearson blocks to `legacy/`.
- [ ] Strip Pearson references from `glmcc_analyses.R`, plotting scripts and the figure inventory.
- [ ] Search & confirm: `grep -rni "pearson" *.R *.py` should return only `legacy/` matches.

### 0.2 Canonical region table (single source of truth)
- [ ] Pick `regions_match_revision_March2024.xlsx` as canonical → export to `regions.csv`.
- [ ] Move the other three xlsx (`regions_match.xlsx`, `regions_match-DESKTOP-72JS1J9.xlsx`, `regions_match_revision_May2023.xlsx`) into `legacy/regions/`.
- [ ] Produce `neurons.csv` with columns: `animal, electrode_id, neuron_id, region6, region9, day_night, isi_violation_pct, snr, retained`.
- [ ] Acceptance: every neuron in every `DATA*/cell*.txt` has a row in `neurons.csv`; row counts match per-experiment file counts.

### 0.3 Filename hygiene (uniform 13-animal coverage)
- [ ] Reconcile `result_GLMCC_lightON_*.csv` so all 13 animals are present at top level (currently 9/13). For the 4 missing (`171208, 180111, 180228, 180302, 180419` — verify), copy/symlink from `GLMCC/` subfolder.
- [ ] Reconcile `result_GLMCC_ongoing_bis_*.csv` (currently 6/13 — verify).
- [ ] Acceptance: `for d in $DAYS; do test -f result_GLMCC_lightON_${d}.csv || echo MISSING $d; done` prints nothing.

### 0.4 Repository scaffold
- [ ] Create directory layout:
  ```
  data/raw/           # *.mat, DATA*/, *_put_inh.txt
  data/processed/     # spike trains, region tables
  data/surrogates/    # all *_surrogates/
  results/glmcc/      # adjacency, node-metrics, surrogate-validated
  results/te/         # adjacency, node-metrics
  src/r/              # R scripts
  src/py/             # Python scripts
  figures/main/
  figures/supp/
  legacy/             # Pearson, old xlsx, old workspaces
  ```
- [ ] Replace hard-coded absolute paths in `*.R` and `*.py` with `here::here()` (R) and `pathlib.Path(__file__).resolve()` (Py).
- [ ] Move `workspace_*.RData` into `legacy/workspaces/`. Keep only the most recent for reference.

### 0.5 Environment lock
- [ ] R: `renv::init()`; pin `lme4`, `lmerTest`, `emmeans`, `igraph`, `network`, `circlize`, `ComplexHeatmap`, `pzfx`, `here`.
- [ ] Python: `uv init` (or conda env); pin `elephant`, `neo`, `quantities`, `pyinform`, `scipy`, `numpy`, `pandas`.
- [ ] Record GLMCC commit hash (Kobayashi et al. 2019 implementation) in `src/py/GLMCC_VERSION.txt`.

---

## Sprint 1 — Complete the GLMCC pipeline across all 13 animals (week 2)

### 1.1 Surrogate generation (currently 2/13)
- [ ] Refactor `surrogate_validated_connectivity.py`:
  - Remove the hard-coded `/home/antonio/onedrive_cnr/` prefix.
  - Parametrise `--animal`, `--condition {ongoing,lightON}`, `--n-surrogates 100`, `--dt-ms 0.1`.
  - Write into `data/surrogates/<animal>_<condition>/<surr_id>/cell_<i>.txt`.
- [ ] Run for the 11 animals missing a surrogate folder × 2 conditions = 22 jobs.
- [ ] Acceptance: `data/surrogates/<animal>_<condition>/{0..99}/cell_*.txt` exists for all 13 × 2.

### 1.2 GLMCC inference on real + surrogate spike trains
- [ ] Wrap GLMCC call into `src/py/run_glmcc.py` with the same `--animal --condition` interface.
- [ ] For each animal × condition: produce `results/glmcc/adj_<animal>_<condition>.csv` (real) and `results/glmcc/surrogate_adj/<animal>_<condition>/<surr_id>.csv` (100 surrogates).
- [ ] Acceptance: 13 × 2 = 26 real adjacency files + 26 × 100 = 2600 surrogate adjacency files.

### 1.3 Edge-level validation
- [ ] New script `src/py/validate_edges.py`:
  - Per (i,j) pair: empirical p = fraction of surrogates whose |J_ij| ≥ |J_ij,real|.
  - FDR correction (Benjamini–Hochberg, q = 0.05) across all pairs within an animal × condition.
  - Output `results/glmcc/validated_adj_<animal>_<condition>.csv` (zeroed where not significant).
- [ ] Acceptance: 26 validated adjacencies. Report per-animal edge survival rate.

### 1.4 GLMCC sensitivity sweep (Methods supplement)
- [ ] Pick 2 representative animals (1 day, 1 night). Sweep `bin_width ∈ {0.5, 1, 2} ms` and `delay window ∈ {±10, ±25, ±50} ms`.
- [ ] Plot node-metric correlation across parameter settings → `figures/supp/glmcc_sensitivity.png`.

---

## Sprint 2 — Complete the Transfer Entropy pipeline (week 3)

### 2.1 TE recompute, parametrised and batched
- [ ] Rewrite `transfer_entropy.py`:
  - CLI: `--animal --condition --k {1,5,10} --binsize-ms 5`.
  - Replace `sio.toarray()` with explicit binning of spike trains from `data/processed/spiketrains_*.txt` (consistent with GLMCC inputs).
  - Output `results/te/te_<animal>_<condition>_k<k>.csv` (NxN matrix).
- [ ] Decide on `k` by AIC/BIC on a held-out animal (sketch in `src/py/te_choose_k.py`). Default kept at `k=10` if no better choice.
- [ ] Acceptance: 26 files at the chosen `k`; the file `180131_ongoing` (currently missing in `TE/`) is regenerated.

### 2.2 TE significance
- [ ] Run TE on the same 100 dither surrogates from §1.1 → `results/te/surrogate_te_<animal>_<condition>/<surr_id>.csv`.
- [ ] Edge-level p-values + BH-FDR → `results/te/validated_te_<animal>_<condition>.csv`.

### 2.3 Cross-check: same spike-train inputs feed GLMCC and TE
- [ ] Add `tests/test_inputs.py` asserting that the SpikeTrain object passed to GLMCC and to TE for a given animal × condition is byte-identical (same binning, same epoching).

---

## Sprint 3 — Region table, graph metrics and aggregation (week 4)

### 3.1 Region scheme — pre-register the 6-region set
- [ ] Lock the 6 macro-regions in `regions.csv`: **PVH, PH, ZI, VM-thalamus, DMH, ARH**.
- [ ] Also keep a 9-region column for the supplement (adds Anterior HT, Mammillary Complex, Ventromedial HT).
- [ ] Acceptance: per-region neuron count table (animal × region) generated as `results/region_yield.csv` and `figures/main/fig1_yield.png`.

### 3.2 Replace hand-coded region indices in `glmcc_analyses.R`
- [ ] All `regions_<animal>[seq(...)] = "..."` blocks → replaced by a single `region <- neurons %>% filter(animal == X) %>% pull(region6)` lookup.
- [ ] Acceptance: `grep -n 'regions_[0-9]' glmcc_analyses.R` returns nothing.

### 3.3 Per-animal graph metrics (both estimators)
- [ ] For each (animal, condition, estimator ∈ {glmcc, te}): node strength, clustering coefficient, local efficiency, hub score (z_i + participation coefficient).
- [ ] Output `results/<estimator>/metrics_<animal>_<condition>.csv` with `neuron_id, region, metric, value`.

### 3.4 Region-level aggregation — fix the count bias
- [ ] New script `src/r/aggregate_density.R`:
  - For each animal × condition × region-pair (r,s): density = (#significant edges between r and s) / (#possible pairs between r and s).
  - Output `results/<estimator>/density_<animal>_<condition>.csv`.
- [ ] Aggregate across animals by *averaging densities* (not summing counts). Output `results/<estimator>/density_meanCI_<day|night>_<ongoing|lightON>.csv`.

---

## Sprint 4 — Statistics, refit with mixed-effects (week 5)

### 4.1 First-level: neuron-as-observation models
- [ ] In `src/r/stats_first_level.R`, replace every `lm(metric ~ condition*region + (1|nodeid))` with:
  ```r
  lmer(metric ~ condition * region + (1 | animal/neuron),
       data = df, REML = TRUE)
  ```
  using `lmerTest` + Kenward–Roger.
- [ ] Outputs per metric × estimator: fixed-effect estimates, 95 % CIs, Type-III ANOVA, ICC, residual diagnostics.

### 4.2 Second-level: day × night × condition × region
- [ ] Model:
  ```r
  lmer(density ~ condition * region * day_night + (1 | animal),
       data = df_density)
  ```
- [ ] Pre-specified contrasts (single family, BH-FDR or Šidák):
  1. LightON − Ongoing | Day
  2. LightON − Ongoing | Night
  3. (LightON − Ongoing | Night) − (LightON − Ongoing | Day)
  4. Region-pair-specific tests for the two reported absences: PV-HT↔DM-HT, PV-HT↔VM-T, Post-HT↔Arc-HT.

### 4.3 Hub stability (quantitative replacement of alluvials)
- [ ] Per animal: Spearman ρ of hub score, Ongoing vs LightON. Aggregate with Wilcoxon-signed-rank and BCa bootstrap CI.
- [ ] Permutation null: shuffle condition labels within animal × neuron, n_perm = 10,000.
- [ ] Output: `figures/main/fig4_hub_stability.png` and CSV with ρ, p, CI per animal.

### 4.4 Convergent validity (GLMCC ↔ TE)
- [ ] Per animal × condition: Spearman ρ of node strength between GLMCC and TE. Same for clustering coefficient, local efficiency.
- [ ] Aggregate ρ across animals; report mean + CI. Flag region-pairs where GLMCC and TE disagree (sign or magnitude).
- [ ] Output: `figures/main/fig5_convergence.png`.

---

## Sprint 5 — Temporal structure of the LightON paradigm (week 6)

### 5.1 Adaptation check
- [ ] Split the 5 h LightON session into hour-1 vs hour-5. Recompute firing rate and graph metrics per half.
- [ ] If first vs last hour differ significantly, restrict main analyses to the stationary regime.

### 5.2 Early vs late evoked window
- [ ] Within each 10 s ON pulse, partition into early `[0, 500 ms]` and late `[500 ms, 10 s]`.
- [ ] Re-run §1–§4 on each window; report in supplement.

---

## Sprint 6 — Workflow engine + figures (week 7)

### 6.1 Snakemake DAG
- [ ] Write `Snakefile` keyed by wildcards `{animal}_{condition}_{estimator}`. Targets:
  - `data/processed/spiketrains_{animal}_{condition}.txt`
  - `data/surrogates/{animal}_{condition}/{surr_id}/cell_*.txt`
  - `results/{estimator}/adj_{animal}_{condition}.csv`
  - `results/{estimator}/validated_adj_{animal}_{condition}.csv`
  - `results/{estimator}/metrics_{animal}_{condition}.csv`
  - `results/{estimator}/density_{animal}_{condition}.csv`
- [ ] `snakemake -n` shows 13 × 2 × 2 = 52 leaf jobs covered. CI runs on a small subset.

### 6.2 Deterministic figure regeneration
- [ ] `src/r/make_figures.R` and `src/py/make_figures.py` produce every panel from `results/` with a fixed RNG seed.
- [ ] Retire all `*_v2.png`, `*_fixed.png`, `*_reduced.png` from `figures/` (move to `legacy/figures/`).

### 6.3 Five-figure layout
- [ ] **Fig 1** — Paradigm: stereotaxis + probe (Allen atlas), neuron yield heatmap (13 × 6), stimulus protocol + example raster.
- [ ] **Fig 2** — Firing rate + three graph metrics, Ongoing vs LightON, mixed-effect estimates + CIs.
- [ ] **Fig 3** — 2 × 2 (Day/Night × Ongoing/LightON) region-level GLMCC circos; FDR-significant edges marked; companion contrast heatmap.
- [ ] **Fig 4** — Hub stability: per-animal ρ, group summary, permutation null.
- [ ] **Fig 5** — Convergent validity GLMCC ↔ TE on node strength, CC, LE.
- [ ] **Supplement** — Per-animal full circos, neuron-count table, GLMCC sensitivity sweep, early/late evoked windows, TE k-selection.

---

## Sprint 7 — Manuscript and reporting (week 8)

### 7.1 Methods parameter table (anchor)
| Layer | Value |
|---|---|
| Bin width (GLMCC) | TBD from §1.4 |
| Delay window (GLMCC) | ±25 ms (default; sweep in supp) |
| Surrogate method | dither_spike_train, Δt = 0.1 ms, n = 100 |
| Edge threshold | BH-FDR q = 0.05 (per animal × condition) |
| TE history (k) | TBD from §2.1 |
| TE binsize | 5 ms |
| Region scheme | 6 macro-regions (pre-registered); 9 in supp |
| First-level model | lmer(metric ~ condition*region + (1\|animal/neuron)) |
| Second-level model | lmer(density ~ condition*region*day_night + (1\|animal)) |
| Multiple comparisons | BH-FDR within pre-specified contrast family |

### 7.2 Draft and circulate
- [ ] Methods (rewritten with the table above).
- [ ] Results section reorganised around the 5 main figures.
- [ ] Send to Storchi (Manchester) and Brown (Manchester) with the new mixed-effect estimates, CIs, and the GLMCC ↔ TE convergence figure.

### 7.3 Pre-submission checks
- [ ] Reproducibility audit: clone repo to a clean machine, run `snakemake all`, regenerate every figure with no manual step.
- [ ] Statistics audit: every p-value in the manuscript traces back to a single function call in `src/r/stats_*.R`.
- [ ] Data deposit plan: spike trains + region table + adjacency matrices → Zenodo / G-Node.

---

## Out of scope but worth tracking for the discussion
- Wavelength control (>500 nm, melanopsin-insensitive) — needs new recordings.
- Behavioural correlate (core temperature, corticosterone) — needs new cohort.
- Chronic Neuropixels follow-up for within-animal hub stability across cycles.
