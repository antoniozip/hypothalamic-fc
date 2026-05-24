# Implementation Plan — Hypothalamic FC Project Review

> Branch: `review-improvements` (forked from `main`)
> Created: 2026-05-24
> Based on: Project review of `tim_brown/`

---

## Phase 0 — Git setup & hygiene (one-time)

### 0.1 — Proper .gitignore
Add patterns and remove tracked violations:
```
**/.DS_Store
__pycache__/
*.py[oc]
.venv/
graphify-out/
*.RData
*.Rhistory
data/raw/*.mat        # large MAT files — Git LFS or exclude
*.mp4                 # video files
*.tiff                # large TIFFs
```

### 0.2 — Initial commit on main
Commit only source code, config, and small data files. Exclude:
- `data/raw/*.mat` (>100 MB each)
- `*.mp4`, `*.tiff`
- `graphify-out/` (generated)
- `.Rhistory`, `__pycache__/`, `.DS_Store`

### 0.3 — Branch
```bash
git add -A  # after .gitignore is correct
git commit -m "Initial commit: hypothalamic FC project"
git checkout -b review-improvements
```

**Verify:**
- `git status` is clean on `main`
- `review-improvements` branch exists
- `.DS_Store` and `__pycache__` are NOT tracked

---

## Phase 1 — CRITICAL fixes (statistical validity)

### 1.1 — Fix `src/r/stats_first_level.R`: add region to model

**Current (broken):**
```r
lmer(metric ~ condition + (1 | animal/neuron), data = all_data)
```

**Target:**
```r
lmer(metric ~ condition * region + (1 | animal/neuron), data = all_data)
```

**Steps:**
1. Update formula on line 68
2. Update the fallback model (lines 71-76) to also include region
3. Add `region` to the factor conversion (already done on line 58)
4. Add `condition:region` interaction terms to the output table
5. Re-run: `Rscript src/r/stats_first_level.R --estimator glmcc`
6. Diff old vs new `results/stats_first_level_glmcc.csv`
7. Update manuscript Methods section if estimates change materially

**Risk:** Including `condition * region` with 6+ regions and only 13 animals may cause convergence failures. Mitigation: implement `lmerTest::lmer(..., control = lmerControl(optimizer = "bobyqa"))` with a tryCatch fallback to `condition + region` (main effects only, no interaction).

**Files:** `src/r/stats_first_level.R`

---

### 1.2 — Replace heuristic p-values with surrogate-based validation

**Problem:** `validate_edges.py` lines 80-95 assign fabricated p-values.

**Solution A (preferred):** Merge `validate_ccg.py` logic into `validate_edges.py`:
1. Remove the heuristic `else` branch entirely
2. When `--surrogates` is not set, auto-detect if surrogate files exist; if not, raise an error
3. Add `--method {ccg,glmcc}` flag — CCG peak method for lightON (already working in validate_ccg.py), GLMCC for ongoing
4. Use `validate_ccg.py`'s shuffle-ISI surrogate generation (elephant library)

**Solution B (quick):** Rename heuristic outputs to `validated_adj_*_heuristic.csv` and add a metadata column flagging each edge as `empirical` vs `heuristic`. Update all downstream consumers to filter on this flag.

**Files:** `src/py/validate_edges.py`, `src/py/validate_ccg.py` (merge)

---

### 1.3 — Fix `run_glmcc.py` silent failure and path handling

**Steps:**
1. Replace hardcoded `~/Documents/GLMCC/` with environment variable `GLMCC_HOME` or config file lookup
2. Remove `shell=True` from subprocess call — use list form:
   ```python
   subprocess.run([sys.executable, str(glmcc_path), str(spike_path_abs), "exp", "GLM"])
   ```
3. On failure, write a metadata JSON alongside the placeholder CSV with `{"status": "glmcc_unavailable", "error": "..."}` 
4. Print clear warning to stderr: `"WARNING: GLMCC not available — wrote ZERO placeholder. Results WILL be invalid."`
5. Remove unused `spike_trains` list accumulation (lines 50-60) — replace with lightweight line count

**Files:** `src/py/run_glmcc.py`

---

## Phase 2 — HIGH priority (pipeline & reproducibility)

### 2.1 — Complete the Snakefile DAG

Add these missing rules:

```python
rule extract_spike_trains:
    """MAT → processed spike train txt"""
    input:  "data/raw/DATA{animal}_{condition}.mat"
    output: "data/processed/spiketrains_DATA{animal}_{condition}.txt"
    shell:  "python src/py/extract_spike_trains.py --animal {wildcards.animal} --condition {wildcards.condition}"

rule generate_surrogates:
    """Generate 100 shuffle-ISI surrogates"""
    input:  "data/processed/spiketrains_DATA{animal}_{condition}.txt"
    output: expand("data/surrogates/{animal}_{condition}/{sid}/cell_{i}.txt", sid=range(100), i=range(64))
    shell:  "python src/py/generate_surrogates.py --animal {wildcards.animal} --condition {wildcards.condition} --n 100"

rule run_glmcc_adjacency:
    input:  "data/processed/spiketrains_DATA{animal}_{condition}.txt"
    output: "results/glmcc/adj_{animal}_{condition}.csv"
    shell:  "python src/py/run_glmcc.py --animal {wildcards.animal} --condition {wildcards.condition}"

rule run_te_adjacency:
    input:  "data/processed/spiketrains_DATA{animal}_{condition}.txt"
    output: "results/te/te_{animal}_{condition}_k10.csv"
    shell:  "python src/py/run_transfer_entropy.py --animal {wildcards.animal} --condition {wildcards.condition} --k 10"
```

Also add `rule all` targets for figures 1 and 4.

**Files:** `Snakefile`

---

### 2.2 — Containerize GLMCC

**Steps:**
1. Identify exact GLMCC commit (check `src/py/GLMCC_VERSION.txt`)
2. Write `Dockerfile`:
   ```dockerfile
   FROM python:3.12-slim
   RUN apt-get update && apt-get install -y gcc g++
   RUN git clone https://github.com/.../GLMCC.git /opt/GLMCC
   WORKDIR /opt/GLMCC
   RUN git checkout <COMMIT_HASH>
   RUN pip install -r requirements.txt
   COPY entrypoint.sh /entrypoint.sh
   ENTRYPOINT ["/entrypoint.sh"]
   ```
3. Write `scripts/run_glmcc_container.sh` wrapper
4. Update `run_glmcc.py` to detect container vs local installation
5. Add to Snakefile: `container: "docker://glmcc:latest"` on the glmcc rule

**Files:** `Dockerfile`, `scripts/run_glmcc_container.sh`, `src/py/run_glmcc.py`

---

### 2.3 — Address TE zero-output on sparse data

**Investigation steps:**
1. Run TE on a dense synthetic spike train (Poisson 10 Hz) → should produce non-zero
2. Run TE on real data with bin sizes {1, 5, 10, 20, 50} ms → find threshold where signal emerges
3. Test alternative implementations: pyinform vs numpy vs JIDT (Java) via subprocess
4. Document minimum firing rate / bin size requirements
5. If pyinform results are used, add a supplementary note explaining the discrepancy

**Files:** `src/py/te_choose_k.py` (extend), new `src/py/te_validation.py`

---

## Phase 3 — MEDIUM priority (code quality)

### 3.1 — Fix BH-FDR test count

In both `validate_edges.py` (line 98) and `validate_ccg.py` (line 174):
```python
# Change:
n_tests = len(p_flat)
# To:
n_tests = n_units * (n_units - 1)  # exclude diagonal
```
Also exclude diagonal p-values from the flat array before sorting.

**Files:** `src/py/validate_edges.py`, `src/py/validate_ccg.py`

### 3.2 — Repo hygiene sweep

1. Remove tracked `__pycache__/` and `.DS_Store`:
   ```bash
   git rm --cached -r src/py/__pycache__/
   git rm --cached src/.DS_Store results/.DS_Store
   ```
2. Populate `README.md` with:
   - Project overview (1 paragraph from manuscript abstract)
   - Setup instructions (renv + uv)
   - Quick start (snakemake -n)
   - Parameter table from TODO.md Sprint 7

**Files:** `README.md`, `.gitignore`, tracked cache files

### 3.3 — Remove dead code from stats_first_level.R

Line 46 loads `day_night` but it's never used. Either:
- Remove the column assignment (if second-level handles it), OR
- Add a supplementary model: `lmer(metric ~ condition * day_night + (1|animal/neuron))`

**Files:** `src/r/stats_first_level.R`

### 3.4 — Add region data to metrics CSVs

`stats_first_level.R` reads `metrics_*.csv` expecting a `region` column. Verify that `compute_graph_metrics.R` actually writes a `region` column to its output. If not, merge region info from `neurons.csv`.

**Files:** `src/r/compute_graph_metrics.R`

---

## Phase 4 — LOW priority (polish & graphify insights)

### 4.1 — Wire missing figures into Snakefile

```python
rule fig1_yield:
    input:  "data/processed/neurons.csv"
    output: "figures/main/fig1_yield.png"
    shell:  "Rscript src/r/make_figures.R --figure fig1"

rule fig4_hub_stability:
    input:  "results/hub_stability.csv"
    output: "figures/main/fig4_hub_stability.png"
    shell:  "Rscript src/r/hub_stability.R"
```

### 4.2 — Clean up root-level large files

Move to `data/raw/` (as TODO.md Sprint 0.4 planned):
- `DATA*.mat` → `data/raw/`
- `*.mp4` → `data/raw/video/`
- `*.tiff` → `figures/source/`

Update all scripts that reference root-level paths.

### 4.3 — Act on graphify isolated nodes

Review 289 isolated nodes for:
- Dead code candidates (delete or move to legacy)
- Undocumented internal functions (add docstrings)
- Missing cross-references (add imports/calls)

### 4.4 — Pin Python dependencies

Replace `>=` with `==` in `pyproject.toml`:
```toml
dependencies = [
    "elephant==1.2.0",
    "neo==0.14.4",
    ...
]
```
Or add `uv.lock` with exact hashes.

---

## Execution order

| Step | Phase | Effort | Depends on |
|------|-------|--------|------------|
| 0.1-0.3 | Git setup | 30 min | nothing |
| 1.1 | Fix stats model | 1 hr | 0.3 |
| 1.2 | Fix p-values | 2 hr | 0.3 |
| 1.3 | Fix run_glmcc | 45 min | 0.3 |
| 2.1 | Snakefile DAG | 1.5 hr | 1.3 |
| 3.2 | Repo hygiene | 30 min | 0.3 |
| 3.1 | Fix FDR count | 15 min | 1.2 |
| 3.3 | Dead code removal | 15 min | 1.1 |
| 3.4 | Region column check | 30 min | 1.1 |
| 2.2 | GLMCC container | 2 hr | 1.3 |
| 2.3 | TE investigation | 2 hr | — |
| 4.1 | Wire figures | 30 min | 2.1 |
| 4.2 | Move large files | 1 hr | 2.1 |
| 4.3-4.4 | Polish | 1 hr | — |

---

## Success criteria

- [ ] `snakemake -n` shows complete DAG with no missing inputs
- [ ] `Rscript src/r/stats_first_level.R --estimator glmcc` produces `condition:region` interaction terms
- [ ] `python src/py/validate_edges.py --surrogates --animal 171019 --condition lightON` runs without heuristic fallback
- [ ] All `*.csv` in `results/glmcc/validated_adj_*` are surrogate-validated (no heuristic edges)
- [ ] `grep -r "heuristic\|placeholder\|ZERO" results/` returns nothing
- [ ] `.DS_Store` and `__pycache__` are not tracked by git
- [ ] `README.md` is non-empty and contains setup instructions
- [ ] GLMCC runs via `docker run` or `GLMCC_HOME` env var
