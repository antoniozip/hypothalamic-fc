# Hypothalamic Functional Connectivity Project

Light-evoked reorganization of functional connectivity across the mouse hypothalamus.

## Overview

Recordings from 682 neurons across 7 hypothalamic nuclei in 13 freely behaving mice
during a 5-hour light-ON paradigm. Functional connectivity is inferred using two
complementary directed estimators:

- **GLMCC** — Generalized Linear Model Cross-Correlation (Kobayashi et al. 2019)
- **Transfer Entropy** — Information-theoretic directed connectivity (Schreiber 2000)

Edges are validated against shuffle-ISI surrogates with Benjamini-Hochberg FDR
correction (q = 0.05). Graph-theoretic metrics (node strength, clustering coefficient,
local efficiency, hub score) are fitted with linear mixed-effects models.

## Quick Start

```bash
# Dry run — verify pipeline
snakemake -n

# Full pipeline (4 cores)
snakemake all -j4

# Single animal, single condition
python src/py/validate_edges.py --animal 171019 --condition lightON --method ccg
Rscript src/r/stats_first_level.R --estimator glmcc
```

## Setup

### R dependencies
```r
install.packages("renv")
renv::restore()
```

### Python dependencies
```bash
uv sync
# or: pip install -r <(uv export --no-hashes)
```

### GLMCC (external dependency)
Set `GLMCC_HOME` environment variable pointing to the GLMCC installation:
```bash
export GLMCC_HOME=~/Documents/GLMCC
```

## Repository Structure

```
.
├── data/
│   ├── raw/               # Raw .mat, cell files (not in git)
│   ├── processed/         # Extracted spike trains, neurons.csv
│   └── surrogates/        # Shuffle-ISI surrogate spike trains
├── src/
│   ├── py/                # Python pipeline scripts (14 files)
│   └── r/                 # R analysis scripts (6 files)
├── results/
│   ├── glmcc/             # GLMCC adjacencies, metrics, densities
│   └── te/                # Transfer entropy matrices
├── figures/main/          # Manuscript figures 1-5
├── legacy/                # Archived Pearson pipeline, old scripts
├── Snakefile              # Complete workflow DAG
└── manuscript_jneurosci.md # Manuscript draft
```

## Key Parameters

| Layer | Value |
|---|---|
| Bin width (GLMCC) | TBD from sensitivity sweep |
| Delay window (GLMCC) | ±25 ms |
| Surrogate method | shuffle-ISI, Δt = 0.1 ms, n = 100 |
| Edge threshold | BH-FDR q = 0.05 |
| TE history (k) | 10 |
| TE binsize | 5 ms |
| Region scheme | 6 macro-regions |
| First-level model | lmer(metric ~ condition * region + (1\|animal/neuron)) |
| Second-level model | lmer(density ~ condition × day_night + (1\|animal)) |

## Animals

| ID | Phase | Neurons |
|----|-------|---------|
| 171019 | Day | 126 |
| 171207 | Night | 24 |
| 171208 | Night | 32 |
| ... | ... | ... |

Full table: `data/processed/neurons.csv`

## References

- Kobayashi et al. (2019) — GLMCC method
- Schreiber (2000) — Transfer entropy
- Saper & Lowell (2014) — Hypothalamus review
- Oh et al. (2014) — Mouse connectome
- Swanson (2000) — Hypothalamus anatomy
