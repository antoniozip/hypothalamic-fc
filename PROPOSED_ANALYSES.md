# Proposed New Analyses — Hypothalamic FC Project

Based on the regenerated pipeline and unexploited data.

---

## Tier 1 — High Impact, Low Effort (use existing CSVs)

### 1. Rich-Club Analysis
**What:** Do high-degree nodes preferentially connect to each other?
**Why:** The manuscript references van den Heuvel & Sporns 2011 (Rich Club) in the methods but never computed it. A significant rich club would mean the hypothalamic network has a core of densely interconnected hubs — this is a hallmark of integrative networks and a stronger claim than "PVH and DMH are hubs."
**Data needed:** `validated_adj_*.csv` (already have)
**Output:** Rich-club coefficient curve per animal, normalized against 1000 degree-preserving random networks. Group-level rich-club significance.
**Script:** `src/r/rich_club.R`

### 2. Edge Consistency / Core Connectome
**What:** Which specific neuron pairs are consistently connected across animals?
**Why:** Identifies a "canonical hypothalamic connectome" — edges present in most animals are likely anatomical, while variable edges may be state-dependent. Directly addresses the question "is hypothalamic connectivity stereotyped or variable?"
**Data needed:** `validated_adj_*.csv` + region mapping from `neurons.csv`
**Approach:** Binarize each adjacency → for each region pair, compute fraction of animals where an edge exists. Rank edges by consistency.
**Output:** Consistency matrix (region × region), per-region-pair consistency score, top-10 most consistent edges.
**Script:** `src/r/edge_consistency.R`

### 3. Within vs Between-Region Connectivity Decomposition
**What:** What fraction of total connectivity is within-region vs between-region?
**Why:** Differentiates local processing (within-region) from long-range integration (between-region). The high connector-hub proportion (94.4%) suggests between-region connectivity dominates — let's quantify it.
**Data needed:** `validated_adj_*.csv` + `neurons.csv` region labels
**Output:** Per-animal within/between ratio, region-pair density matrix normalized by neuron count, group comparison (Day vs Night).
**Script:** modify existing `aggregate_density.R` or new `src/r/decompose_regions.R`

### 4. Centrality Convergence
**What:** Do strength, betweenness, eigenvector, and PageRank centrality agree on hub identity?
**Why:** If all centrality measures converge on the same neurons, hub status is robust. If they diverge, different neurons play different hub roles (e.g., connector vs integrator).
**Data needed:** `validated_adj_*.csv`
**Output:** Spearman correlation matrix of 4 centrality measures per animal. Per-neuron centrality profile. Group-level convergence score.
**Script:** `src/r/centrality_convergence.R`

### 5. Firing Rate × Connectivity Coupling
**What:** Does a neuron's firing rate predict its node strength or hub score?
**Why:** Distinguishes "rate-driven" connectivity (neurons appear connected simply because they fire more) from genuine coupling. If FR and strength are uncorrelated, the detected edges reflect genuine temporal coupling rather than rate artifacts.
**Data needed:** `validated_adj_*.csv` + `.mat` spike times (for FR computation)
**Output:** Per-animal FR × strength scatter with Spearman ρ. Group-level summary. Split by E/I classification.
**Script:** `src/py/fr_connectivity_coupling.py`

---

## Tier 2 — High Impact, Medium Effort

### 6. Surrogate Null Models for Graph Metrics
**What:** Are the observed modularity, clustering, and efficiency higher than expected by chance?
**Why:** Edge validation tests individual edges — but graph metrics can be significant even with random edge patterns. Comparing observed graph metrics against degree-preserving random networks tells us whether the network topology itself is non-random.
**Data needed:** `validated_adj_*.csv` + igraph rewiring functions
**Approach:** For each animal, generate 1000 degree-preserving random networks (rewire(..., keeping_degseq=TRUE)). Compute modularity, clustering, path length for each. Compare observed to null distribution.
**Output:** Z-scores per metric per animal. Group-level significance.
**Script:** `src/r/null_models.R`

### 7. Global Efficiency & Small-World Propensity
**What:** Is the hypothalamic network a small-world network?
**Why:** Small-world networks balance local clustering with global integration — perfect for a structure that needs both specialized nuclei and whole-hypothalamus coordination. We have local efficiency; adding global efficiency and small-world indices (σ, ω) completes the picture.
**Data needed:** `validated_adj_*.csv`
**Output:** Sigma and omega per animal. Comparison against equivalent random + lattice networks. Group summary.
**Script:** `src/r/small_world.R`

---

## Tier 3 — Exploratory / Requires Additional Data

### 8. Temporal Dynamics: Early vs Late LightON
**What:** Does the network reorganize over the 5-hour recording?
**Why:** The manuscript mentions adaptation (Sprint 5) and reports firing rate decreases in 7/13 animals. Does connectivity also change? Split lightON into hour-1 vs hour-5, recompute adjacencies, compare.
**Data needed:** Raw spike times from `.mat` files, split by time window
**Output:** Per-animal early-vs-late adjacency comparison, Jaccard similarity of edge sets, metric change.
**Script:** `src/py/temporal_dynamics.py`

### 9. Motif Analysis
**What:** Are specific 3-node circuit motifs over-represented?
**Why:** Network motifs are the building blocks of computation. Over-represented feedforward loops suggest hierarchical processing; over-represented feedback loops suggest recurrent dynamics. This is the most mechanistic analysis possible from the connectivity data.
**Data needed:** Binarized `validated_adj_*.csv`
**Approach:** Count all 13 possible 3-node motifs in each adjacency. Compare against 1000 degree-preserving random networks (mfinder or in-house implementation). Z-score per motif.
**Output:** Motif significance profile per animal. Group-level over-represented motifs.
**Script:** `src/r/motif_analysis.R` (or Python with graph-tool)

### 10. Distance-Dependent Connectivity
**What:** Does connection probability decay with physical distance between recording sites?
**Why:** If connectivity is distance-independent, the network is truly functional (not just volume conduction or proximity artifact). If it decays with distance, spatial organization matters. Requires probe geometry (channel coordinates).
**Data needed:** Probe layout specs (A4x64-Poly2-7mm-23s-200-160) + `validated_adj_*.csv`
**Output:** Distance-decay curve per animal. Exponential fit (λ = space constant). E vs I comparison.
**Script:** `src/py/distance_decay.py`

---

## Recommended Execution Order

| # | Analysis | Time | Dependencies |
|---|----------|------|-------------|
| 1 | Rich-club | 1h | validated_adj only |
| 2 | Edge consistency | 1h | validated_adj + neurons.csv |
| 3 | Within/between decomposition | 30min | Same |
| 4 | Centrality convergence | 45min | validated_adj only |
| 5 | FR × connectivity | 1h | validated_adj + .mat |
| 6 | Surrogate null models | 2h | validated_adj + igraph |
| 7 | Small-world | 1h | validated_adj + igraph |
| 8 | Temporal early/late | 3h | Raw .mat spike times |
| 9 | Motif analysis | 3h | validated_adj + null model |
| 10 | Distance decay | 2h | Probe geometry spec |

## Papers These Analyses Would Support

- Rich-club → "Core-periphery organization of the hypothalamic functional connectome"
- Edge consistency → "A stereotyped core connectome across individual mice"
- Motif analysis → "Canonical circuit motifs in hypothalamic processing"
- Within/between → "Local processing vs long-range integration in the hypothalamus"
- FR × connectivity → "Dissociating rate-driven from coupling-driven connectivity"
- Small-world → "Small-world topology of light-evoked hypothalamic networks"
