#!/usr/bin/env python3
"""Tier 3 — Motif Analysis + Distance Decay (temporal dynamics deferred).

Temporal dynamics requires recomputing full CCG adjacency matrices for
early/late windows, which is CPU-intensive (O(n^2 × spike_count) per animal).
Deferred to batch execution with more compute.

Output:
  figures/main/figT3_motif_analysis.png
  figures/main/figT3_distance_decay.png
  results/glmcc/motifs.csv
"""

from pathlib import Path
from io import StringIO
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FIGURES_DIR = PROJECT_ROOT / "figures" / "main"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

ANIMALS = [
    "day1", "night1", "night2", "night3", "night4", "night5",
    "day2", "night6", "night7", "day3", "day4", "day5", "day6",
]
DAYNIGHT = {
    "day1": "D", "night1": "N", "night2": "N", "night3": "N",
    "night4": "N", "night5": "N", "day2": "D", "night6": "N",
    "night7": "N", "day3": "D", "day4": "D", "day5": "D",
    "day6": "D",
}


def load_adjacency_binary(animal):
    adj_path = PROJECT_ROOT / "results" / "glmcc" / f"validated_adj_{animal}_lightON.csv"
    if not adj_path.exists():
        return None
    raw = adj_path.read_text()
    if raw.startswith("\ufeff"):
        raw = raw[1:]
    adj = np.loadtxt(StringIO(raw), delimiter=",", dtype=float)
    return (np.abs(adj) > 0).astype(int)


# ═══════════════════════════════════════════════════════════════════
# 9. MOTIF ANALYSIS
# ═══════════════════════════════════════════════════════════════════
print("=== 9. Motif Analysis ===\n")

N_NULL_MOTIF = 50
motif_results = []

for animal in ANIMALS:
    adj_bin = load_adjacency_binary(animal)
    if adj_bin is None or adj_bin.shape[0] < 3:
        print(f"  SKIP {animal}")
        continue

    n = adj_bin.shape[0]
    np.fill_diagonal(adj_bin, 0)
    rng = np.random.RandomState(42)

    # Count connected triplets in observed
    n_sample = min(10000, n * (n - 1) * (n - 2) // 6)
    connected_obs = 0
    for _ in range(n_sample):
        i, j, k = rng.choice(n, 3, replace=False)
        sub = adj_bin[np.ix_([i, j, k], [i, j, k])]
        if int(sub.sum()) >= 2:
            connected_obs += 1
    frac_obs = connected_obs / n_sample

    # Null: shuffle edges preserving in/out degree
    in_deg = adj_bin.sum(axis=0).astype(int)
    out_deg = adj_bin.sum(axis=1).astype(int)
    null_fracs = []

    for _ in range(N_NULL_MOTIF):
        rand_adj = np.zeros_like(adj_bin)
        out_stubs = np.repeat(np.arange(n), out_deg)
        in_stubs = np.repeat(np.arange(n), in_deg)
        rng.shuffle(out_stubs)
        rng.shuffle(in_stubs)
        for si in range(min(len(out_stubs), len(in_stubs))):
            oi, ii = out_stubs[si], in_stubs[si]
            if oi != ii:
                rand_adj[oi, ii] = 1

        conn_null = 0
        for _ in range(n_sample // 10):
            i, j, k = rng.choice(n, 3, replace=False)
            sub = rand_adj[np.ix_([i, j, k], [i, j, k])]
            if int(sub.sum()) >= 2:
                conn_null += 1
        null_fracs.append(conn_null / max(n_sample // 5, 1))

    null_mean = np.mean(null_fracs)
    null_std = np.std(null_fracs)
    z_score = (frac_obs - null_mean) / null_std if null_std > 0 else 0

    # Also count fully-connected triplets (3 edges = triangle)
    triangles_obs = 0
    for _ in range(n_sample):
        i, j, k = rng.choice(n, 3, replace=False)
        sub = adj_bin[np.ix_([i, j, k], [i, j, k])]
        if int(sub.sum()) == 6:  # all 6 possible directed edges
            triangles_obs += 1
    frac_tri = triangles_obs / n_sample

    motif_results.append({
        "animal": animal, "day_night": DAYNIGHT.get(animal, "?"),
        "n_neurons": n,
        "frac_connected": frac_obs, "frac_triangles": frac_tri,
        "null_mean": null_mean, "null_std": null_std,
        "z_score_connected": z_score,
    })
    print(f"  {animal}: connected={frac_obs:.3f}, triangles={frac_tri:.3f}, z={z_score:.2f}")

mdf = pd.DataFrame(motif_results)
mdf.to_csv(PROJECT_ROOT / "results" / "glmcc" / "motifs.csv", index=False)
print(f"\n  Mean connected frac: {mdf['frac_connected'].mean():.3f} ± {mdf['frac_connected'].std():.3f}")
print(f"  Mean Z: {mdf['z_score_connected'].mean():.2f}")

# --- Motif figure ---
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle("T9: Motif Analysis — 3-Node Connectivity Patterns",
             fontsize=14, fontweight="bold")

# A: Connected triplet fraction per animal
ax = axes[0, 0]
mdf_sorted = mdf.sort_values("frac_connected")
colors = ["#2166ac" if DAYNIGHT.get(a) == "D" else "#b2182b" for a in mdf_sorted["animal"]]
ax.barh(range(len(mdf_sorted)), mdf_sorted["frac_connected"], color=colors, alpha=0.85)
ax.set_yticks(range(len(mdf_sorted)))
ax.set_yticklabels(mdf_sorted["animal"], fontsize=8)
ax.set_xlabel("Fraction of Triplets with >=2 edges")
ax.set_title(f"A: 3-Node Connectivity (mean={mdf['frac_connected'].mean():.3f})")

# B: Z-scores vs null
ax = axes[0, 1]
z_vals = mdf_sorted["z_score_connected"].reindex(mdf_sorted.index)
bar_colors = ["#2166ac" if z > 0 else "#b2182b" for z in z_vals]
ax.barh(range(len(mdf_sorted)), z_vals, color=bar_colors, alpha=0.85)
ax.set_yticks(range(len(mdf_sorted)))
ax.set_yticklabels(mdf_sorted["animal"], fontsize=8)
ax.axvline(0, color="black", linestyle="--", linewidth=0.7)
ax.set_xlabel("Z-score vs Degree-Preserving Null")
ax.set_title(f"B: Motif Enrichment (mean Z={mdf['z_score_connected'].mean():.2f})")

# C: Connected vs Triangle fraction
ax = axes[1, 0]
ax.scatter(mdf["frac_connected"], mdf["frac_triangles"], c=colors, s=60, alpha=0.8)
ax.set_xlabel("Fraction Connected (>=2 edges)")
ax.set_ylabel("Fraction Fully Connected (6 edges)")
ax.set_title("C: Connected vs Fully-Connected Triplets")
if len(mdf) > 3:
    rho, p = stats.spearmanr(mdf["frac_connected"], mdf["frac_triangles"])
    ax.text(0.05, 0.95, f"rho={rho:.3f}", transform=ax.transAxes, fontsize=9)

# D: Motif diagram explanation (text panel)
ax = axes[1, 1]
ax.axis("off")
explanation = (
    "Motif Analysis Explanation\n"
    "==========================\n\n"
    "Counts 3-node subgraph patterns in\n"
    "the directed adjacency matrix.\n\n"
    f"Connected triplets (>=2 edges):\n"
    f"  Observed: {mdf['frac_connected'].mean():.3f}\n"
    f"  Null mean: {mdf['null_mean'].mean():.3f}\n\n"
    f"Mean Z = {mdf['z_score_connected'].mean():.2f}\n"
    f"(Z>1.96 = significant enrichment)\n\n"
    f"Triangles (all 6 edges):\n"
    f"  Mean: {mdf['frac_triangles'].mean():.3f}\n\n"
    "High connectivity means most\n"
    "triplets are connected by chance\n"
    "due to dense overall connectivity."
)
ax.text(0.1, 0.9, explanation, transform=ax.transAxes, fontsize=9,
        va="top", family="monospace")

plt.tight_layout()
fig.savefig(FIGURES_DIR / "figT3_motif_analysis.png", dpi=150, bbox_inches="tight")
print(f"  Saved: figT3_motif_analysis.png\n")


# ═══════════════════════════════════════════════════════════════════
# 10. DISTANCE-DEPENDENT CONNECTIVITY
# ═══════════════════════════════════════════════════════════════════
print("=== 10. Distance-Dependent Connectivity ===\n")

# A4x64-Poly2 probe geometry (approximate):
# 4 shanks, 16 channels per shank (2 cols x 8 rows)
CH_PER_SHANK = 16
N_SHANKS = 4
COLS = 2
ROWS = 8
COL_SPACING = 200  # um
ROW_SPACING = 160  # um
SHANK_SPACING = 250  # um

positions = []
for shank in range(N_SHANKS):
    sx = shank * SHANK_SPACING
    for row in range(ROWS):
        for col in range(COLS):
            positions.append((sx + col * COL_SPACING, row * ROW_SPACING))
positions = np.array(positions)
n_ch = len(positions)

# Pairwise distance matrix
dists = np.zeros((n_ch, n_ch))
for i in range(n_ch):
    for j in range(n_ch):
        dists[i, j] = np.sqrt(np.sum((positions[i] - positions[j]) ** 2))

dist_results = []

for animal in ANIMALS:
    adj_bin = load_adjacency_binary(animal)
    if adj_bin is None:
        continue

    n_u = min(adj_bin.shape[0], n_ch)
    adj_sub = adj_bin[:n_u, :n_u]
    np.fill_diagonal(adj_sub, 0)
    d_sub = dists[:n_u, :n_u]

    d_flat = d_sub[np.triu_indices(n_u, 1)]
    c_flat = adj_sub[np.triu_indices(n_u, 1)]

    if len(d_flat) > 10:
        rho, p = stats.spearmanr(d_flat, c_flat)
    else:
        rho, p = np.nan, np.nan

    dist_results.append({
        "animal": animal, "day_night": DAYNIGHT.get(animal, "?"),
        "n_neurons": n_u, "spearman_rho": rho, "p_value": p,
        "mean_dist_um": d_flat.mean(),
    })
    print(f"  {animal}: rho(dist,conn)={rho:.3f}, p={p:.4f}")

ddf = pd.DataFrame(dist_results)
print(f"\n  Mean rho: {ddf['spearman_rho'].mean():.3f} ± {ddf['spearman_rho'].std():.3f}")
n_sig = (ddf['p_value'] < 0.05).sum()
n_total = ddf['p_value'].notna().sum()
print(f"  Significant (p<0.05): {n_sig}/{n_total}")

# --- Distance decay figure ---
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle("T10: Distance-Dependent Connectivity (A4x64-Poly2 probe)",
             fontsize=14, fontweight="bold")

# A: Scatter + binned curve for largest animal (day1)
ax = axes[0, 0]
animal_show = "day1"
adj_show = load_adjacency_binary(animal_show)
if adj_show is not None:
    n_s = min(adj_show.shape[0], n_ch)
    d_flat = dists[:n_s, :n_s][np.triu_indices(n_s, 1)]
    c_flat = adj_show[:n_s, :n_s][np.triu_indices(n_s, 1)]

    # Subsample for scatter
    n_pts = min(3000, len(d_flat))
    idx = np.random.RandomState(0).choice(len(d_flat), n_pts, replace=False)
    ax.scatter(d_flat[idx], c_flat[idx] + np.random.normal(0, 0.01, n_pts),
               alpha=0.3, s=4, c="#444444")

    # Binned
    bins = np.arange(0, d_flat.max() + 50, 100)
    bc, bp = [], []
    for k in range(len(bins) - 1):
        m = (d_flat >= bins[k]) & (d_flat < bins[k + 1])
        if m.sum() >= 5:
            bc.append((bins[k] + bins[k + 1]) / 2)
            bp.append(c_flat[m].mean())
    ax.plot(bc, bp, "r-o", linewidth=2, markersize=4, label="Binned mean")
    rho_val = ddf[ddf["animal"] == animal_show]["spearman_rho"].values[0]
    ax.set_xlabel("Inter-Channel Distance (um)")
    ax.set_ylabel("Connection Probability")
    ax.set_title(f"A: Distance Decay ({animal_show}, rho={rho_val:.3f})")
    ax.legend(fontsize=8)

# B: Per-animal rho
ax = axes[0, 1]
ddf_v = ddf.dropna(subset=["spearman_rho"]).sort_values("spearman_rho")
colors_b = ["#2166ac" if DAYNIGHT.get(a) == "D" else "#b2182b" for a in ddf_v["animal"]]
ax.barh(range(len(ddf_v)), ddf_v["spearman_rho"], color=colors_b, alpha=0.85)
ax.set_yticks(range(len(ddf_v)))
ax.set_yticklabels(ddf_v["animal"], fontsize=8)
ax.axvline(0, color="black", linestyle="--", linewidth=0.7)
ax.set_xlabel("Spearman rho (Distance vs Connection)")
ax.set_title(f"B: Distance-Connectivity Correlation\n(mean rho={ddf['spearman_rho'].mean():.3f})")

# C: Same-shank vs cross-shank
ax = axes[1, 0]
same_conns_all = []
cross_conns_all = []
for animal in ANIMALS:
    adj = load_adjacency_binary(animal)
    if adj is None:
        continue
    n_u = min(adj.shape[0], n_ch)
    for shank in range(N_SHANKS):
        start = shank * CH_PER_SHANK
        end = start + CH_PER_SHANK
        for i in range(start, min(end, n_u)):
            for j in range(n_u):
                if i == j:
                    continue
                if start <= j < end:
                    same_conns_all.append(adj[i, j])
                else:
                    cross_conns_all.append(adj[i, j])

same_mean = np.mean(same_conns_all) if same_conns_all else 0
cross_mean = np.mean(cross_conns_all) if cross_conns_all else 0
ax.bar(["Same Shank", "Cross-Shank"], [same_mean, cross_mean],
       color=["#2166ac", "#b2182b"], alpha=0.85)
ax.set_ylabel("Mean Connection Probability")
ax.set_title(f"C: Within vs Cross-Shank\n(Same={same_mean:.3f}, Cross={cross_mean:.3f})")

# D: Distance distribution
ax = axes[1, 1]
all_d = []
for animal in ANIMALS:
    adj = load_adjacency_binary(animal)
    if adj is None:
        continue
    n_u = min(adj.shape[0], n_ch)
    all_d.extend(dists[:n_u, :n_u][np.triu_indices(n_u, 1)])
ax.hist(all_d, bins=30, color="#6baed6", edgecolor="white", alpha=0.85)
ax.set_xlabel("Inter-Channel Distance (um)")
ax.set_ylabel("Count")
ax.set_title(f"D: Pairwise Distance Distribution\n(n={len(all_d):,} channel pairs)")

plt.tight_layout()
fig.savefig(FIGURES_DIR / "figT3_distance_decay.png", dpi=150, bbox_inches="tight")
print(f"  Saved: figT3_distance_decay.png")

print("\n=== Tier 3 (motifs + distance) complete ===")
print("Note: Temporal dynamics deferred — requires recomputing CCG per window (CPU-heavy).")
