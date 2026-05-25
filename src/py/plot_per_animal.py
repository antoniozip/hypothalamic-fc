#!/usr/bin/env python3
"""Per-animal validated adjacency heatmaps.

Reads CCG-validated adjacency matrices and generates a multi-panel figure
showing connectivity patterns for each animal × condition.

Output: figures/main/figS1_per_animal_adjacencies.png
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "glmcc"
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

NEURON_COUNTS = {
    "day1": 126, "night1": 24, "night2": 32, "night3": 33,
    "night4": 35, "night5": 24, "day2": 31, "night6": 54,
    "night7": 79, "day3": 47, "day4": 61, "day5": 71,
    "day6": 56,
}


def load_adjacency(animal: str, condition: str) -> np.ndarray | None:
    """Load validated adjacency, trying multiple possible paths."""
    candidates = [
        RESULTS_DIR / f"validated_adj_{animal}_{condition}.csv",
        RESULTS_DIR / f"validated_adj_ccg_{animal}_{condition}.csv",
        RESULTS_DIR / f"validated_adj_{animal}_ongoing_bis.csv",
    ]
    for p in candidates:
        if p.exists():
            raw = p.read_text()
            if raw.startswith("\ufeff"):
                raw = raw[1:]
            from io import StringIO
            arr = np.loadtxt(StringIO(raw), delimiter=",", dtype=float, ndmin=2)
            return arr
    return None


def compute_edge_stats(adj: np.ndarray) -> dict:
    """Compute edge survival stats."""
    n = adj.shape[0]
    n_possible = n * (n - 1)
    n_sig = int((adj != 0).sum())  # off-diagonal non-zero
    # Subtract diagonal if any
    n_sig -= int(np.diag(adj).astype(bool).sum())
    survival = n_sig / n_possible if n_possible > 0 else 0
    return {"n_units": n, "n_significant": n_sig, "survival": survival}


def plot_all_heatmaps():
    """Generate multi-panel per-animal adjacency heatmaps."""
    # Collect all available adjacencies
    lightON_available = []
    ongoing_available = []
    for animal in ANIMALS:
        adj_l = load_adjacency(animal, "lightON")
        if adj_l is not None:
            lightON_available.append((animal, adj_l))
        adj_o = load_adjacency(animal, "ongoing_bis")
        if adj_o is None:
            adj_o = load_adjacency(animal, "ongoing")
        if adj_o is not None:
            ongoing_available.append((animal, adj_o))

    n_light = len(lightON_available)
    n_ongoing = len(ongoing_available)
    total_panels = n_light + n_ongoing

    # Layout: 5 columns, enough rows for all panels
    ncols = 5
    nrows_light = (n_light + ncols - 1) // ncols
    nrows_ongoing = (n_ongoing + ncols - 1) // ncols
    nrows_total = nrows_light + nrows_ongoing + 1  # +1 for separator row

    # Find global color scale (use 95th percentile to avoid outlier saturation)
    all_vals = []
    for _, adj in lightON_available + ongoing_available:
        all_vals.extend(adj[adj != 0].flatten())
    if all_vals:
        vmax = np.percentile(np.abs(all_vals), 95)
    else:
        vmax = 1.0

    fig = plt.figure(figsize=(20, 4.5 * nrows_total))
    fig.suptitle(
        "Per-Animal Validated Adjacency Matrices (GLMCC, CCG + shuffle-ISI surrogates, BH-FDR q=0.05)",
        fontsize=16, fontweight="bold", y=0.995,
    )

    panel_idx = 0

    # --- LightON panels ---
    for row_idx in range(nrows_light):
        for col_idx in range(ncols):
            panel_idx += 1
            data_idx = row_idx * ncols + col_idx
            if data_idx >= n_light:
                break
            animal, adj = lightON_available[data_idx]
            ax = fig.add_subplot(nrows_total, ncols,
                                  row_idx * ncols + col_idx + 1)
            stats = compute_edge_stats(adj)
            dn = DAYNIGHT.get(animal, "?")
            im = ax.imshow(adj, cmap="RdBu_r", aspect="auto",
                           vmin=-vmax, vmax=vmax, interpolation="none")
            ax.set_title(
                f"{animal} ({dn}) — LightON\n"
                f"N={stats['n_units']}  {stats['n_significant']}/{stats['n_units']*(stats['n_units']-1)} edges ({stats['survival']:.1%})",
                fontsize=9,
            )
            ax.set_xticks([])
            ax.set_yticks([])

    # --- Separator row ---
    sep_row = nrows_light * ncols
    ax_sep = fig.add_subplot(nrows_total, ncols,
                              nrows_light * ncols + 1)
    ax_sep.text(0.5, 0.5, "▼ Ongoing (pre-stimulus baseline) ▼",
                ha="center", va="center", fontsize=13, fontweight="bold",
                color="#333333", transform=ax_sep.transAxes)
    ax_sep.axis("off")

    # --- Ongoing panels ---
    panel_offset = nrows_light * ncols + ncols  # skip separator row cells
    for row_idx in range(nrows_ongoing):
        for col_idx in range(ncols):
            data_idx = row_idx * ncols + col_idx
            if data_idx >= n_ongoing:
                break
            animal, adj = ongoing_available[data_idx]
            ax = fig.add_subplot(nrows_total, ncols,
                                  panel_offset + row_idx * ncols + col_idx + 1)
            stats = compute_edge_stats(adj)
            dn = DAYNIGHT.get(animal, "?")
            im = ax.imshow(adj, cmap="RdBu_r", aspect="auto",
                           vmin=-vmax, vmax=vmax, interpolation="none")
            ax.set_title(
                f"{animal} ({dn}) — Ongoing\n"
                f"N={stats['n_units']}  {stats['n_significant']}/{stats['n_units']*(stats['n_units']-1)} edges ({stats['survival']:.1%})",
                fontsize=9,
            )
            ax.set_xticks([])
            ax.set_yticks([])

    # Colorbar
    cbar_ax = fig.add_axes([0.92, 0.08, 0.01, 0.84])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Connection strength (CCG peak)", fontsize=11)

    # Legend
    fig.text(0.5, 0.01,
             "Red = excitatory (positive CCG peak)  |  Blue = inhibitory (negative CCG peak)  |  "
             "5 animals missing ongoing condition (night7,day3,day4,day5,day6)",
             ha="center", fontsize=9, color="#888888")

    out_path = FIGURES_DIR / "figS1_per_animal_adjacencies.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    print(f"  LightON panels: {n_light}/{len(ANIMALS)}")
    print(f"  Ongoing panels: {n_ongoing}/{len(ANIMALS)}")

    return out_path


def plot_aggregate_summary():
    """Single-panel summary: edge survival rate across animals."""
    animals_list = []
    survivals_l = []
    survivals_o = []
    n_neurons = []

    for animal in ANIMALS:
        adj_l = load_adjacency(animal, "lightON")
        adj_o = load_adjacency(animal, "ongoing_bis")
        if adj_o is None:
            adj_o = load_adjacency(animal, "ongoing")

        animals_list.append(animal)
        n_neurons.append(NEURON_COUNTS.get(animal, 0))
        survivals_l.append(compute_edge_stats(adj_l)["survival"] if adj_l is not None else 0)
        survivals_o.append(compute_edge_stats(adj_o)["survival"] if adj_o is not None else 0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Edge Survival Rate Across Animals", fontsize=14, fontweight="bold")

    x = range(len(animals_list))
    colors = ["#2166ac" if DAYNIGHT.get(a) == "D" else "#b2182b" for a in animals_list]

    # Bar chart: survival rate
    width = 0.35
    bars_l = ax1.bar([i - width/2 for i in x], survivals_l, width,
                      label="LightON", color="#2166ac", alpha=0.8)
    bars_o = ax1.bar([i + width/2 for i in x], survivals_o, width,
                      label="Ongoing", color="#b2182b", alpha=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(animals_list, rotation=45, ha="right", fontsize=8)
    ax1.set_ylabel("Edge Survival Rate")
    ax1.set_ylim(0, 1.0)
    ax1.legend()
    ax1.set_title("Edge Survival Rate (BH-FDR q=0.05)")

    # Bar chart: neuron count
    bars_n = ax2.bar(x, n_neurons, color=colors, alpha=0.8)
    ax2.set_xticks(x)
    ax2.set_xticklabels(animals_list, rotation=45, ha="right", fontsize=8)
    ax2.set_ylabel("Number of Neurons")
    ax2.set_title("Neuron Yield per Animal")
    # Add day/night legend
    d_patch = mpatches.Patch(color="#2166ac", label="Day")
    n_patch = mpatches.Patch(color="#b2182b", label="Night")
    ax2.legend(handles=[d_patch, n_patch])

    out_path = FIGURES_DIR / "figS1_edge_survival_summary.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")

    return out_path


if __name__ == "__main__":
    plot_all_heatmaps()
    plot_aggregate_summary()
    print("\nDone.")
