#!/usr/bin/env python3
"""E-I classification analysis from put_inh.txt files.

Reads putative inhibitory neuron indices, merges with neuron metadata,
and generates figures showing E-I distribution across regions and conditions.

Output:
  data/processed/neurons_ei.csv       — neurons.csv + EI column
  figures/main/figS6_ei_distribution.png
  figures/main/figS6_ei_metrics.png
"""

import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FIGURES_DIR = PROJECT_ROOT / "figures" / "main"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

ANIMALS = [
    "171019", "171207", "171208", "171213", "180110", "180111",
    "180131", "180221", "180228", "180302", "180419", "180420", "180423",
]
DAYNIGHT = {
    "171019": "D", "171207": "N", "171208": "N", "171213": "N",
    "180110": "N", "180111": "N", "180131": "D", "180221": "N",
    "180228": "N", "180302": "D", "180419": "D", "180420": "D",
    "180423": "D",
}


def parse_put_inh(animal: str) -> set[int]:
    """Parse put_inh file, return set of 1-indexed inhibitory neuron IDs."""
    path = PROJECT_ROOT / f"{animal}_put_inh.txt"
    if not path.exists():
        return set()
    text = path.read_text().strip()
    if not text:
        return set()
    ids = set()
    for part in re.split(r'[,\s]+', text):
        part = part.strip()
        if part.isdigit():
            ids.add(int(part))
    return ids


def main():
    # ─── Load neurons.csv ────────────────────────────────────────
    neurons_path = PROJECT_ROOT / "data" / "processed" / "neurons.csv"
    df = pd.read_csv(neurons_path)

    # ─── Add EI column ───────────────────────────────────────────
    ei_map = {}
    for animal in ANIMALS:
        inh_ids = parse_put_inh(animal)
        animal_neurons = df[df["animal"] == int(animal)]
        for _, row in animal_neurons.iterrows():
            nid = int(row["neuron_id"])
            key = (animal, nid)
            ei_map[key] = "I" if nid in inh_ids else "E"

    df["ei"] = df.apply(
        lambda row: ei_map.get((str(row["animal"]), int(row["neuron_id"])), "E"),
        axis=1,
    )

    # Save updated table
    out_csv = PROJECT_ROOT / "data" / "processed" / "neurons_ei.csv"
    df.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")
    print(f"  E: {(df['ei'] == 'E').sum()} neurons")
    print(f"  I: {(df['ei'] == 'I').sum()} neurons")

    # ─── Figure S6a: E-I distribution by region ──────────────────
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle("E-I Classification Analysis",
                 fontsize=15, fontweight="bold")

    # Panel A: Stacked bar — E/I count per region
    ax = axes[0, 0]
    ei_by_region = df.groupby(["region6", "ei"]).size().unstack(fill_value=0)
    # Order regions by total count
    region_order = ei_by_region.sum(axis=1).sort_values(ascending=False).index
    ei_by_region = ei_by_region.loc[region_order]

    x = np.arange(len(ei_by_region))
    w = 0.6
    ax.bar(x, ei_by_region.get("E", 0), w, label="Excitatory",
           color="#2166ac", alpha=0.85)
    ax.bar(x, ei_by_region.get("I", 0), w,
           bottom=ei_by_region.get("E", 0), label="Inhibitory",
           color="#b2182b", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(ei_by_region.index, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Neuron Count")
    ax.set_title("A: E-I distribution by hypothalamic region")
    ax.legend(fontsize=9)

    # Add percentage labels
    for i, region in enumerate(ei_by_region.index):
        total = ei_by_region.loc[region].sum()
        inh = ei_by_region.loc[region].get("I", 0)
        pct = inh / total * 100 if total > 0 else 0
        ax.text(i, total + 0.5, f"{pct:.0f}% I", ha="center",
                fontsize=7, color="#b2182b")

    # Panel B: E-I proportion by animal
    ax = axes[0, 1]
    ei_by_animal = df.groupby(["animal", "ei"]).size().unstack(fill_value=0)
    animal_order_sorted = sorted(ei_by_animal.index,
                                  key=lambda a: ei_by_animal.loc[a].get("I", 0)
                                  / ei_by_animal.loc[a].sum(),
                                  reverse=True)
    ei_by_animal = ei_by_animal.loc[animal_order_sorted]

    x = np.arange(len(ei_by_animal))
    inh_pct = [ei_by_animal.loc[a].get("I", 0) / ei_by_animal.loc[a].sum() * 100
               for a in ei_by_animal.index]
    colors = ["#2166ac" if DAYNIGHT.get(str(a)) == "D" else "#b2182b"
              for a in ei_by_animal.index]
    ax.bar(x, inh_pct, color=colors, alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(ei_by_animal.index, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Inhibitory Neurons (%)")
    ax.set_title("B: Inhibitory neuron proportion by animal (Day/Night)")
    ax.axhline(y=np.mean(inh_pct), color="#888888", linestyle="--",
               label=f"Mean = {np.mean(inh_pct):.1f}%")
    ax.legend(fontsize=8)

    # Panel C: E-I metrics — load first-level stats and split by EI
    ax = axes[1, 0]

    # Load metrics data for animals with both ongoing and lightON
    metrics_dir = PROJECT_ROOT / "results" / "glmcc"
    ei_metrics = []

    for animal in ANIMALS:
        inh_ids = parse_put_inh(animal)
        if not inh_ids:
            continue

        # Load metrics for lightON
        mpath = metrics_dir / f"metrics_{animal}_lightON.csv"
        if not mpath.exists():
            continue
        mdf = pd.read_csv(mpath)
        mdf["ei"] = "E"
        mdf.loc[mdf["neuron_id"].isin(inh_ids), "ei"] = "I"

        for _, row in mdf.iterrows():
            ei_metrics.append({
                "animal": animal,
                "day_night": DAYNIGHT.get(animal, "?"),
                "ei": row["ei"],
                "node_strength": row.get("node_strength", np.nan),
                "clustering_coefficient": row.get("clustering_coefficient", np.nan),
                "local_efficiency": row.get("local_efficiency", np.nan),
                "hub_score": row.get("hub_score", np.nan),
            })

    ei_df = pd.DataFrame(ei_metrics)
    if len(ei_df) > 0:
        # Compare E vs I for node_strength
        ei_summary = ei_df.groupby(["ei"]).agg(
            mean_ns=("node_strength", "mean"),
            sem_ns=("node_strength", "sem"),
            mean_cc=("clustering_coefficient", "mean"),
            sem_cc=("clustering_coefficient", "sem"),
            mean_le=("local_efficiency", "mean"),
            sem_le=("local_efficiency", "sem"),
        )

        metrics_names = ["Node Strength", "Clustering Coeff", "Local Efficiency"]
        e_means = [ei_summary.loc["E", "mean_ns"],
                    ei_summary.loc["E", "mean_cc"],
                    ei_summary.loc["E", "mean_le"]]
        i_means = [ei_summary.loc["I", "mean_ns"],
                    ei_summary.loc["I", "mean_cc"],
                    ei_summary.loc["I", "mean_le"]]
        e_sems = [ei_summary.loc["E", "sem_ns"],
                   ei_summary.loc["E", "sem_cc"],
                   ei_summary.loc["E", "sem_le"]]
        i_sems = [ei_summary.loc["I", "sem_ns"],
                   ei_summary.loc["I", "sem_cc"],
                   ei_summary.loc["I", "sem_le"]]

        x = np.arange(len(metrics_names))
        w = 0.35
        ax.bar(x - w/2, e_means, w, yerr=e_sems, capsize=3,
               label="Excitatory (E)", color="#2166ac", alpha=0.85)
        ax.bar(x + w/2, i_means, w, yerr=i_sems, capsize=3,
               label="Inhibitory (I)", color="#b2182b", alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(metrics_names, fontsize=9)
        ax.set_ylabel("Mean Value (± SEM)")
        ax.set_title("C: E vs I graph metrics (lightON)")
        ax.legend(fontsize=8)

    # Panel D: Per-region E-I breakdown
    ax = axes[1, 1]
    ei_region_detail = df.groupby(["region6", "ei"]).size().unstack(fill_value=0)
    ei_region_detail["pct_I"] = (
        ei_region_detail.get("I", 0)
        / (ei_region_detail.get("E", 0) + ei_region_detail.get("I", 0))
        * 100
    )
    ei_region_detail = ei_region_detail.loc[
        ei_region_detail.sum(axis=1).sort_values(ascending=False).index
    ]

    x = np.arange(len(ei_region_detail))
    ax.bar(x, ei_region_detail["pct_I"], color="#b2182b", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(ei_region_detail.index, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Inhibitory Neurons (%)")
    ax.set_title("D: Inhibitory proportion by region")
    ax.axhline(y=ei_region_detail["pct_I"].mean(), color="#888888",
               linestyle="--", linewidth=0.8,
               label=f"Mean = {ei_region_detail['pct_I'].mean():.1f}%")
    ax.legend(fontsize=8)

    plt.tight_layout()
    out_path = FIGURES_DIR / "figS6_ei_distribution.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")

    # ─── Print EI summary ─────────────────────────────────────────
    print("\nE-I Summary:")
    for region in ei_by_region.index:
        total = ei_by_region.loc[region].sum()
        inh = ei_by_region.loc[region].get("I", 0)
        exc = ei_by_region.loc[region].get("E", 0)
        print(f"  {region:20s}: {exc:3d} E + {inh:3d} I = {total:3d} total  "
              f"({inh/total*100:5.1f}% I)")


if __name__ == "__main__":
    main()
