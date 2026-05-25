#!/usr/bin/env python3
"""Tier 1 #5: Firing Rate × Connectivity Coupling.

Does a neuron's firing rate predict its node strength or hub score?
If FR and strength are uncorrelated, detected edges reflect genuine
temporal coupling rather than rate artifacts.

Output: figures/main/figT1_fr_connectivity.png
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
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

TS_KEY_MAP = {"lightON": "evoked_ts", "ongoing": "ongoing_ts"}


def load_spike_times(animal: str, condition: str) -> list | None:
    """Load spike times (seconds) from .mat file."""
    mat_path = PROJECT_ROOT / f"DATA{animal}_{condition}.mat"
    if not mat_path.exists():
        if condition == "ongoing":
            alts = sorted(PROJECT_ROOT.glob(f"DATA{animal}*ongoing*.mat"))
            if alts:
                mat_path = alts[0]
    if not mat_path.exists():
        return None

    m = loadmat(str(mat_path), simplify_cells=False)
    ts_key = TS_KEY_MAP.get(condition, "evoked_ts")
    if ts_key not in m:
        ts_keys = [k for k in m if k.endswith("_ts")]
        if not ts_keys:
            return None
        ts_key = ts_keys[0]

    ts = m[ts_key]
    n = min(ts.shape[1], 64)
    times = []
    for i in range(n):
        try:
            t = np.atleast_1d(ts[0, i][0]).astype(float)
            if len(t) > 0 and t.max() > 100:
                t = t / 1000.0
            times.append(t)
        except (IndexError, ValueError):
            times.append(np.array([]))
    return times


def load_adjacency(animal: str) -> np.ndarray | None:
    """Load validated adjacency matrix."""
    adj_path = PROJECT_ROOT / "results" / "glmcc" / f"validated_adj_{animal}_lightON.csv"
    if not adj_path.exists():
        return None
    raw = adj_path.read_text()
    if raw.startswith("\ufeff"):
        raw = raw[1:]
    from io import StringIO
    return np.loadtxt(StringIO(raw), delimiter=",", dtype=float)


def load_metrics(animal: str) -> pd.DataFrame | None:
    """Load graph metrics CSV."""
    mpath = PROJECT_ROOT / "results" / "glmcc" / f"metrics_{animal}_lightON.csv"
    if not mpath.exists():
        return None
    return pd.read_csv(mpath)


def load_ei(animal: str) -> dict:
    """Load E-I classification. Returns {neuron_id: 'E' or 'I'}."""
    path = PROJECT_ROOT / f"{animal}_put_inh.txt"
    if not path.exists():
        return {}
    text = path.read_text().strip()
    if not text:
        return {}
    inh_ids = set()
    for part in text.replace(",", " ").split():
        if part.strip().isdigit():
            inh_ids.add(int(part.strip()))
    return {i: "I" if i in inh_ids else "E" for i in range(1, 65)}


def main():
    results = []

    for animal in ANIMALS:
        # Load data
        times = load_spike_times(animal, "lightON")
        adj = load_adjacency(animal)
        metrics_df = load_metrics(animal)
        ei_map = load_ei(animal)

        if times is None or adj is None:
            print(f"  SKIP {animal} — missing data")
            continue

        n_units = min(len(times), adj.shape[0])
        duration = times[0][-1] - times[0][0] if len(times[0]) > 1 else 1.0

        for i in range(n_units):
            fr = len(times[i]) / max(duration, 1.0)
            ns = float(np.sum(np.abs(adj[i, :])))  # node strength

            # Region and EI from metrics CSV if available
            region = "unknown"
            ei = ei_map.get(i + 1, "E")
            if metrics_df is not None and i < len(metrics_df):
                row = metrics_df.iloc[i]
                region = row.get("region", "unknown")
                cc = row.get("clustering_coefficient", np.nan)
                le = row.get("local_efficiency", np.nan)
                hs = row.get("hub_score", np.nan)
            else:
                cc = np.nan
                le = np.nan
                hs = np.nan

            results.append({
                "animal": animal,
                "day_night": DAYNIGHT.get(animal, "?"),
                "neuron_id": i + 1,
                "region": region,
                "ei": ei,
                "firing_rate": fr,
                "node_strength": ns,
                "clustering_coefficient": cc,
                "local_efficiency": le,
                "hub_score": hs,
            })

    df = pd.DataFrame(results)
    print(f"  {len(df)} neurons from {df['animal'].nunique()} animals")

    # --- Per-animal FR × Strength Spearman ρ ---
    rho_values = []
    for animal in df["animal"].unique():
        adf = df[df["animal"] == animal]
        if len(adf) < 5:
            continue
        rho, p = stats.spearmanr(adf["firing_rate"], adf["node_strength"])
        rho_values.append({
            "animal": animal,
            "day_night": DAYNIGHT.get(animal, "?"),
            "n": len(adf),
            "rho": rho,
            "p_value": p,
        })

    rho_df = pd.DataFrame(rho_values)
    print(f"\n  Mean FR×Strength ρ: {rho_df['rho'].mean():.3f} ± {rho_df['rho'].std():.3f}")
    print(f"  Significant animals (p<0.05): {(rho_df['p_value'] < 0.05).sum()}/{len(rho_df)}")

    # --- Figure ---
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle("T5: Firing Rate × Connectivity Coupling",
                 fontsize=15, fontweight="bold")

    # Panel A: FR vs Strength scatter (all neurons, log-log)
    ax = axes[0, 0]
    valid = df[(df["firing_rate"] > 0) & (df["node_strength"] > 0)]
    ax.scatter(valid["firing_rate"], valid["node_strength"],
               alpha=0.2, s=5, c="#444444")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Firing Rate (Hz, log)")
    ax.set_ylabel("Node Strength (log)")
    ax.set_title(f"A: FR vs Strength (n={len(valid)} neurons)")

    # Add regression line
    if len(valid) > 10:
        log_fr = np.log10(valid["firing_rate"].values)
        log_ns = np.log10(valid["node_strength"].values)
        slope, intercept, r_val, p_val, _ = stats.linregress(log_fr, log_ns)
        x_line = np.logspace(np.log10(valid["firing_rate"].min()),
                             np.log10(valid["firing_rate"].max()), 100)
        y_line = 10 ** (intercept + slope * np.log10(x_line))
        ax.plot(x_line, y_line, "r-", linewidth=1.5, alpha=0.7,
                label=f"R²={r_val**2:.3f}, p={p_val:.2e}")
        ax.legend(fontsize=8)

    # Panel B: Per-animal ρ bar chart
    ax = axes[0, 1]
    rho_df_sorted = rho_df.sort_values("rho")
    colors = ["#2166ac" if r["day_night"] == "D" else "#b2182b"
              for _, r in rho_df_sorted.iterrows()]
    ax.barh(range(len(rho_df_sorted)), rho_df_sorted["rho"].values,
            color=colors, alpha=0.85)
    ax.set_yticks(range(len(rho_df_sorted)))
    ax.set_yticklabels(rho_df_sorted["animal"].values, fontsize=8)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.7)
    ax.set_xlabel("Spearman ρ (FR vs Node Strength)")
    ax.set_title(f"B: Per-animal FR×Strength ρ\nMean ρ = {rho_df['rho'].mean():.3f}")

    # Panel C: FR distribution by EI
    ax = axes[0, 2]
    for ei_val, color, label in [("E", "#2166ac", "Excitatory"),
                                   ("I", "#b2182b", "Inhibitory")]:
        subset = df[df["ei"] == ei_val]
        ax.hist(subset["firing_rate"], bins=30, alpha=0.5, color=color,
                label=label, density=True)
    ax.set_xlabel("Firing Rate (Hz)")
    ax.set_ylabel("Density")
    ax.set_title("C: FR Distribution by E-I Type")
    ax.legend(fontsize=8)

    # Panel D: FR vs Clustering Coefficient
    ax = axes[1, 0]
    valid_cc = df.dropna(subset=["clustering_coefficient", "firing_rate"])
    valid_cc = valid_cc[(valid_cc["firing_rate"] > 0)]
    ax.scatter(valid_cc["firing_rate"], valid_cc["clustering_coefficient"],
               alpha=0.2, s=5, c="#444444")
    ax.set_xscale("log")
    ax.set_xlabel("Firing Rate (Hz, log)")
    ax.set_ylabel("Clustering Coefficient")
    ax.set_title("D: FR vs Clustering Coefficient")
    if len(valid_cc) > 10:
        rho_cc, p_cc = stats.spearmanr(valid_cc["firing_rate"],
                                        valid_cc["clustering_coefficient"])
        ax.text(0.95, 0.95, f"ρ={rho_cc:.3f}, p={p_cc:.2e}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=9, bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    # Panel E: FR vs Local Efficiency
    ax = axes[1, 1]
    valid_le = df.dropna(subset=["local_efficiency", "firing_rate"])
    valid_le = valid_le[valid_le["firing_rate"] > 0]
    ax.scatter(valid_le["firing_rate"], valid_le["local_efficiency"],
               alpha=0.2, s=5, c="#444444")
    ax.set_xscale("log")
    ax.set_xlabel("Firing Rate (Hz, log)")
    ax.set_ylabel("Local Efficiency")
    ax.set_title("E: FR vs Local Efficiency")
    if len(valid_le) > 10:
        rho_le, p_le = stats.spearmanr(valid_le["firing_rate"],
                                        valid_le["local_efficiency"])
        ax.text(0.95, 0.95, f"ρ={rho_le:.3f}, p={p_le:.2e}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=9, bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    # Panel F: FR × Strength by region
    ax = axes[1, 2]
    region_rhos = []
    for region in df["region"].unique():
        rdf = df[(df["region"] == region) & (df["firing_rate"] > 0)]
        if len(rdf) < 5:
            continue
        rho, p = stats.spearmanr(rdf["firing_rate"], rdf["node_strength"])
        region_rhos.append({"region": region, "rho": rho, "n": len(rdf), "p": p})

    region_rhos = sorted(region_rhos, key=lambda x: x["rho"])
    ax.barh([r["region"] for r in region_rhos],
            [r["rho"] for r in region_rhos],
            color="#6baed6", alpha=0.85)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.7)
    ax.set_xlabel("Spearman ρ (FR vs Node Strength)")
    ax.set_title("F: FR×Strength ρ by Region")

    plt.tight_layout()
    out_path = FIGURES_DIR / "figT1_fr_connectivity.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\nSaved: {out_path}")

    # --- Summary statistics ---
    print(f"\n=== FR × Connectivity Summary ===")
    print(f"  Mean FR: {df['firing_rate'].mean():.3f} ± {df['firing_rate'].std():.3f} Hz")
    print(f"  Mean Strength: {df['node_strength'].mean():.1f} ± {df['node_strength'].std():.1f}")

    # E vs I comparison
    for metric in ["firing_rate", "node_strength"]:
        e_vals = df[df["ei"] == "E"][metric].dropna()
        i_vals = df[df["ei"] == "I"][metric].dropna()
        if len(e_vals) > 0 and len(i_vals) > 0:
            t_stat, t_p = stats.ttest_ind(e_vals, i_vals)
            print(f"  {metric}: E={e_vals.mean():.3f}, I={i_vals.mean():.3f}, "
                  f"t={t_stat:.2f}, p={t_p:.3f}")


if __name__ == "__main__":
    main()
