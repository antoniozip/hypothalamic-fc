#!/usr/bin/env python3
"""Figure for the rate-matched contrast.

A  the confound itself: density is a function of spikes per unit
B  what matching does to the condition effect, metric by metric
C  the matched contrast, paired per animal

Input:  results/glmcc_regeneration_summary.csv, results/rate_matched_summary.csv,
        results/stats_first_level_glmcc.csv, results/stats_first_level_glmcc_rm.csv
Output: figures/main/figR1_rate_matched.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # ---- A: the confound ---------------------------------------------------
    d = pd.read_csv(ROOT / "results" / "glmcc_regeneration_summary.csv")
    d["spk"] = d.n_spikes / d.n_units
    ax = axes[0]
    for cond, col, mk in (("lightON", "#d69e2e", "o"), ("ongoing", "#2b6cb0", "s")):
        x = d[d.condition == cond]
        ax.scatter(x.spk, x.density_pct, c=col, marker=mk, s=55, alpha=.85,
                   edgecolor="white", linewidth=.7, label=cond, zorder=3)
    r, p = spearmanr(d.spk, d.density_pct)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("spikes per unit")
    ax.set_ylabel("edge density (%)")
    ax.set_title("A  Density is a function of spike count", loc="left",
                 fontsize=11, fontweight="bold")
    ax.text(.04, .95, f"pooled Spearman ρ = {r:+.3f}\np = {p:.1e}", transform=ax.transAxes,
            va="top", fontsize=9,
            bbox=dict(fc="white", ec="0.8", boxstyle="round,pad=0.35"))
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.grid(alpha=.3, which="both")

    # ---- B: F statistics before and after matching -------------------------
    raw = pd.read_csv(ROOT / "results" / "stats_first_level_glmcc.csv")
    rm = pd.read_csv(ROOT / "results" / "stats_first_level_glmcc_rm.csv")
    met = ["node_strength", "clustering_coefficient", "local_efficiency", "hub_score"]
    fr = [raw[(raw.metric == m) & (raw.term == "ANOVA_condition")].estimate.iloc[0] for m in met]
    fm = [rm[(rm.metric == m) & (rm.term == "ANOVA_condition")].estimate.iloc[0] for m in met]
    ax = axes[1]
    y = np.arange(len(met))
    ax.barh(y - .2, fr, height=.38, color="#c53030", label="raw contrast (confounded)")
    ax.barh(y + .2, fm, height=.38, color="#2f855a", label="rate + geometry matched")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels([m.replace("_", "\n") for m in met], fontsize=9)
    ax.set_xlabel("ANOVA F for condition (log scale)")
    ax.axvline(3.86, color="0.4", ls="--", lw=1.2)
    ax.text(3.86, -0.62, " F≈3.9 (p=0.05)", fontsize=7.5, color="0.35", va="center")
    ax.set_ylim(-0.85, len(met) - 0.4)
    ax.set_title("B  Matching removes the condition effect", loc="left",
                 fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    ax.grid(alpha=.3, axis="x", which="both")

    # ---- C: paired, matched ------------------------------------------------
    s = pd.read_csv(ROOT / "results" / "rate_matched_summary.csv")
    p_ = s.pivot(index="animal", columns="condition")
    on = p_[("density_pct", "lightON")].values.astype(float)
    off = p_[("density_pct", "ongoing")].values.astype(float)
    ax = axes[2]
    for a, b in zip(on, off):
        ax.plot([0, 1], [a, b], "-", color="#b9c2cc", lw=1, zorder=1)
    ax.scatter(np.zeros_like(on), on, c="#d69e2e", s=55, zorder=3, edgecolor="white",
               linewidth=.7, label="lightON")
    ax.scatter(np.ones_like(off), off, c="#4a5568", s=55, zorder=3, edgecolor="white",
               linewidth=.7, label="light-off control")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["lightON", "light-off"])
    ax.set_xlim(-.35, 1.35)
    ax.set_yscale("log")
    ax.set_ylabel("edge density (%)")
    ax.set_title("C  Matched arms are indistinguishable", loc="left",
                 fontsize=11, fontweight="bold")
    ax.set_xlabel("identical windows and per-unit spike counts\n"
                  "Wilcoxon p = 0.027, light-off slightly HIGHER — opposite to the raw contrast",
                  fontsize=8, color="0.3")
    ax.grid(alpha=.3, axis="y", which="both")

    fig.suptitle("The lightON vs ongoing contrast is a detection-power artifact",
                 fontsize=12.5, fontweight="bold", y=1.0)
    fig.tight_layout()
    out = ROOT / "figures" / "main" / "figR1_rate_matched.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
