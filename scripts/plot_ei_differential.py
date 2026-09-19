#!/usr/bin/env python3
"""Figure for the E/I differential analysis.

A  who is classifiable at all
B  the expected firing-rate confound, which is absent
C  why the sign validation cannot work: the sign distribution is degenerate
D  edge-type densities, paired per animal -- pooled vs paired disagree

Output: figures/main/figE1_ei_differential.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.3))

    # ---- A: classifiability ------------------------------------------------
    t = pd.read_csv(ROOT / "data" / "processed" / "neurons_ei_strict.csv")
    frac = t.groupby("animal").pct_inh_animal.first().sort_values()
    ax = axes[0]
    cols = ["#c53030" if v < 10 else "#2f855a" for v in frac.values]
    ax.barh(range(len(frac)), frac.values, color=cols)
    ax.set_yticks(range(len(frac)))
    ax.set_yticklabels(frac.index, fontsize=8)
    ax.axvline(10, color="0.4", ls="--", lw=1.2)
    ax.set_xlabel("% of units labelled inhibitory")
    ax.set_title("A  Only 7 of 13 were classified", loc="left", fontsize=11,
                 fontweight="bold")
    ax.text(10.8, 0.2, "threshold", fontsize=7.5, color="0.35")
    ax.text(.97, .45, "classification\nnot run", transform=ax.transAxes, ha="right",
            fontsize=8, color="#c53030")
    ax.grid(alpha=.3, axis="x")

    # ---- B: the confound that isn't ---------------------------------------
    r = pd.read_csv(ROOT / "results" / "ei" / "rates.csv")
    p = r.pivot_table(index=["animal", "condition"], columns="ei", values="median_spikes")
    ax = axes[1]
    for a, b in zip(p["E"], p["I"]):
        ax.plot([0, 1], [a, b], "-", color="#b9c2cc", lw=1, zorder=1)
    ax.scatter(np.zeros(len(p)), p["E"], c="#2b6cb0", s=45, zorder=3, edgecolor="white",
               linewidth=.6)
    ax.scatter(np.ones(len(p)), p["I"], c="#d69e2e", s=45, zorder=3, edgecolor="white",
               linewidth=.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["excitatory", "inhibitory"])
    ax.set_xlim(-.35, 1.35)
    ax.set_ylabel("median spikes per unit")
    ax.set_title("B  No firing-rate difference", loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("ratio I/E = 1.01x, paired p = 0.33\n"
                  "narrow-spiking units normally fire faster", fontsize=8, color="0.3")
    ax.grid(alpha=.3, axis="y")

    # ---- C: degenerate sign distribution -----------------------------------
    s = pd.read_csv(ROOT / "results" / "ei" / "outgoing_sign.csv")
    ax = axes[2]
    bins = np.linspace(0, 1, 21)
    for k, c, lab in (("E", "#2b6cb0", "excitatory"), ("I", "#d69e2e", "inhibitory")):
        ax.hist(s[s.ei == k].frac_neg, bins=bins, alpha=.6, color=c, label=lab)
    ax.set_xlabel("fraction of a unit's outgoing couplings that are negative")
    ax.set_ylabel("units")
    ax.set_title("C  The sign carries no information", loc="left", fontsize=11,
                 fontweight="bold")
    ax.text(.03, .95, "94.9% of units emit\nonly negative couplings\n\nE 0.988 vs I 0.984\np = 0.93",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox=dict(fc="white", ec="0.8", boxstyle="round,pad=0.35"))
    ax.legend(frameon=False, fontsize=9, loc="center left")

    # ---- D: edge types, pooled vs paired -----------------------------------
    e = pd.read_csv(ROOT / "results" / "ei" / "edge_types.csv")
    e["type"] = e.pre + "→" + e.post
    order = ["E→E", "E→I", "I→E", "I→I"]
    pooled = (e.groupby("type").apply(lambda g: 100 * g.edges.sum() / g.possible.sum())
              .reindex(order))
    pv = e.pivot_table(index=["animal", "condition"], columns="type",
                       values="density_pct").dropna()[order]
    ax = axes[3]
    x = np.arange(4)
    for i, row in pv.iterrows():
        ax.plot(x, row.values, "-", color="#b9c2cc", lw=.9, zorder=1)
    ax.plot(x, pooled.values, "o-", color="#c53030", lw=2.2, ms=7, zorder=3,
            label="pooled over animals")
    ax.plot(x, pv.median().values, "s-", color="#2f855a", lw=2.2, ms=7, zorder=3,
            label="median of per-animal")
    ax.set_xticks(x)
    ax.set_xticklabels(order)
    ax.set_ylabel("edge density (%)")
    ax.set_title("D  The gradient is a pooling artifact", loc="left", fontsize=11,
                 fontweight="bold")
    ax.set_xlabel("pooled shows E→E > I→I; paired does not\n"
                  "E→E > I→I in only 6 of 14, all Holm p = 1.0", fontsize=8, color="0.3")
    ax.legend(frameon=False, fontsize=8.5)
    ax.grid(alpha=.3, axis="y")

    fig.suptitle("Putative E/I class explains nothing in this connectivity",
                 fontsize=12.5, fontweight="bold", y=1.0)
    fig.tight_layout()
    out = ROOT / "figures" / "main" / "figE1_ei_differential.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
