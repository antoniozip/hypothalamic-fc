#!/usr/bin/env python3
"""Figure for the jitter sensitivity sweep.

Panel A is the result: edge survival against displacement, with the full within-window
shuffle as the floor and the expectation for a genuinely synaptic edge set drawn for
contrast. Panel B shows the same thing in the weights rather than the edge set. Panel C
shows what the edge set does depend on.

Input:  results/jitter_sensitivity.csv
Output: figures/main/figJ1_jitter_sensitivity.png
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

LEVELS = ["0", "5", "10", "25", "50", "100", "shuffle"]
XPOS = [0, 1, 2, 3, 4, 5, 6.8]          # shuffle set apart: it is not on the ms axis
LABELS = ["0", "5", "10", "25", "50", "100", "full\nshuffle"]


def mechanism_panel(ax):
    """Spike-count product for pairs that do and do not get an edge."""
    import shutil
    from gap_spike_trains import gap_directory_windows, stimulus_windows

    work = Path("/tmp/claude-1000/mech_fig")
    shutil.rmtree(work, ignore_errors=True)
    data = []
    for animal, old in [("night5", "180111"), ("day3", "180302"), ("night6", "180221")]:
        mat = PROJECT_ROOT / f"DATA{old}_ongoing.mat"
        adj = PROJECT_ROOT / "results" / "glmcc" / "sensitivity" / f"adj_{animal}_baseline.csv"
        if not (mat.exists() and adj.exists()):
            continue
        s, e = stimulus_windows(mat)
        st = work / animal
        gap_directory_windows(PROJECT_ROOT / f"DATA{old}_ongoing", st, s + 45.0,
                              s + 45.0 + (e - s))
        counts = np.array([np.atleast_1d(np.loadtxt(p)).size
                           for p in sorted(st.glob("cell*.txt"),
                                           key=lambda p: int(p.stem[4:]))])
        A = np.loadtxt(adj, delimiter=",")
        n = len(counts)
        off = ~np.eye(n, dtype=bool)
        prod = np.outer(counts, counts)[off].astype(float)
        present = A[off] != 0
        data.append((animal, prod[present], prod[~present]))
    shutil.rmtree(work, ignore_errors=True)

    pos, ticks = [], []
    for i, (animal, yes, no) in enumerate(data):
        for j, (vals, col) in enumerate(((no, "#9aa5b1"), (yes, "#2b6cb0"))):
            v = vals[vals > 0]
            if not v.size:
                continue
            p = i * 2.6 + j
            bp = ax.boxplot([np.log10(v)], positions=[p], widths=0.7,
                            patch_artist=True, showfliers=False)
            bp["boxes"][0].set_facecolor(col)
            bp["boxes"][0].set_alpha(0.85)
            bp["medians"][0].set_color("black")
            pos.append(p)
        ticks.append((i * 2.6 + 0.5, animal))
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([t[1] for t in ticks])
    ax.set_ylabel("log$_{10}$ ( $n_i \\times n_j$ )")
    ax.set_title("C  Spike count still gates detection", loc="left", fontsize=11,
                 fontweight="bold")
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc="#9aa5b1", alpha=.85),
                       plt.Rectangle((0, 0), 1, 1, fc="#2b6cb0", alpha=.85)],
              labels=["no edge", "edge"], frameon=False, fontsize=9, loc="lower right")


def main() -> int:
    df = pd.read_csv(PROJECT_ROOT / "results" / "jitter_sensitivity.csv")
    df["level"] = df["level"].astype(str)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # ---- A: edge survival --------------------------------------------------
    ax = axes[0]
    for animal, g in df.groupby("animal"):
        g = g.set_index("level").reindex(LEVELS)
        ax.plot(XPOS, g.survival_pct.values, "-", color="#b9c2cc", lw=1, zorder=1)
    m = [df[df.level == lv].survival_pct.mean() for lv in LEVELS]
    sd = [df[df.level == lv].survival_pct.std() for lv in LEVELS]
    ax.errorbar(XPOS, m, yerr=sd, fmt="o-", color="#2b6cb0", lw=2.2, ms=6,
                capsize=3, zorder=3, label="mean ± SD (n=13)")
    ax.axhspan(0, 5, color="#c53030", alpha=0.07, zorder=0)
    ax.text(5.6, 7.5, "floor: full shuffle", fontsize=8, color="#c53030", ha="right")
    ax.axvline(6.0, color="0.8", lw=1, ls=":")
    ax.set_xticks(XPOS)
    ax.set_xticklabels(LABELS)
    ax.set_xlabel("spike displacement (ms, uniform ±)")
    ax.set_ylabel("edges surviving (% of baseline)")
    ax.set_ylim(-4, 106)
    ax.set_title("A  Edges collapse between 5 and 25 ms", loc="left",
                 fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=9, loc="upper right")

    ax.grid(alpha=.3)

    # ---- B: how many edges the jittered data produces at all ---------------
    # Not the weight correlation: past 10 ms so few edges survive that r is computed
    # on a handful of pairs per animal and is pure noise. Absolute counts are stable.
    ax = axes[1]
    base_n = [df[df.level == lv].n_base.sum() for lv in LEVELS]
    jit_n = [df[df.level == lv].n_jit.sum() for lv in LEVELS]
    ov_n = [df[df.level == lv].overlap.sum() for lv in LEVELS]
    ax.plot(XPOS, base_n, "o--", color="#718096", lw=1.6, ms=5, label="baseline edges")
    ax.plot(XPOS, jit_n, "o-", color="#c05621", lw=2.2, ms=6,
            label="edges found after jitter")
    ax.plot(XPOS, ov_n, "o-", color="#2b6cb0", lw=2.2, ms=6,
            label="of those, also in baseline")
    ax.axvline(6.0, color="0.8", lw=1, ls=":")
    ax.set_xticks(XPOS)
    ax.set_xticklabels(LABELS)
    ax.set_xlabel("spike displacement (ms, uniform ±)")
    ax.set_ylabel("edges, summed over 13 animals")
    ax.set_title("B  Jittered data yields far fewer edges", loc="left", fontsize=11,
                 fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=.3)

    # ---- C: mechanism ------------------------------------------------------
    try:
        mechanism_panel(axes[2])
    except Exception as exc:                       # data may be unavailable
        axes[2].text(.5, .5, f"panel unavailable:\n{exc}", ha="center", va="center",
                     transform=axes[2].transAxes, fontsize=8)

    fig.suptitle("Jitter sensitivity (correct millisecond time base): edges are timing-dependent",
                 fontsize=12.5, fontweight="bold", y=1.0)
    fig.tight_layout()
    out = PROJECT_ROOT / "figures" / "main" / "figJ1_jitter_sensitivity.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
