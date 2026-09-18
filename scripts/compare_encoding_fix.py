#!/usr/bin/env python3
"""Compare validated adjacencies before and after the encoding fix.

Diagnostic only -- the authoritative metrics come from src/r/compute_graph_metrics.R.
This exists to answer one question: is the lightON/ongoing difference in node
strength still an artefact of how the matrices were encoded?

Usage: python scripts/compare_encoding_fix.py --old <dir with results/glmcc>
"""

import argparse
import statistics
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ANIMALS = ["day1", "day2", "day3", "day4", "day5", "day6",
           "night1", "night2", "night3", "night4", "night5", "night6", "night7"]
CONDITIONS = ["lightON", "ongoing"]


def _encoding(matrix: np.ndarray) -> str:
    return "binary" if set(np.unique(matrix).tolist()) <= {0.0, 1.0} else "weighted"


def _summarise(root: Path) -> dict:
    """Node strength per condition, the quantity stats_first_level.R models."""
    out: dict[str, dict] = {}
    for condition in CONDITIONS:
        strengths, encodings, ranks = [], set(), []
        for animal in ANIMALS:
            path = root / "glmcc" / f"validated_adj_{animal}_{condition}.csv"
            if not path.exists():
                continue
            adj = np.loadtxt(path, delimiter=",", ndmin=2)
            encodings.add(_encoding(adj))
            ranks.append(adj.shape[0])
            # igraph strength(mode="all") on a directed weighted graph
            strengths.extend((np.abs(adj).sum(axis=0) + np.abs(adj).sum(axis=1)).tolist())
        if strengths:
            out[condition] = {
                "n_neurons": len(strengths),
                "mean": statistics.mean(strengths),
                "median": statistics.median(strengths),
                "max": max(strengths),
                "encodings": sorted(encodings),
                "n_matrices": len(ranks),
            }
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old", type=Path, required=True,
                        help="directory containing the pre-fix results/ tree")
    args = parser.parse_args()

    for label, root in (("BEFORE", args.old / "results"), ("AFTER", PROJECT_ROOT / "results")):
        summary = _summarise(root)
        print(f"\n{label}")
        for condition, stats in summary.items():
            print(f"  {condition:8s} matrices={stats['n_matrices']:2d} "
                  f"neurons={stats['n_neurons']:4d} encoding={','.join(stats['encodings']):16s} "
                  f"mean_strength={stats['mean']:12.2f} median={stats['median']:10.2f}")
        if len(summary) == 2:
            a, b = (summary[c]["mean"] for c in CONDITIONS)
            ratio = a / b if b else float("inf")
            print(f"  -> lightON/ongoing mean node strength ratio: {ratio:.2f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
