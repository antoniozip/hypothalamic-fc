"""Regression tests for validated adjacency matrix consistency.

Guards the defects found on 2026-09-18:
  1. lightON matrices held weighted CCG peak counts while ongoing held binary
     masks, so `node_strength ~ condition` compared encodings, not conditions.
  2. Producers disagreed on spike-time units (.mat in ms vs DATA*/cell*.txt in
     seconds) while sharing WINDOW_MS = 50.0, giving a +/-50 s window on one path.
  3. Empty cell files were dropped silently, shifting every later neuron index
     relative to data/processed/neurons.csv.

Run: python tests/test_validated_consistency.py
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "glmcc"
NEURONS_CSV = PROJECT_ROOT / "data" / "processed" / "neurons.csv"

FNAME_RE = re.compile(
    r"^validated_adj_(?!ccg_)(?P<animal>.+?)_(?P<condition>lightON|ongoing|ongoing_bis)\.csv$"
)


def _encoding(matrix: np.ndarray) -> str:
    """Classify a validated adjacency as 'binary' or 'weighted'."""
    values = set(np.unique(matrix).tolist())
    return "binary" if values <= {0.0, 1.0} else "weighted"


def _load_validated() -> dict[tuple[str, str], np.ndarray]:
    out: dict[tuple[str, str], np.ndarray] = {}
    for path in sorted(RESULTS_DIR.glob("validated_adj_*.csv")):
        match = FNAME_RE.match(path.name)
        if match is None:  # skips the validated_adj_ccg_* family
            continue
        out[(match["animal"], match["condition"])] = np.loadtxt(path, delimiter=",")
    return out


def _neuron_counts() -> Counter:
    counts: Counter = Counter()
    with NEURONS_CSV.open() as handle:
        for row in csv.DictReader(handle):
            counts[row["animal"]] += 1
    return counts


def test_single_encoding_across_conditions() -> list[str]:
    """Every validated matrix must use the same encoding."""
    matrices = _load_validated()
    by_encoding: dict[str, list[str]] = {}
    for (animal, condition), matrix in matrices.items():
        by_encoding.setdefault(_encoding(matrix), []).append(f"{animal}_{condition}")

    if len(by_encoding) <= 1:
        return []
    return [
        "mixed encodings across validated matrices: "
        + "; ".join(f"{kind}={len(names)}" for kind, names in sorted(by_encoding.items()))
    ]


def test_dimensions_match_neuron_table() -> list[str]:
    """Matrix rank must equal the neuron count in neurons.csv for that animal.

    Positional region mapping in src/r/compute_graph_metrics.R assumes
    matrix row i corresponds to neuron_id i.
    """
    matrices = _load_validated()
    expected = _neuron_counts()
    failures = []
    for (animal, condition), matrix in sorted(matrices.items()):
        want = expected.get(animal)
        if want is None:
            failures.append(f"{animal}_{condition}: animal absent from neurons.csv")
        elif matrix.shape[0] != want:
            failures.append(
                f"{animal}_{condition}: matrix n={matrix.shape[0]} but neurons.csv has {want}"
            )
    return failures


def test_paired_conditions_same_rank() -> list[str]:
    """lightON and ongoing for one animal must index the same neurons.

    stats_first_level.R nests (1 | animal/neuron); the pairing is meaningless
    if neuron_id refers to different units in the two conditions.
    """
    matrices = _load_validated()
    animals = {animal for animal, _ in matrices}
    failures = []
    for animal in sorted(animals):
        ranks = {
            condition: matrices[(animal, condition)].shape[0]
            for condition in ("lightON", "ongoing")
            if (animal, condition) in matrices
        }
        if len(set(ranks.values())) > 1:
            failures.append(f"{animal}: rank differs across conditions {ranks}")
    return failures


def main() -> int:
    checks = [
        ("single encoding across conditions", test_single_encoding_across_conditions),
        ("dimensions match neurons.csv", test_dimensions_match_neuron_table),
        ("paired conditions same rank", test_paired_conditions_same_rank),
    ]
    total_failures = 0
    for name, check in checks:
        failures = check()
        if failures:
            total_failures += len(failures)
            print(f"FAIL: {name}")
            for failure in failures:
                print(f"    {failure}")
        else:
            print(f"PASS: {name}")
    print(f"\n{total_failures} failure(s)")
    return 1 if total_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
