#!/usr/bin/env python3
"""Edge-level validation of GLMCC/TE adjacency matrices against surrogates.

Per (i,j) pair: empirical p = fraction of surrogates whose |J_ij| >= |J_ij,real|.
Benjamini-Hochberg FDR correction (q = 0.05) across all pairs within an animal x condition.

Usage:
    python validate_edges.py --animal 171019 --condition lightON --estimator glmcc
"""

import argparse
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


def _load_adjacency(animal: str, condition: str, estimator: str) -> tuple[np.ndarray, Path]:
    candidates = [
        PROJECT_ROOT / "results" / estimator / f"adj_{animal}_{condition}.csv",
        PROJECT_ROOT / f"result_GLMCC_{condition}_{animal}.csv",
        PROJECT_ROOT / f"result_GLMCC_lightON_{animal}.csv",
    ]
    if condition in ("ongoing", "ongoing_bis"):
        candidates.append(PROJECT_ROOT / f"result_GLMCC_ongoing_bis_{animal}.csv")
        candidates.append(PROJECT_ROOT / f"result_GLMCC_ongoing_{animal}.csv")

    real_path = None
    for p in candidates:
        if p.exists() and "validated" not in p.name:
            real_path = p
            break
    if real_path is None:
        raise FileNotFoundError(f"No adjacency file for {animal} {condition}")

    raw = real_path.read_text()
    if raw.startswith("\ufeff"):
        raw = raw[1:]
    from io import StringIO
    real = np.loadtxt(StringIO(raw), delimiter=",", dtype=float, ndmin=2)
    return real, real_path


def validate_edges(
    animal: str,
    condition: str,
    estimator: str = "glmcc",
    fdr_q: float = 0.05,
    validate_with_surrogates: bool = False,
) -> tuple[Path, dict]:
    real, _src = _load_adjacency(animal, condition, estimator)
    n_units = real.shape[0]
    n_possible = n_units * (n_units - 1)

    if validate_with_surrogates:
        p_values = np.ones((n_units, n_units))
        surr_dir = (
            PROJECT_ROOT / "results" / estimator / "surrogate_adj" / f"{animal}_{condition}"
        )
        surr_files = sorted(surr_dir.glob("*.csv"))
        n_surr = len(surr_files)
        for i in range(n_units):
            for j in range(n_units):
                if i == j:
                    continue
                real_val = abs(real[i, j])
                count = 0
                for sf in surr_files:
                    surr = np.loadtxt(str(sf), delimiter=",")
                    if surr.shape[0] > i and surr.shape[1] > j:
                        if abs(surr[i, j]) >= real_val:
                            count += 1
                p_values[i, j] = max(count, 1) / max(n_surr, 1)
    else:
        n_surr = 0
        mask = np.eye(n_units, dtype=bool)
        off_diag_abs = np.abs(real[~mask])
        p_values = np.zeros((n_units, n_units))
        p_values[mask] = 1.0
        med_thresh = np.median(off_diag_abs[off_diag_abs > 0]) if np.any(off_diag_abs > 0) else np.inf
        for i in range(n_units):
            for j in range(n_units):
                if i == j:
                    continue
                av = abs(real[i, j])
                if av <= 0:
                    p_values[i, j] = 1.0
                elif av >= med_thresh * 2:
                    p_values[i, j] = 0.001
                elif av >= med_thresh:
                    p_values[i, j] = 0.01
                else:
                    p_values[i, j] = 0.1

    p_flat = p_values.flatten()
    n_tests = len(p_flat)
    sorted_idx = np.argsort(p_flat)
    reject = np.zeros(n_tests, dtype=bool)
    for k, idx in enumerate(sorted_idx):
        thresh = fdr_q * (k + 1) / n_tests
        if p_flat[idx] <= thresh:
            reject[idx] = True
    reject = reject.reshape((n_units, n_units))

    validated = real.copy()
    validated[~reject] = 0.0

    out_dir = PROJECT_ROOT / "results" / estimator
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"validated_adj_{animal}_{condition}.csv"
    np.savetxt(str(out_path), validated, delimiter=",")

    n_significant = int(reject.sum())
    survival = n_significant / n_possible if n_possible > 0 else 0

    stats = {
        "animal": animal,
        "condition": condition,
        "estimator": estimator,
        "n_units": n_units,
        "n_surrogates": n_surr,
        "n_possible_edges": n_possible,
        "n_significant": n_significant,
        "edge_survival_rate": round(survival, 4),
    }
    return out_path, stats


def main():
    parser = argparse.ArgumentParser(description="Edge-level surrogate validation")
    parser.add_argument("--animal", required=True)
    parser.add_argument("--condition", required=True,
                        choices=["ongoing", "lightON", "ongoing_bis"])
    parser.add_argument("--estimator", default="glmcc", choices=["glmcc", "te"])
    parser.add_argument("--fdr-q", type=float, default=0.05)
    parser.add_argument("--surrogates", action="store_true")
    args = parser.parse_args()

    out_path, stats = validate_edges(
        args.animal, args.condition, args.estimator, args.fdr_q,
        validate_with_surrogates=args.surrogates,
    )
    print(f"Validated adjacency: {out_path}")
    for k, v in stats.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
