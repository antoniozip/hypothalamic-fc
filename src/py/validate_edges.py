#!/usr/bin/env python3
"""Edge-level validation of adjacency matrices against surrogates.

Two methods:
  --method ccg    Use CCG peaks with shuffle-ISI surrogates (PREFERRED for lightON).
                  Calls into validate_ccg module. Generates surrogates on-the-fly.
  --method glmcc  Validate against pre-computed GLMCC surrogates on disk.
                  Requires results/glmcc/surrogate_adj/<animal>_<condition>/*.csv

Benjamini-Hochberg FDR correction (q = 0.05) applied to OFF-DIAGONAL pairs only.

Usage:
    python validate_edges.py --animal day1 --condition lightON --method ccg
    python validate_edges.py --animal day1 --condition lightON --method glmcc --surrogates
"""

import argparse
import sys
from pathlib import Path

import numpy as np

# Ensure src/py is on the path for cross-module imports
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

# Import CCG validation functions (shared with validate_ccg.py)
from validate_ccg import (
    load_real_spike_times,
    ccg_peak_matrix,
    compute_empirical_p_values,
    apply_bh_fdr,
)


class RankMismatchError(ValueError):
    """Raised when an adjacency matrix cannot be aligned to the neuron table."""


PROJECT_ROOT = _SCRIPT_DIR.parent.parent


def _load_adjacency(animal: str, condition: str, estimator: str) -> tuple[np.ndarray, Path]:
    """Load the raw adjacency matrix, trying multiple possible paths."""
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


def _align_to_units(
    adj: np.ndarray, spike_times: list[np.ndarray], source: Path,
) -> np.ndarray:
    """Expand an adjacency matrix so row i corresponds to unit i.

    GLMCC was run on the spike-bearing units only, so a recording containing
    spike-free units yields a matrix smaller than the unit count. Reinserting a
    zero row and column at each spike-free index restores the correspondence
    with neuron_id in data/processed/neurons.csv; a unit with no spikes has no
    connections, so zero is the correct value rather than a placeholder.
    """
    n_units = len(spike_times)
    if adj.shape[0] == n_units:
        return adj

    empty_idx = [i for i, t in enumerate(spike_times) if len(t) == 0]
    if adj.shape[0] + len(empty_idx) != n_units:
        raise RankMismatchError(
            f"{source.name}: matrix is {adj.shape[0]}x{adj.shape[0]} but the "
            f"recording has {n_units} units and {len(empty_idx)} of them are "
            f"spike-free. The matrix does not describe this unit set."
        )

    keep = [i for i in range(n_units) if i not in set(empty_idx)]
    aligned = np.zeros((n_units, n_units), dtype=float)
    aligned[np.ix_(keep, keep)] = adj
    print(
        f"  Aligned {source.name}: {adj.shape[0]} -> {n_units} units "
        f"(zero rows inserted at unit index {empty_idx})",
        flush=True,
    )
    return aligned


def _validate_with_surrogates(
    real: np.ndarray, animal: str, condition: str, estimator: str,
    fdr_q: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """Validate against pre-computed surrogate adjacency files on disk."""
    n_units = real.shape[0]
    p_values = np.ones((n_units, n_units))
    surr_dir = (
        PROJECT_ROOT / "results" / estimator / "surrogate_adj"
        / f"{animal}_{condition}"
    )
    surr_files = sorted(surr_dir.glob("*.csv"))
    n_surr = len(surr_files)

    if n_surr == 0:
        raise FileNotFoundError(
            f"No surrogate files found in {surr_dir}. "
            f"Generate surrogates first or use --method ccg."
        )

    print(f"  Using {n_surr} pre-computed surrogates from {surr_dir}", flush=True)
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

    return p_values, n_units


def _validate_with_ccg(
    animal: str, condition: str, fdr_q: float = 0.05, n_surrogates: int = 100,
) -> tuple[np.ndarray, int, list[np.ndarray]]:
    """Validate using CCG peaks + shuffle-ISI surrogates (on-the-fly)."""
    print("  Loading spike trains from .mat files...", flush=True)
    real_times = load_real_spike_times(animal, condition)

    print("  Computing CCG peak matrix...", flush=True)
    real_peaks = ccg_peak_matrix(real_times)

    print(f"  Generating {n_surrogates} shuffle-ISI surrogates and computing p-values...",
          flush=True)
    p_values, n_units = compute_empirical_p_values(
        real_times, real_peaks, n_surrogates=n_surrogates
    )

    return p_values, n_units, real_times


def validate_edges(
    animal: str,
    condition: str,
    estimator: str = "glmcc",
    method: str = "ccg",
    fdr_q: float = 0.05,
    use_surrogates: bool = False,
    n_surrogates: int = 100,
) -> tuple[Path, dict]:
    """Validate edges and write validated adjacency + stats.

    Args:
        method: 'ccg' (CCG peaks + shuffle-ISI) or 'glmcc' (pre-computed surrogates)
        use_surrogates: when method='glmcc', use surrogate files on disk

    Returns:
        (output_path, stats_dict)
    """
    # --- Get p-values ---
    if method == "ccg":
        print(f"Validating {animal}/{condition} with CCG + shuffle-ISI surrogates...", flush=True)
        p_values, n_units, real_times = _validate_with_ccg(
            animal, condition, fdr_q, n_surrogates
        )
        real, real_path = _load_adjacency(animal, condition, estimator)
        real = _align_to_units(real, real_times, real_path)
        if real.shape[0] != n_units:
            raise RankMismatchError(
                f"{real_path.name}: {real.shape[0]} units after alignment but "
                f"the CCG mask covers {n_units}"
            )

    elif method == "glmcc":
        if not use_surrogates:
            print(
                "ERROR: --method glmcc requires --surrogates flag.\n"
                "  Pre-computed surrogate files are needed for GLMCC validation.\n"
                "  Use --method ccg for on-the-fly surrogate generation instead.",
                file=sys.stderr,
            )
            sys.exit(1)

        real, real_path = _load_adjacency(animal, condition, estimator)
        n_units = real.shape[0]
        p_values, _ = _validate_with_surrogates(
            real, animal, condition, estimator, fdr_q
        )

    else:
        raise ValueError(f"Unknown method: {method}. Use 'ccg' or 'glmcc'.")

    # --- BH-FDR correction (off-diagonal only) ---
    reject = apply_bh_fdr(p_values, fdr_q)

    # --- Write validated adjacency ---
    n_possible = n_units * (n_units - 1)

    # Both methods store the same quantity: the estimator's own coupling weights,
    # masked by whichever significance test was requested. Substituting the CCG
    # peak count here instead made lightON a spike count while ongoing stayed a
    # binary mask, so `node_strength ~ condition` compared encodings.
    validated = real.copy()
    validated[~reject] = 0.0

    out_dir = PROJECT_ROOT / "results" / estimator
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"validated_adj_{animal}_{condition}.csv"
    np.savetxt(str(out_path), validated, delimiter=",", fmt="%.6f")

    # Provenance travels with the matrix. A run summary is a shared file that a
    # later run can overwrite; a sidecar cannot be invalidated by an unrelated run.
    import json
    from datetime import datetime, timezone
    meta_path = out_path.with_suffix(".meta.json")
    meta_path.write_text(json.dumps({
        "encoding": f"{estimator}_weight_masked",
        "method": method,
        "n_units": int(n_units),
        "n_surrogates": int(n_surrogates),
        "fdr_q": fdr_q,
        "fdr_procedure": "benjamini_hochberg_step_up",
        "written_at": datetime.now(timezone.utc).isoformat(),
    }, indent=2))

    # Write p-values for transparency. These are the expensive part of the
    # pipeline; keeping them on disk means the FDR level or the encoding can be
    # revisited without regenerating surrogates.
    p_path = out_dir / f"p_values_{animal}_{condition}.csv"
    np.savetxt(str(p_path), p_values, delimiter=",", fmt="%.6f")

    n_significant = int(reject.sum())
    survival = n_significant / n_possible if n_possible > 0 else 0

    stats = {
        "animal": animal,
        "condition": condition,
        "estimator": estimator,
        "method": method,
        "encoding": f"{estimator}_weight_masked",
        "n_surrogates": n_surrogates,
        "n_units": n_units,
        "n_possible_edges": n_possible,
        "n_significant": n_significant,
        "edge_survival_rate": round(survival, 4),
        "fdr_q": fdr_q,
        "p_values_path": str(p_path),
    }
    return out_path, stats


def main():
    parser = argparse.ArgumentParser(
        description="Edge-level validation with empirical p-values",
        epilog=(
            "Examples:\n"
            "  # CCG method (recommended for lightON):\n"
            "  python validate_edges.py --animal day1 --condition lightON --method ccg\n"
            "  # GLMCC method (requires pre-computed surrogates):\n"
            "  python validate_edges.py --animal night2 --condition lightON --method glmcc --surrogates\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--animal", required=True)
    parser.add_argument("--condition", required=True,
                        choices=["ongoing", "lightON", "ongoing_bis"])
    parser.add_argument("--estimator", default="glmcc", choices=["glmcc", "te"])
    parser.add_argument("--method", default="ccg", choices=["ccg", "glmcc"],
                        help="Validation method: 'ccg' (shuffle-ISI, on-the-fly) "
                             "or 'glmcc' (pre-computed surrogates on disk)")
    parser.add_argument("--surrogates", action="store_true",
                        help="Required when --method glmcc: use pre-computed surrogate files")
    parser.add_argument("--fdr-q", type=float, default=0.05)
    parser.add_argument("--n-surrogates", type=int, default=100,
                        help="Shuffle-ISI surrogates per unit (default: 100)")
    args = parser.parse_args()

    out_path, stats = validate_edges(
        args.animal, args.condition, args.estimator,
        method=args.method,
        fdr_q=args.fdr_q,
        n_surrogates=args.n_surrogates,
        use_surrogates=args.surrogates,
    )
    print(f"\nValidated adjacency: {out_path}")
    for k, v in stats.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
