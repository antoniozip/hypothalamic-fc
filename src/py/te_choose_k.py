#!/usr/bin/env python3
"""Pick TE history length k by AIC/BIC on a held-out animal.

Sweeps k ∈ {1, 3, 5, 10, 15, 20} and reports AIC/BIC per k.
"""

import argparse
from pathlib import Path

import numpy as np
import pyinform

from run_transfer_entropy import (
    _load_spike_train_from_file,
    _load_spike_train_from_mat,
    _bin_spike_trains,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


def compute_te_aic_bic(
    animal: str,
    condition: str,
    k_values: list[int],
    binsize_ms: float = 5.0,
) -> dict:
    txt_path = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / f"spiketrains_DATA{animal}_{condition}.txt"
    )
    if txt_path.exists():
        spike_times = _load_spike_train_from_file(txt_path)
    else:
        mat_path = PROJECT_ROOT / f"DATA{animal}_{condition}.mat"
        st_list = _load_spike_train_from_mat(mat_path, condition)
        n_units = st_list.shape[0]
        spike_times = [list(st_list[i, :]) for i in range(n_units)]
        spike_times = [st for st in spike_times if len(st) > 0]

    binned = _bin_spike_trains(spike_times, binsize_ms)
    n_units = binned.shape[0]
    n_bins = binned.shape[1]

    results = {"k": [], "aic": [], "bic": []}
    for k in k_values:
        aic_total = 0.0
        bic_total = 0.0
        n_pairs = 0
        for i in range(n_units):
            for j in range(n_units):
                if i == j:
                    continue
                te = pyinform.transfer_entropy(binned[i, :], binned[j, :], k=k)
                n_params = k * 2
                aic = -2 * te * n_bins + 2 * n_params
                bic = -2 * te * n_bins + n_params * np.log(n_bins)
                aic_total += aic
                bic_total += bic
                n_pairs += 1
        results["k"].append(k)
        results["aic"].append(aic_total / n_pairs)
        results["bic"].append(bic_total / n_pairs)
        print(f"  k={k}: mean AIC={results['aic'][-1]:.2f}, mean BIC={results['bic'][-1]:.2f}")

    best_aic = results["k"][np.argmin(results["aic"])]
    best_bic = results["k"][np.argmin(results["bic"])]
    print(f"\nBest k by AIC: {best_aic}, by BIC: {best_bic}")

    return results


def main():
    parser = argparse.ArgumentParser(description="TE history-length selection")
    parser.add_argument("--animal", required=True)
    parser.add_argument("--condition", default="ongoing",
                        choices=["ongoing", "lightON", "evoked"])
    parser.add_argument("--binsize-ms", type=float, default=5.0)
    parser.add_argument("--k-values", type=int, nargs="+",
                        default=[1, 3, 5, 10, 15, 20])
    args = parser.parse_args()

    compute_te_aic_bic(args.animal, args.condition, args.k_values, args.binsize_ms)


if __name__ == "__main__":
    main()
