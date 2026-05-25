# -*- coding: utf-8 -*-
"""Transfer Entropy connectivity inference for hypothalamic spike trains.

Bins spike trains from data/processed/ into a discrete time series,
then computes pairwise transfer entropy using a pure numpy implementation.

Usage:
    python run_transfer_entropy.py --animal day1 --condition ongoing
    python run_transfer_entropy.py --animal day2 --condition lightON --k 10 --binsize-ms 5
"""

import argparse
from pathlib import Path

import numpy as np
import scipy.io as sio

from transfer_entropy_np import transfer_entropy, transfer_entropy_matrix


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "te"


def _load_spike_train_from_file(filepath: Path) -> np.ndarray:
    with open(filepath) as f:
        raw = f.read().strip()
    rows = raw.split("\n")
    times = []
    for row in rows:
        parts = row.strip().split()
        if parts:
            times.append([float(x) for x in parts])
    return times


def _load_spike_train_from_mat(matpath: Path, condition: str) -> np.ndarray:
    var = sio.loadmat(str(matpath))
    key = "evoked_st" if "evoked_st" in var else "ongoing_st"
    st = var[key]
    if hasattr(st, "toarray"):
        st = st.toarray()
    return st


def _bin_spike_trains(
    spike_times: list[list[float]],
    binsize_ms: float,
    t_start: float | None = None,
    t_stop: float | None = None,
) -> np.ndarray:
    if t_start is None:
        t_start = min(
            min(st) for st in spike_times if len(st) > 0
        )
    if t_stop is None:
        t_stop = max(
            max(st) for st in spike_times if len(st) > 0
        )
    binsize_s = binsize_ms / 1000.0
    n_bins = int(np.ceil((t_stop - t_start) / binsize_s)) + 1
    n_units = len(spike_times)
    binned = np.zeros((n_units, n_bins), dtype=int)
    for i, st in enumerate(spike_times):
        for t in st:
            idx = int((t - t_start) / binsize_s)
            if 0 <= idx < n_bins:
                binned[i, idx] = 1
    return binned


def run_transfer_entropy(
    animal: str,
    condition: str,
    k: int = 10,
    binsize_ms: float = 5.0,
) -> Path:
    txt_path = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / f"spiketrains_DATA{animal}_{condition}.txt"
    )
    if txt_path.exists():
        spike_times = _load_spike_train_from_file(txt_path)
    else:
        mat_path = PROJECT_ROOT / "data" / "raw" / f"DATA{animal}_{condition}.mat"
        if not mat_path.exists():
            mat_path = PROJECT_ROOT / f"DATA{animal}_{condition}.mat"
        if mat_path.exists():
            spike_times_list = _load_spike_train_from_mat(mat_path, condition)
            spike_times = [list(spike_times_list[i, :]) for i in range(spike_times_list.shape[0])]
            spike_times = [st for st in spike_times if len(st) > 0]
        else:
            raise FileNotFoundError(
                f"No spike data for animal={animal}, condition={condition}"
            )

    binned = _bin_spike_trains(spike_times, binsize_ms)
    n_units = binned.shape[0]

    te = transfer_entropy_matrix(binned, k=k)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RESULTS_DIR / f"te_{animal}_{condition}_k{k}.csv"
    np.savetxt(str(out_path), te, delimiter=",")
    return out_path


def main():
    parser = argparse.ArgumentParser(description="Transfer Entropy connectivity inference")
    parser.add_argument("--animal", required=True)
    parser.add_argument("--condition", required=True,
                        choices=["ongoing", "lightON", "ongoing_bis", "evoked"])
    parser.add_argument("--k", type=int, default=10,
                        help="History length for TE (default: 10)")
    parser.add_argument("--binsize-ms", type=float, default=5.0,
                        help="Bin size in ms (default: 5)")
    args = parser.parse_args()

    out_path = run_transfer_entropy(args.animal, args.condition, args.k, args.binsize_ms)
    print(f"TE matrix written to {out_path}")


if __name__ == "__main__":
    main()
