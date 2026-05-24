#!/usr/bin/env python3
"""Compute surrogate adjacency matrices using cross-correlogram peak (CCG).

Bypasses GLMCC for edge validation by using the CCG peak as a connectivity
measure. This is thousands of times faster than GLMCC and robust to
surrogate spike trains that would cause GLMCC's log(rate) to fail.

Output: results/{estimator}/surrogate_adj/{animal}_{condition}/{surr_id}.csv
"""

import argparse
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
N_WORKERS = 6
WINDOW_MS = 50.0
BIN_MS = 1.0

ANIMALS = ["171019","171207","171208","171213","180110","180111","180131","180221",
           "180228","180302","180419","180420","180423"]
CONDITIONS = ["lightON", "ongoing"]


def ccg_peak_matrix(cell_times: list[np.ndarray]) -> np.ndarray:
    """Compute N×N CCG peak adjacency matrix from list of spike time arrays."""
    n = len(cell_times)
    adj = np.zeros((n, n))
    n_bins = int(2 * WINDOW_MS / BIN_MS) + 1

    for i in range(n):
        a = cell_times[i]
        if len(a) < 2:
            continue
        for j in range(i + 1, n):
            b = cell_times[j]
            if len(b) < 2:
                continue

            hist = np.zeros(n_bins, dtype=int)
            chunk_size = max(1, len(b) // 10)
            for i0 in range(0, len(b), chunk_size):
                chunk = b[i0:i0 + chunk_size]
                diffs = chunk[:, None] - a[None, :]
                mask = np.abs(diffs) <= WINDOW_MS
                valid = diffs[mask]
                if len(valid) > 0:
                    idx = ((valid + WINDOW_MS) / BIN_MS).astype(int)
                    np.add.at(hist, idx, 1)

            peak = int(hist.max())
            adj[i, j] = peak
            adj[j, i] = peak

    return adj


def compute_surrogate_adj(animal: str, condition: str, surr_id: int) -> bool:
    src_dir = PROJECT_ROOT / "data" / "surrogates" / f"{animal}_{condition}" / str(surr_id)
    if not src_dir.exists():
        return False

    out_dir = PROJECT_ROOT / "results" / "glmcc" / "surrogate_adj" / f"{animal}_{condition}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{surr_id}.csv"
    if out_path.exists():
        return True

    cell_files = sorted(src_dir.glob("cell_*.txt"))
    times = []
    for cf in cell_files:
        t = np.atleast_1d(np.loadtxt(str(cf)))
        if len(t) > 0:
            times.append(t * 1000.0)

    if len(times) < 2:
        return False

    adj = ccg_peak_matrix(times)
    np.savetxt(str(out_path), adj, delimiter=",", fmt="%d")
    return True


def batch_run(animal: str, condition: str):
    surr_dir = PROJECT_ROOT / "data" / "surrogates" / f"{animal}_{condition}"
    if not surr_dir.exists():
        print(f"No surrogates for {animal}_{condition}")
        return

    surr_ids = sorted(
        int(d.name) for d in surr_dir.iterdir() if d.is_dir() and d.name.isdigit()
    )
    total = len(surr_ids)
    done = [0]

    def process(sid):
        return compute_surrogate_adj(animal, condition, sid)

    with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {pool.submit(process, sid): sid for sid in surr_ids}
        for f in as_completed(futures):
            if f.result():
                done[0] += 1

    print(f"{animal}_{condition}: {done[0]}/{total} done")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--animal")
    parser.add_argument("--condition", choices=["lightON", "ongoing"])
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()

    global N_WORKERS
    N_WORKERS = args.workers

    if args.all:
        jobs = [(a, c) for a in ANIMALS for c in CONDITIONS]
    elif args.animal and args.condition:
        jobs = [(args.animal, args.condition)]
    else:
        parser.error("Specify --animal + --condition or --all")

    for a, c in jobs:
        batch_run(a, c)


if __name__ == "__main__":
    main()
