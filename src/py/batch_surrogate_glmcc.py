#!/usr/bin/env python3
"""Batch-run GLMCC on surrogate spike-train sets to generate surrogate adjacency matrices.

For each animal x condition x surr_id:
  1. Reads cell spike trains from data/surrogates/{animal}_{condition}/{surr_id}/
  2. Shifts spike times to start from 0 (GLMCC expects 0-based ms times)
  3. Writes temporary cell files
  4. Calls Est_Data.py (Kobayashi et al. 2019 GLMCC implementation)
  5. Collects the W_py_T.csv output
  6. Saves to results/{estimator}/surrogate_adj/{animal}_{condition}/{surr_id}.csv

Usage:
    python src/py/batch_surrogate_glmcc.py --animal 171019 --condition lightON
    python src/py/batch_surrogate_glmcc.py --all
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
GLMCC_DIR = Path(os.path.expanduser("~/Documents/documents_old/HypoNetworkAntonio/GLMCC"))
EST_DATA = GLMCC_DIR / "Est_Data.py"
N_WORKERS = 6

ANIMALS = ["171019","171207","171208","171213","180110","180111","180131","180221",
           "180228","180302","180419","180420","180423"]
CONDITIONS = ["lightON", "ongoing"]


def process_surrogate_set(animal: str, condition: str, surr_id: int) -> bool:
    src_dir = PROJECT_ROOT / "data" / "surrogates" / f"{animal}_{condition}" / str(surr_id)
    if not src_dir.exists():
        print(f"SKIP: {animal}_{condition}/{surr_id} — no source")
        return False

    cell_files = sorted(src_dir.glob("cell_*.txt"))
    n_units = len(cell_files)
    if n_units == 0:
        return False

    out_dir = PROJECT_ROOT / "results" / "glmcc" / "surrogate_adj" / f"{animal}_{condition}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{surr_id}.csv"
    if out_path.exists():
        return True

    with tempfile.TemporaryDirectory(prefix=f"glmcc_{animal}_{condition}_{surr_id}_") as tmp_dir:
        tmp_path = Path(tmp_dir)

        max_time = 0.0
        min_time = float("inf")
        all_times = []
        for cell_file in cell_files:
            times = np.atleast_1d(np.loadtxt(str(cell_file)))
            if len(times) == 0:
                continue
            all_times.append(times)
            min_time = min(min_time, float(times[0]))
            max_time = max(max_time, float(times[-1]))

        if not all_times:
            return False

        for idx, (cell_file, times) in enumerate(zip(cell_files, all_times)):
            times_ms = (times - min_time) * 1000.0
            np.savetxt(tmp_path / f"cell{idx}.txt", times_ms, fmt="%.6f")

        t_val = int(np.ceil((max_time - min_time)))
        if t_val < 10:
            t_val = 5400

        for mod_file in ["glmcc.py", "glmcc_fitting.py"]:
            src_mod = GLMCC_DIR / mod_file
            if src_mod.exists():
                shutil.copy2(str(src_mod), str(tmp_path / mod_file))

        try:
            result = subprocess.run(
                [sys.executable, str(EST_DATA), str(tmp_path), str(n_units), "exp", "GLM"],
                cwd=str(tmp_path),
                capture_output=True,
                text=True,
                timeout=600,
            )
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT: {animal}_{condition}/{surr_id} (>{n_units}s)")
            return False
        except Exception as e:
            print(f"  ERROR: {animal}_{condition}/{surr_id} — {e}")
            return False

        w_file = tmp_path / f"W_py_{t_val}.csv"
        if not w_file.exists():
            alt = sorted(tmp_path.glob("W_py_*.csv"))
            if alt:
                w_file = alt[0]
            else:
                return False

        try:
            adj = np.loadtxt(str(w_file), delimiter=",")
            np.savetxt(str(out_path), adj, delimiter=",")
        except Exception:
            return False

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
        if done[0] % 10 == 0:
            print(f"  {animal}_{condition}: {done[0]}/{total}", flush=True)
        return process_surrogate_set(animal, condition, sid)

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
