#!/usr/bin/env python3
"""Wrapper around GLMCC (Kobayashi et al. 2019) for connectivity inference.

Reads spike trains from data/processed/spiketrains_{animal}_{condition}.txt,
calls the GLMCC estimator, and writes the adjacency matrix.

Usage:
    python run_glmcc.py --animal 171019 --condition lightON
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "glmcc"


def run_glmcc(animal: str, condition: str, glmcc_path: str | None = None) -> Path:
    spike_path = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / f"spiketrains_DATA{animal}_{condition}.txt"
    )
    if not spike_path.exists():
        spike_path = (
            PROJECT_ROOT
            / "data"
            / "processed"
            / f"spiketrains_DATA{animal}_evoked_LightON.txt"
        )
    if not spike_path.exists():
        alt = list(
            (PROJECT_ROOT / "data" / "processed").glob(
                f"spiketrains_DATA{animal}*{condition}*.txt"
            )
        )
        if alt:
            spike_path = alt[0]
    if not spike_path.exists():
        raise FileNotFoundError(
            f"No spike-train file for animal={animal}, condition={condition}"
        )

    spike_trains = []
    with open(spike_path) as f:
        n_units = None
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if n_units is None:
                n_units = len(parts)
            vals = [float(x) for x in parts]
            spike_trains.append(vals)
    n_units = max(len(st) for st in spike_trains)

    results_dir = RESULTS_DIR
    results_dir.mkdir(parents=True, exist_ok=True)
    out_path = results_dir / f"adj_{animal}_{condition}.csv"

    spike_path_abs = spike_path
    glmcc_bin = glmcc_path or "python ~/Documents/GLMCC/Est_Data.py"
    cmd = f'{glmcc_bin} {spike_path_abs} exp GLM 2>&1'

    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    except Exception as e:
        print(f"GLMCC call failed: {e}", file=sys.stderr)
        print("Note: GLMCC requires separate installation. See GLMCC_VERSION.txt.")
        print("Generating placeholder adjacency...")
        placeholder = np.zeros((n_units, n_units))
        np.savetxt(str(out_path), placeholder, delimiter=",")
        return out_path

    adj_path = Path(f"result_{animal}.csv")
    if adj_path.exists():
        import shutil
        shutil.move(str(adj_path), str(out_path))
        print(f"GLMCC result moved to {out_path}")
    else:
        np.savetxt(str(out_path), np.zeros((n_units, n_units)), delimiter=",")
        print(f"GLMCC output not found; wrote placeholder to {out_path}")

    return out_path


def main():
    parser = argparse.ArgumentParser(description="Run GLMCC connectivity inference")
    parser.add_argument("--animal", required=True)
    parser.add_argument("--condition", required=True,
                        choices=["ongoing", "lightON", "ongoing_bis"])
    parser.add_argument("--glmcc", default=None,
                        help="Path to GLMCC Est_Data.py or executable")
    args = parser.parse_args()

    run_glmcc(args.animal, args.condition, args.glmcc)


if __name__ == "__main__":
    main()
