#!/usr/bin/env python3
"""Generate dither-spike-train surrogates for GLMCC/TE significance testing.

Spike times in .mat files are in milliseconds; converted to seconds for Neo/elephant.
Uses cell-specific t_stop from Units['ts'] to match original pipeline.

Usage:
    python generate_surrogates.py --animal 171019 --condition lightON
    python generate_surrogates.py --animal 180131 --condition ongoing --n-surrogates 100 --dt-ms 0.1
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from neo.core import SpikeTrain
from quantities import s, ms
from elephant.spike_train_surrogates import surrogates


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

TS_KEY_MAP = {
    "lightON": "evoked_ts",
    "ongoing": "ongoing_ts",
    "evoked": "evoked_ts",
}


def generate_surrogates(
    animal: str,
    condition: str,
    n_surrogates: int = 100,
    dt_ms: float = 0.1,
) -> int:
    mat_path = PROJECT_ROOT / f"DATA{animal}_{condition}.mat"
    if not mat_path.exists():
        if condition == "lightON":
            mat_path = PROJECT_ROOT / f"DATA{animal}_lightON.mat"
        elif condition in ("ongoing", "ongoing_bis"):
            candidates = sorted(PROJECT_ROOT.glob(f"DATA{animal}*ongoing*.mat"))
            if candidates:
                mat_path = candidates[0]
    if not mat_path.exists():
        raise FileNotFoundError(f"No .mat file for animal={animal}, condition={condition}")

    raw = loadmat(str(mat_path), simplify_cells=False)
    ts_key = TS_KEY_MAP.get(condition, "evoked_ts")
    if ts_key not in raw:
        alt_keys = [k for k in raw if k.endswith("_ts")]
        ts_key = alt_keys[0] if alt_keys else "evoked_ts"
    spike_times = raw[ts_key]
    units_s = raw["Units"]

    n_units = spike_times.shape[1]
    out_dir = PROJECT_ROOT / "data" / "surrogates" / f"{animal}_{condition}"
    out_dir.mkdir(parents=True, exist_ok=True)

    for cell_id in range(n_units):
        cell_ts = spike_times[0, cell_id][0]
        if len(cell_ts) == 0:
            continue

        cell_spikes_s = cell_ts / 1e3

        unit_ts = units_s["ts"][0, cell_id]
        t_stop = float(np.atleast_1d(unit_ts).flatten()[-1]) + 1.0 if len(unit_ts) > 0 else cell_spikes_s[-1] + 1.0

        try:
            st = SpikeTrain(cell_spikes_s * s, t_stop=t_stop)
            all_surr = surrogates(
                st,
                n_surrogates=n_surrogates,
                method="dither_spike_train",
                dt=dt_ms * ms,
            )
        except Exception as exc:
            print(f"  WARNING cell {cell_id}: surrogate failed ({exc}), skipping", flush=True)
            continue

        for surr_id in range(n_surrogates):
            surr_dir = out_dir / str(surr_id)
            surr_dir.mkdir(parents=True, exist_ok=True)
            out_path = surr_dir / f"cell_{cell_id}.txt"
            with open(out_path, "w") as f:
                for spike in all_surr[surr_id]:
                    print(float(spike.magnitude), file=f)

    return n_surrogates


def main():
    parser = argparse.ArgumentParser(description="Generate dither spike-train surrogates")
    parser.add_argument("--animal", required=True, help="Animal ID (e.g. 171019)")
    parser.add_argument(
        "--condition", required=True,
        choices=["ongoing", "lightON", "evoked", "ongoing_bis"],
        help="Recording condition",
    )
    parser.add_argument("--n-surrogates", type=int, default=100)
    parser.add_argument("--dt-ms", type=float, default=0.1)
    args = parser.parse_args()

    n = generate_surrogates(args.animal, args.condition, args.n_surrogates, args.dt_ms)
    print(f"Generated {n} surrogate sets for {args.animal}_{args.condition}")


if __name__ == "__main__":
    main()
