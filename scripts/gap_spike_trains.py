#!/usr/bin/env python3
"""Impose a stimulus-window duty cycle on a continuous spike-train directory.

Used to test whether GLMCC's sign inversion between conditions is caused by the
lightON epoch being a ~10% duty-cycle concatenation of stimulus windows rather
than a continuous recording. Gapping an `ongoing` recording this way, with the
within-window firing rate preserved, reproduces the lightON sign signature.

Usage:
    python scripts/gap_spike_trains.py --src DATA171213_ongoing --dst /tmp/gapped \
        --on 9.59 --period 89.9
"""

import argparse
import os
from glob import glob
from pathlib import Path

import numpy as np


def measure_windows(src: Path, gap_s: float = 1.0) -> tuple[float, float, float]:
    """Return (median window length, median onset period, duty cycle) of a directory."""
    pooled = []
    for path in sorted(src.glob("cell*.txt"), key=lambda p: int(p.stem[4:])):
        v = np.atleast_1d(np.loadtxt(path))
        if v.size:
            pooled.append(v)
    allt = np.sort(np.concatenate(pooled))
    idx = np.flatnonzero(np.diff(allt) > gap_s)
    starts = np.concatenate([[allt[0]], allt[idx + 1]])
    ends = np.concatenate([allt[idx], [allt[-1]]])
    lengths = ends - starts
    duty = 100 * lengths.sum() / (allt[-1] - allt[0])
    return float(np.median(lengths)), float(np.median(np.diff(starts))), float(duty)


def gap_directory(src: Path, dst: Path, on_s: float, period_s: float) -> dict:
    """Keep only spikes inside periodic windows; preserve the original time base."""
    dst.mkdir(parents=True, exist_ok=True)
    cells = sorted(src.glob("cell*.txt"), key=lambda p: int(p.stem[4:]))

    pooled = [np.atleast_1d(np.loadtxt(c)) for c in cells]
    allt = np.concatenate([v for v in pooled if v.size])
    t0, t1 = allt.min(), allt.max()
    starts = np.arange(t0, t1, period_s)
    ends = starts + on_s

    kept = total = 0
    for path, v in zip(cells, pooled):
        total += v.size
        if v.size:
            i = np.searchsorted(starts, v, side="right") - 1
            inside = (i >= 0) & (v <= ends[np.clip(i, 0, len(ends) - 1)])
            w = v[inside]
        else:
            w = v
        kept += w.size
        np.savetxt(dst / path.name, w, fmt="%.6f")

    observed = len(starts) * on_s
    return {
        "n_windows": len(starts), "spikes_kept": kept, "spikes_total": total,
        "duty_cycle_pct": 100 * observed / (t1 - t0),
        "rate_hz_per_neuron": kept / observed / len(cells),
    }


def stimulus_windows(mat_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read true light-pulse edges (seconds) from a .mat file's `mask_stimulus`.

    The mask is a per-millisecond uint8 flag over the whole session.
    """
    from scipy.io import loadmat

    m = loadmat(str(mat_path), simplify_cells=False)
    if "mask_stimulus" not in m:
        raise KeyError(f"{mat_path.name} has no mask_stimulus")
    on = m["mask_stimulus"].flatten() > 0
    trans = np.diff(on.astype(np.int8))
    starts = np.flatnonzero(trans == 1) / 1000.0
    ends = np.flatnonzero(trans == -1) / 1000.0
    k = min(len(starts), len(ends))
    return starts[:k], ends[:k]


def gap_directory_windows(src: Path, dst: Path, starts: np.ndarray,
                          ends: np.ndarray) -> dict:
    """Keep only spikes inside the given explicit windows; preserve the time base.

    Unlike gap_directory(), which imposes a regular fitted period, this takes real
    window edges so the mask can be phase-aligned to the actual stimulus.
    """
    dst.mkdir(parents=True, exist_ok=True)
    cells = sorted(src.glob("cell*.txt"), key=lambda p: int(p.stem[4:]))
    order = np.argsort(starts)
    starts, ends = starts[order], ends[order]

    kept = total = 0
    for path in cells:
        v = np.atleast_1d(np.loadtxt(path))
        total += v.size
        if v.size:
            i = np.searchsorted(starts, v, side="right") - 1
            inside = (i >= 0) & (v <= ends[np.clip(i, 0, len(ends) - 1)])
            w = v[inside]
        else:
            w = v
        kept += w.size
        np.savetxt(dst / path.name, w, fmt="%.6f")

    observed = float((ends - starts).sum())
    epoch = float(ends[-1] - starts[0])
    return {
        "n_windows": len(starts), "spikes_kept": kept, "spikes_total": total,
        "observed_s": observed, "epoch_s": epoch,
        "duty_cycle_pct": 100 * observed / epoch,
        "rate_hz_per_neuron": kept / max(observed, 1) / max(len(cells), 1),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", type=Path)
    parser.add_argument("--dst", type=Path)
    parser.add_argument("--on", type=float, default=9.59, help="window length (s)")
    parser.add_argument("--period", type=float, default=89.9, help="onset period (s)")
    parser.add_argument("--measure", type=Path,
                        help="instead of gapping, report the window structure of this directory")
    args = parser.parse_args()

    if args.measure:
        on, period, duty = measure_windows(args.measure)
        print(f"{args.measure}: window={on:.2f}s period={period:.1f}s duty={duty:.1f}%")
        return 0

    if not args.src or not args.dst:
        parser.error("--src and --dst are required unless --measure is given")
    stats = gap_directory(args.src, args.dst, args.on, args.period)
    for key, value in stats.items():
        print(f"  {key}: {value:,.3f}" if isinstance(value, float) else f"  {key}: {value:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
