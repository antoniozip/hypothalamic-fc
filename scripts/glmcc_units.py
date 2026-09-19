#!/usr/bin/env python3
"""Run GLMCC on second-valued spike files with the correct time base.

GLMCC reads `WIN = 50` and `DELTA = 1` in whatever unit the spike files happen to be
written in, and applies `tau = 4` and its 1-4 delay scan the same way. The repository's
`DATA*/cell*.txt` files, and everything derived from them by `gap_spike_trains.py` and
`jitter_spike_trains.py`, are in **seconds**. Handing those to the binary directly builds
a +/-50 SECOND cross-correlogram at 1 s resolution with a 4 s synaptic time constant --
not a synaptic measurement.

That is what the duty-cycle, jitter and sensitivity controls did until 2026-09-19, because
they called the binary with four arguments and no explicit T. It made every
millisecond-scale manipulation invisible by construction: a +/-25 ms displacement is 1/40
of one bin. The main connectivity matrices were never affected -- `regenerate_glmcc.py`
converts to milliseconds and passes the true duration, which is what this module does too.

Use `run_glmcc()` for any spike directory written in seconds.
"""

from __future__ import annotations

import subprocess
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
GLMCC_BIN = PROJECT_ROOT / "vendor" / "glmcc-c" / "glmcc"


def to_milliseconds(src: Path, dst: Path) -> float:
    """Rewrite second-valued cell files as 0-based milliseconds. Returns the span in seconds."""
    dst.mkdir(parents=True, exist_ok=True)
    cells = sorted(src.glob("cell*.txt"), key=lambda p: int(p.stem[4:]))
    trains = [np.atleast_1d(np.loadtxt(c)) for c in cells]
    nonempty = [t for t in trains if t.size]
    if not nonempty:
        return 0.0
    t0 = min(t.min() for t in nonempty)
    t1 = max(t.max() for t in nonempty)
    for path, t in zip(cells, trains):
        np.savetxt(dst / path.name, np.sort((t - t0) * 1000.0) if t.size else t, fmt="%.3f")
    return float(t1 - t0)


def run_glmcc(work: Path, n_cells: int, mode: str = "exp", method: str = "GLM",
              keep_ms: bool = False) -> tuple[np.ndarray | None, float, str]:
    """Convert `work` to milliseconds, run GLMCC with the true duration, return (W, secs, err)."""
    import shutil

    ms = work.parent / (work.name + "_ms")
    shutil.rmtree(ms, ignore_errors=True)
    dur_s = to_milliseconds(work, ms)
    if dur_s <= 0:
        shutil.rmtree(ms, ignore_errors=True)
        return None, 0.0, "no spikes in window"

    started = time.time()
    proc = subprocess.run(
        [str(GLMCC_BIN), ".", str(n_cells), mode, method, f"{dur_s * 1.01:.1f}"],
        cwd=ms, capture_output=True, text=True)
    elapsed = round(time.time() - started, 1)

    produced = sorted(ms.glob("W_py_*.csv"))
    if proc.returncode != 0 or not produced:
        err = (proc.stderr or "no W_py output").strip()[:180]
        shutil.rmtree(ms, ignore_errors=True)
        return None, elapsed, err
    W = np.loadtxt(produced[0], delimiter=",")
    if not keep_ms:
        shutil.rmtree(ms, ignore_errors=True)
    return W, elapsed, ""
