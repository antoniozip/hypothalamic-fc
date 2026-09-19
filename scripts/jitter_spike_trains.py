#!/usr/bin/env python3
"""Jitter spike times within their observation windows.

Separates fine-timing coupling (synaptic, 1-5 ms lags) from stimulus-locked rate
co-modulation (common drive, ~100 ms-1 s). Displacing each spike uniformly within
+/- JITTER ms destroys structure below that scale while leaving the slower rate
envelope intact: GLMCC uses a +/-50 ms CCG window with 1 ms bins, so a 25 ms jitter
smears any synaptic peak across the window but preserves the 10 s stimulus envelope.

Spikes displaced outside their originating window are dropped rather than clipped,
since clipping piles spikes at window edges and creates artificial synchrony.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def jitter_directory(src: Path, dst: Path, starts: np.ndarray, ends: np.ndarray,
                     jitter_ms: float = 25.0, seed: int = 42) -> dict:
    """Uniformly displace each spike within +/- jitter_ms, constrained to its window.

    Args:
        src: directory of cell*.txt spike files (times in seconds).
        dst: output directory.
        starts, ends: window edges in seconds; a spike stays only if it remains
            inside the window it started in.
        jitter_ms: displacement half-width in milliseconds.
        seed: RNG seed, so the control is reproducible.

    Returns:
        Summary dict with spike retention statistics.
    """
    dst.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    order = np.argsort(starts)
    starts, ends = starts[order], ends[order]
    jitter_s = jitter_ms / 1000.0

    cells = sorted(src.glob("cell*.txt"), key=lambda p: int(p.stem[4:]))
    kept = total = outside = 0
    for path in cells:
        v = np.atleast_1d(np.loadtxt(path))
        total += v.size
        if v.size:
            idx = np.clip(np.searchsorted(starts, v, side="right") - 1, 0, len(starts) - 1)
            shifted = v + rng.uniform(-jitter_s, jitter_s, size=v.size)
            # Keep only spikes still inside the window they came from.
            inside = (shifted >= starts[idx]) & (shifted <= ends[idx])
            outside += int((~inside).sum())
            w = np.sort(shifted[inside])
        else:
            w = v
        kept += w.size
        np.savetxt(dst / path.name, w, fmt="%.6f")

    return {
        "spikes_kept": kept, "spikes_total": total,
        "dropped_at_edges": outside,
        "retention_pct": 100 * kept / max(total, 1),
        "jitter_ms": jitter_ms, "seed": seed,
    }
