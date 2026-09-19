#!/usr/bin/env python3
"""Characterise WHERE cross-correlogram peaks sit in lag, not just how high they are.

The jitter control shows whether structure lives below ~25 ms, but cannot separate synaptic
coupling from tightly-locked common drive. Peak lag can:

  common drive  -> CCG peaks at lag ~= 0 and is roughly symmetric about it
  synaptic      -> CCG peaks at a short non-zero lag (~1-5 ms) and is asymmetric

Reported per condition:
  zero_lag_pct  share of pairs whose peak sits at |lag| <= 1 ms
  median_abs_lag   median |peak lag| in ms
  asym_index    median |sum(CCG, lag>0) - sum(CCG, lag<0)| / total mass, per pair
                (0 = perfectly symmetric, 1 = all mass on one side)

Usage:
    python scripts/ccg_lag_structure.py --animals day1 night3 --conditions lightON ongoing_aligned
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from gap_spike_trains import gap_directory_windows, stimulus_windows  # noqa: E402

WINDOW_MS = 50.0
BIN_MS = 1.0
N_BINS = int(2 * WINDOW_MS / BIN_MS) + 1

ANIMAL_MAP = {
    "day1": "171019", "day2": "180131", "day3": "180302",
    "day4": "180419", "day5": "180420", "day6": "180423",
    "night1": "171207", "night2": "171208", "night3": "171213",
    "night4": "180110", "night5": "180111", "night6": "180221", "night7": "180228",
}


def ccg_histogram(a_s: np.ndarray, b_s: np.ndarray) -> np.ndarray:
    """Cross-correlogram of two spike trains (seconds in, ms lags out)."""
    a_ms, b_ms = a_s * 1000.0, b_s * 1000.0
    hist = np.zeros(N_BINS, dtype=np.int64)
    if len(a_ms) < 2 or len(b_ms) < 2:
        return hist
    for i0 in range(0, len(b_ms), 200_000):
        sub = b_ms[i0:i0 + 200_000]
        lo = np.searchsorted(a_ms, sub - WINDOW_MS, side="left")
        hi = np.searchsorted(a_ms, sub + WINDOW_MS, side="right")
        counts = hi - lo
        total = int(counts.sum())
        if total == 0:
            continue
        ends = np.cumsum(counts)
        offsets = np.arange(total) - np.repeat(ends - counts, counts)
        a_idx = np.repeat(lo, counts) + offsets
        b_idx = np.repeat(np.arange(len(sub)), counts)
        diffs = sub[b_idx] - a_ms[a_idx]
        diffs = diffs[np.abs(diffs) <= WINDOW_MS]
        if diffs.size:
            idx = ((diffs + WINDOW_MS) / BIN_MS).astype(np.int64)
            hist += np.bincount(idx, minlength=N_BINS)[:N_BINS]
    return hist


def load_cells(path: Path) -> list[np.ndarray]:
    return [np.atleast_1d(np.loadtxt(p))
            for p in sorted(path.glob("cell*.txt"), key=lambda q: int(q.stem[4:]))]


def analyse(trains: list[np.ndarray], min_count: int = 30, z_thresh: float = 5.0) -> dict:
    """Lag structure over all pairs, and over the subset with a prominent peak.

    Prominence is the peak's z-score against the CCG's own flanks (|lag| in 30-50 ms),
    which is the usual way to tell a real short-latency feature from the noisiest bin.
    A peak lag distribution with median |lag| ~= 25 ms is what uniform noise over a
    +/-50 ms window produces, so it means "no peak", not "a peak at 25 ms".
    """
    lag_axis = np.arange(N_BINS) - (N_BINS // 2)
    centre = N_BINS // 2
    flank = np.abs(lag_axis) >= 30

    lags, asyms, n_pairs = [], [], 0
    strong_lags, n_strong = [], 0
    for i in range(len(trains)):
        for j in range(i + 1, len(trains)):
            h = ccg_histogram(trains[i], trains[j])
            if h.sum() < min_count:
                continue
            n_pairs += 1
            k = int(np.argmax(h))
            lags.append(int(lag_axis[k]))
            pos, neg = h[centre + 1:].sum(), h[:centre].sum()
            tot = pos + neg
            asyms.append(abs(pos - neg) / tot if tot else 0.0)

            base, sd = h[flank].mean(), h[flank].std()
            if sd > 0 and (h[k] - base) / sd >= z_thresh:
                n_strong += 1
                strong_lags.append(int(lag_axis[k]))

    if not lags:
        return {"n_pairs": 0}
    lags_arr = np.array(lags)
    out = {
        "n_pairs": n_pairs,
        "zero_lag_pct": round(100 * float((np.abs(lags_arr) <= 1).mean()), 1),
        "median_abs_lag": round(float(np.median(np.abs(lags_arr))), 1),
        "mean_lag": round(float(lags_arr.mean()), 2),
        "asym_index": round(float(np.median(asyms)), 4),
        "n_strong": n_strong,
        "strong_pct": round(100 * n_strong / n_pairs, 1),
    }
    if strong_lags:
        sa = np.array(strong_lags)
        out.update({
            "strong_zero_lag_pct": round(100 * float((np.abs(sa) <= 1).mean()), 1),
            "strong_median_abs_lag": round(float(np.median(np.abs(sa))), 1),
            "strong_short_lag_pct": round(100 * float((np.abs(sa) <= 5).mean()), 1),
        })
    else:
        out.update({"strong_zero_lag_pct": "", "strong_median_abs_lag": "",
                    "strong_short_lag_pct": ""})
    return out


def get_trains(animal: str, condition: str, workroot: Path) -> list[np.ndarray] | None:
    old = ANIMAL_MAP[animal]
    if condition == "lightON":
        return load_cells(PROJECT_ROOT / f"DATA{old}")
    if condition == "ongoing_aligned":
        mat = next((c for c in (PROJECT_ROOT / f"DATA{old}_ongoing.mat",
                                PROJECT_ROOT / f"DATA{old}_lightON.mat") if c.exists()), None)
        if mat is None:
            return None
        starts, ends = stimulus_windows(mat)
        dur = ends - starts
        work = workroot / f"{animal}_aligned"
        gap_directory_windows(PROJECT_ROOT / f"DATA{old}_ongoing", work,
                              starts + 45.0, starts + 45.0 + dur)
        return load_cells(work)
    raise ValueError(condition)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--animals", nargs="+", default=list(ANIMAL_MAP))
    parser.add_argument("--conditions", nargs="+", default=["lightON", "ongoing_aligned"])
    parser.add_argument("--workroot", type=Path, default=Path("/tmp/claude-1000/ccg_lag"))
    args = parser.parse_args()
    args.workroot.mkdir(parents=True, exist_ok=True)

    rows = []
    print(f"{'animal':8s} {'condition':16s} {'pairs':>7s} {'med|lag|':>9s} "
          f"{'strong%':>8s} {'s:0lag%':>8s} {'s:<=5ms%':>9s} {'s:med|lag|':>11s}")
    for animal in args.animals:
        for cond in args.conditions:
            trains = get_trains(animal, cond, args.workroot)
            if trains is None:
                continue
            st = analyse(trains)
            if not st.get("n_pairs"):
                continue
            rows.append({"animal": animal, "condition": cond, **st})
            print(f"{animal:8s} {cond:16s} {st['n_pairs']:7d} {st['median_abs_lag']:9.1f} "
                  f"{st['strong_pct']:8.1f} {str(st['strong_zero_lag_pct']):>8s} "
                  f"{str(st['strong_short_lag_pct']):>9s} {str(st['strong_median_abs_lag']):>11s}",
                  flush=True)

    out = PROJECT_ROOT / "results" / "ccg_lag_structure.csv"
    if rows:
        with out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\nWrote {out}")
        for cond in args.conditions:
            sel = [r for r in rows if r["condition"] == cond]
            if sel:
                strong = [r for r in sel if r.get("strong_short_lag_pct") != ""]
                print(f"  {cond:16s} all pairs: median |lag| "
                      f"{statistics.median(r['median_abs_lag'] for r in sel):4.1f} ms "
                      f"(uniform noise over +/-50 ms gives 25.0)")
                if strong:
                    print(f"  {'':16s} prominent peaks: "
                          f"{statistics.median(r['strong_pct'] for r in strong):4.1f}% of pairs, "
                          f"of those {statistics.median(r['strong_short_lag_pct'] for r in strong):4.1f}% "
                          f"at |lag|<=5 ms, "
                          f"{statistics.median(r['strong_zero_lag_pct'] for r in strong):4.1f}% at |lag|<=1 ms")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
