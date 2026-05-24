#!/usr/bin/env python3
"""Edge validation using CCG peak (cross-correlogram) as connectivity measure.

Generates shuffle_isis surrogates on-the-fly (preserves ISI distribution,
breaks cross-neuron correlations). Compares real vs surrogate CCG peaks
for empirical p-values + BH-FDR.

This is the PREFERRED validation method for lightON conditions.
For ongoing conditions or when pre-computed GLMCC surrogates exist,
use validate_edges.py --method glmcc.

Functions are importable — validate_edges.py calls into this module.

Usage:
    python validate_ccg.py --animal 171019 --condition lightON
    python validate_ccg.py --all
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from neo.core import SpikeTrain
from quantities import s, ms
from elephant.spike_train_surrogates import surrogates


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
N_WORKERS = 6
N_WORKERS_INNER = 4
WINDOW_MS = 50.0
BIN_MS = 1.0
N_SURROGATES = 100
N_BINS = int(2 * WINDOW_MS / BIN_MS) + 1

ANIMALS = ["171019","171207","171208","171213","180110","180111","180131","180221",
           "180228","180302","180419","180420","180423"]
CONDITIONS = ["lightON"]  # ongoing: use heuristic validation (too dense for CCG)

TS_KEY_MAP = {"lightON": "evoked_ts", "ongoing": "ongoing_ts", "evoked": "evoked_ts"}


# ---------------------------------------------------------------------------
# Public API — importable by validate_edges.py and other modules
# ---------------------------------------------------------------------------

def load_real_spike_times(animal: str, condition: str) -> list[np.ndarray]:
    """Load spike times (in ms) from .mat file for a given animal+condition."""
    mat_path = PROJECT_ROOT / f"DATA{animal}_{condition}.mat"
    if not mat_path.exists():
        if condition == "ongoing":
            alts = sorted(PROJECT_ROOT.glob(f"DATA{animal}*ongoing*.mat"))
            if alts:
                mat_path = alts[0]
    if not mat_path.exists():
        raise FileNotFoundError(f"No .mat for {animal}_{condition}")

    m = loadmat(str(mat_path), simplify_cells=False)
    ts_key = TS_KEY_MAP.get(condition, "evoked_ts")
    if ts_key not in m:
        ts_key = [k for k in m if k.endswith("_ts")][0]

    ts = m[ts_key]
    n = ts.shape[1]
    times = []
    for i in range(n):
        t = np.atleast_1d(ts[0, i][0])
        if len(t) > 0:
            times.append(t.astype(float))  # .mat values are in ms
    return times


def _ccg_peak_fast(a_ms: np.ndarray, b_ms: np.ndarray) -> float:
    """Fast CCG peak using binned sliding window. O(N) instead of FFT's O(N log N)."""
    t_min = min(a_ms[0], b_ms[0])
    t_max = max(a_ms[-1], b_ms[-1])
    dur = int(np.ceil(t_max - t_min)) + 1
    bin_a = np.zeros(dur, dtype=np.float32)
    bin_b = np.zeros(dur, dtype=np.float32)
    np.add.at(bin_a, np.floor(a_ms - t_min).astype(int), 1)
    np.add.at(bin_b, np.floor(b_ms - t_min).astype(int), 1)
    hw = int(WINDOW_MS / BIN_MS)
    shifts = np.arange(-hw, hw + 1)
    best = 0
    for s in shifts:
        if s < 0:
            c = np.dot(bin_a[-s:], bin_b[:s])
        elif s > 0:
            c = np.dot(bin_a[:-s], bin_b[s:])
        else:
            c = np.dot(bin_a, bin_b)
        if c > best:
            best = c
    return int(best)


def ccg_peak_one_pair(a_ms: np.ndarray, b_ms: np.ndarray) -> float:
    """CCG peak for a single neuron pair. Uses fast method for large trains."""
    if len(a_ms) < 2 or len(b_ms) < 2:
        return 0.0
    if len(a_ms) > 100000 or len(b_ms) > 100000:
        return _ccg_peak_fast(a_ms, b_ms)
    n_bins = int(2 * WINDOW_MS / BIN_MS) + 1
    hist = np.zeros(n_bins, dtype=int)
    chunk = max(1, len(b_ms) // 20)
    for i0 in range(0, len(b_ms), chunk):
        sub = b_ms[i0:i0 + chunk]
        diffs = sub[:, None] - a_ms[None, :]
        mask = np.abs(diffs) <= WINDOW_MS
        valid = diffs[mask]
        if len(valid) > 0:
            idx = ((valid + WINDOW_MS) / BIN_MS).astype(int)
            np.add.at(hist, idx, 1)
    return int(hist.max())


def ccg_peak_matrix(times_ms: list[np.ndarray]) -> np.ndarray:
    """Full NxN CCG peak matrix. Inputs are spike times in ms."""
    n = len(times_ms)
    adj = np.zeros((n, n))
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]

    def process_pair(pair):
        i, j = pair
        return i, j, ccg_peak_one_pair(times_ms[i], times_ms[j])

    with ThreadPoolExecutor(max_workers=N_WORKERS_INNER) as pool:
        for i, j, p in pool.map(process_pair, pairs):
            adj[i, j] = adj[j, i] = p
    return adj


def compute_empirical_p_values(
    real_times: list[np.ndarray],
    real_peaks: np.ndarray,
    n_surrogates: int = 100,
) -> tuple[np.ndarray, int]:
    """Compute empirical p-values via shuffle-ISI surrogates.

    Returns (p_values, n_units) where p_values[i,i] = 1.0 (diagonal ignored).
    """
    n_units = len(real_times)
    count_exceed = np.zeros((n_units, n_units), dtype=int)
    total_surr = 0

    # Build Neo SpikeTrain objects once
    spike_trains = []
    for t in real_times:
        st = SpikeTrain(t * s, t_stop=t[-1] + 1.0)
        spike_trains.append(st)

    for sid in range(n_surrogates):
        surr_ms_list = []
        for i in range(n_units):
            try:
                surr = surrogates(
                    spike_trains[i], n_surrogates=1, method="shuffle_isis"
                )[0]
                surr_ms_list.append(surr.magnitude.flatten() * 1000.0)
            except Exception:
                surr_ms_list.append(real_times[i].copy())

        surr_peaks = ccg_peak_matrix(surr_ms_list)
        count_exceed += (surr_peaks >= real_peaks).astype(int)
        total_surr += 1

    p_values = np.ones((n_units, n_units))
    for i in range(n_units):
        for j in range(n_units):
            if i != j:
                p_values[i, j] = max(count_exceed[i, j], 1) / max(total_surr, 1)

    return p_values, n_units


def apply_bh_fdr(p_values: np.ndarray, fdr_q: float = 0.05) -> np.ndarray:
    """Apply Benjamini-Hochberg FDR correction. Returns boolean mask.

    Only off-diagonal pairs are tested (n_units * (n_units - 1) tests).
    """
    n_units = p_values.shape[0]
    n_pairs = n_units * (n_units - 1)

    mask = ~np.eye(n_units, dtype=bool)
    p_flat = p_values[mask]  # off-diagonal only
    n_tests = len(p_flat)

    sorted_idx = np.argsort(p_flat)
    reject_flat = np.zeros(n_tests, dtype=bool)
    for k, idx in enumerate(sorted_idx):
        threshold = fdr_q * (k + 1) / n_tests
        if p_flat[idx] <= threshold:
            reject_flat[idx] = True

    reject = np.zeros((n_units, n_units), dtype=bool)
    reject[mask] = reject_flat
    return reject


def validate_ccg(
    animal: str,
    condition: str,
    fdr_q: float = 0.05,
    n_surrogates: int = 100,
) -> dict:
    """Run full CCG-based validation pipeline. Returns stats dict."""
    print(f"Loading real spike trains for {animal}_{condition}...", flush=True)
    real_times = load_real_spike_times(animal, condition)
    n_units = len(real_times)
    n_pairs = n_units * (n_units - 1)
    print(f"  {n_units} neurons, {n_pairs} pairs", flush=True)

    print("Computing real CCG peaks...", flush=True)
    real_ms = real_times  # .mat values are in ms already
    real_peaks = ccg_peak_matrix(real_ms)
    n_real_pos = int((real_peaks > 0).sum())
    print(f"  Real: {n_real_pos} edges with peak >= 1", flush=True)

    print(f"Generating {n_surrogates} shuffle_isis surrogates...", flush=True)
    p_values, _ = compute_empirical_p_values(
        real_times, real_peaks, n_surrogates
    )

    reject = apply_bh_fdr(p_values, fdr_q)

    out_dir = PROJECT_ROOT / "results" / "glmcc"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"validated_adj_ccg_{animal}_{condition}.csv"
    validated = real_peaks.copy()
    validated[~reject] = 0.0
    np.savetxt(str(out_path), validated, delimiter=",", fmt="%d")

    # Also write p-values for transparency
    p_path = out_dir / f"p_values_ccg_{animal}_{condition}.csv"
    np.savetxt(str(p_path), p_values, delimiter=",", fmt="%.6f")

    n_sig = int(reject.sum())
    survival = n_sig / n_pairs if n_pairs > 0 else 0
    stats = {
        "animal": animal,
        "condition": condition,
        "method": "ccg_shuffle_isis",
        "n_units": n_units,
        "n_surrogates": n_surrogates,
        "n_possible_edges": n_pairs,
        "n_significant": n_sig,
        "edge_survival_rate": round(survival, 4),
        "fdr_q": fdr_q,
        "output_path": str(out_path),
    }
    return stats


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--animal")
    parser.add_argument("--condition", choices=["lightON", "ongoing"])
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--n-surrogates", type=int, default=100)
    parser.add_argument("--fdr-q", type=float, default=0.05)
    args = parser.parse_args()

    global N_SURROGATES
    N_SURROGATES = args.n_surrogates

    jobs = []
    if args.all:
        jobs = [(a, c) for a in ANIMALS for c in CONDITIONS]
    elif args.animal and args.condition:
        jobs = [(args.animal, args.condition)]
    else:
        parser.error("Specify --animal + --condition or --all")

    for a, c in jobs:
        out_path = PROJECT_ROOT / "results" / "glmcc" / f"validated_adj_ccg_{a}_{c}.csv"
        if out_path.exists():
            print(f"SKIP {a}_{c} — already exists", flush=True)
            continue
        stats = validate_ccg(a, c, args.fdr_q, args.n_surrogates)
        n_possible = stats['n_units'] * (stats['n_units'] - 1)
        print(f"  {a}_{c}: {stats['n_significant']}/{n_possible} edges significant, "
              f"survival={stats['edge_survival_rate']}", flush=True)


if __name__ == "__main__":
    main()
