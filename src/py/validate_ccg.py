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
    python validate_ccg.py --animal day1 --condition lightON
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

ANIMALS = ["day1","night1","night2","night3","night4","night5","day2","night6",
           "night7","day3","day4","day5","day6"]
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
        # Spike-free units are kept as empty arrays. Dropping them shifts every
        # later matrix index relative to neuron_id in data/processed/neurons.csv,
        # which src/r/compute_graph_metrics.R maps positionally.
        times.append(t.astype(float))  # .mat values are in ms
    return times


def ccg_peak_one_pair(a_ms: np.ndarray, b_ms: np.ndarray) -> float:
    """CCG peak for a single neuron pair: the largest count in any BIN_MS bin
    of the cross-correlogram over lags in [-WINDOW_MS, +WINDOW_MS].

    Only spike pairs that actually fall inside the window are materialised, via
    searchsorted, so the cost tracks the number of coincidences rather than
    len(a) * len(b). Results are identical to the dense-difference form.
    """
    if len(a_ms) < 2 or len(b_ms) < 2:
        return 0.0

    n_bins = int(2 * WINDOW_MS / BIN_MS) + 1
    hist = np.zeros(n_bins, dtype=np.int64)

    # Chunk over b so the coincidence list stays bounded on dense recordings.
    chunk = 200_000
    for i0 in range(0, len(b_ms), chunk):
        sub = b_ms[i0:i0 + chunk]
        lo = np.searchsorted(a_ms, sub - WINDOW_MS, side="left")
        hi = np.searchsorted(a_ms, sub + WINDOW_MS, side="right")
        counts = hi - lo
        total = int(counts.sum())
        if total == 0:
            continue

        # Expand (b spike, matching a spikes) into flat index arrays.
        ends = np.cumsum(counts)
        starts = ends - counts
        offsets = np.arange(total) - np.repeat(starts, counts)
        a_idx = np.repeat(lo, counts) + offsets
        b_idx = np.repeat(np.arange(len(sub)), counts)

        diffs = sub[b_idx] - a_ms[a_idx]
        diffs = diffs[np.abs(diffs) <= WINDOW_MS]
        if diffs.size:
            idx = ((diffs + WINDOW_MS) / BIN_MS).astype(np.int64)
            hist += np.bincount(idx, minlength=n_bins)[:n_bins]

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
    seed: int = 42,
) -> tuple[np.ndarray, int]:
    """Compute empirical p-values via shuffle-ISI surrogates.

    Returns (p_values, n_units) where p_values[i,i] = 1.0 (diagonal ignored).
    """
    n_units = len(real_times)
    count_exceed = np.zeros((n_units, n_units), dtype=int)
    total_surr = 0

    rng = np.random.default_rng(seed)
    np.random.seed(seed)  # elephant's surrogates() draws from the legacy global RNG

    # Build Neo SpikeTrain objects once. A unit with fewer than two spikes has no
    # ISI distribution to shuffle, so it is carried through unchanged.
    t_stop_global = max((t[-1] for t in real_times if len(t) > 0), default=1.0) + 1.0
    spike_trains = []
    for t in real_times:
        if len(t) == 0:
            spike_trains.append(None)
        else:
            spike_trains.append(SpikeTrain(t * s, t_stop=t_stop_global))

    for sid in range(n_surrogates):
        surr_ms_list = []
        for i in range(n_units):
            if spike_trains[i] is None or len(real_times[i]) < 3:
                surr_ms_list.append(real_times[i].copy())
                continue
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

    mask = ~np.eye(n_units, dtype=bool)
    p_flat = p_values[mask]  # off-diagonal only
    n_tests = len(p_flat)

    # Step-up procedure: find the largest rank k with p_(k) <= q*k/m, then reject
    # every hypothesis up to that rank. Testing each p_(k) against its own
    # threshold instead is not BH -- with p-values on a 1/n_surrogates grid the
    # ties are broken by argsort order, so edge selection follows neuron
    # numbering rather than the data.
    sorted_idx = np.argsort(p_flat, kind="stable")
    p_sorted = p_flat[sorted_idx]
    thresholds = fdr_q * np.arange(1, n_tests + 1) / n_tests
    passing = np.flatnonzero(p_sorted <= thresholds)

    reject_flat = np.zeros(n_tests, dtype=bool)
    if passing.size > 0:
        reject_flat[sorted_idx[: passing.max() + 1]] = True

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
