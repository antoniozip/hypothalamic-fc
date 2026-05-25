#!/usr/bin/env python3
"""Sprint 5: Temporal structure of the LightON paradigm.

5.1 Adaptation check: compare first vs last hour of the stimulus session.
5.2 Early vs late evoked window: within each 10 s ON pulse, split into
    early [0, 500 ms] and late [500 ms, 10 s].

Output: results/sprint5/ directory with spiketrains, connectivity, and stats.
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from validate_ccg import ccg_peak_matrix

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
OUT_DIR = PROJECT_ROOT / "results" / "sprint5"

ANIMALS = ["day1","night1","night2","night3","night4","night5","day2","night6",
           "night7","day3","day4","day5","day6"]

TS_KEY_MAP = {"lightON": "evoked_ts", "ongoing": "ongoing_ts"}


def load_data(animal: str, condition: str = "lightON"):
    mat_path = PROJECT_ROOT / f"DATA{animal}_{condition}.mat"
    if not mat_path.exists():
        alts = sorted(PROJECT_ROOT.glob(f"DATA{animal}*{condition}*.mat"))
        mat_path = alts[0]
    m = loadmat(str(mat_path), simplify_cells=False)
    ts_key = TS_KEY_MAP.get(condition, "evoked_ts")
    if ts_key not in m:
        ts_key = [k for k in m if k.endswith("_ts")][0]
    ts = m[ts_key]
    units = m["Units"]
    mask = m["mask_stimulus"].flatten() if "mask_stimulus" in m else None
    max_t = float(m["max_t"][0, 0])
    return ts, units, mask, max_t


def get_pulse_edges(mask):
    """Return (pulse_starts, pulse_ends) in seconds."""
    on = mask > 0
    trans = np.diff(on.astype(int))
    starts = np.where(trans == 1)[0] / 1000.0
    ends = np.where(trans == -1)[0] / 1000.0
    if len(ends) > len(starts):
        ends = ends[:len(starts)]
    if len(starts) > len(ends):
        starts = starts[:len(ends)]
    return starts, ends


def extract_spikes(ts, t_start, t_end, n_units):
    """Extract spikes in [t_start, t_end) for all neurons. Times in seconds."""
    unit_spikes = []
    for i in range(n_units):
        t = np.atleast_1d(ts[0, i][0]) / 1000.0  # ms to s
        mask = (t >= t_start) & (t < t_end)
        unit_spikes.append(t[mask])
    return unit_spikes


def compute_ccg_and_rates(spike_list):
    """Compute firing rates (Hz) and CCG peak matrix."""
    rates = np.array([len(s) / (s[-1] - s[0]) if len(s) > 1 else 0
                      for s in spike_list])
    ms_list = [s * 1000.0 for s in spike_list if len(s) > 1]
    idx_map = [i for i, s in enumerate(spike_list) if len(s) > 1]
    if len(ms_list) < 2:
        return rates, np.zeros((len(spike_list), len(spike_list)))
    ccg = ccg_peak_matrix(ms_list)
    full_ccg = np.zeros((len(spike_list), len(spike_list)))
    for ni, i in enumerate(idx_map):
        for nj, j in enumerate(idx_map):
            full_ccg[i, j] = ccg[ni, nj]
    return rates, full_ccg


def run_adaptation(animal: str):
    """5.1: Compare first vs last hour of the stimulus period."""
    print(f"  {animal}: adaptation check...", flush=True)
    ts, units, mask, max_t = load_data(animal)
    n_units = ts.shape[1]

    if mask is None:
        print(f"  {animal}: no mask, skipping")
        return

    starts_s, ends_s = get_pulse_edges(mask)
    if len(starts_s) == 0:
        print(f"  {animal}: no pulses, skipping")
        return

    stimulus_start = starts_s[0]
    stimulus_end = ends_s[-1]
    stimulus_dur = stimulus_end - stimulus_start
    print(f"    Stimulus span: {stimulus_start:.0f}-{stimulus_end:.0f}s ({stimulus_dur/60:.0f} min)", flush=True)

    # Early epoch: first 25% of stimulus duration
    early_end = stimulus_start + stimulus_dur * 0.25
    # Late epoch: last 25%
    late_start = stimulus_start + stimulus_dur * 0.75

    early_spikes = extract_spikes(ts, stimulus_start, early_end, n_units)
    late_spikes = extract_spikes(ts, late_start, stimulus_end, n_units)

    early_rates, early_ccg = compute_ccg_and_rates(early_spikes)
    late_rates, late_ccg = compute_ccg_and_rates(late_spikes)

    # Firing rate comparison (paired t-test)
    valid = (early_rates > 0) & (late_rates > 0)
    if valid.sum() > 2:
        from scipy.stats import ttest_rel
        _, p_val = ttest_rel(early_rates[valid], late_rates[valid])
    else:
        p_val = 1.0

    early_strength = early_ccg.sum(axis=1)
    late_strength = late_ccg.sum(axis=1)

    out = {
        "animal": animal,
        "n_units": n_units,
        "stimulus_start_s": round(stimulus_start, 0),
        "stimulus_end_s": round(stimulus_end, 0),
        "early_dur_s": round(early_end - stimulus_start, 0),
        "late_dur_s": round(stimulus_end - late_start, 0),
        "early_mean_rate": round(early_rates.mean(), 3),
        "late_mean_rate": round(late_rates.mean(), 3),
        "rate_p_value": round(p_val, 4),
        "early_mean_strength": round(early_strength.mean(), 2),
        "late_mean_strength": round(late_strength.mean(), 2),
    }

    if p_val < 0.05:
        print(f"    ⚠ Rate changed: early={early_rates.mean():.2f} late={late_rates.mean():.2f} Hz (p={p_val:.4f})", flush=True)
    else:
        print(f"    ✓ Rate stable: early={early_rates.mean():.2f} late={late_rates.mean():.2f} Hz (p={p_val:.4f})", flush=True)

    return out


def run_early_late(animal: str):
    """5.2: Within each pulse, partition into early [0, 500ms] and late [500ms, 10s]."""
    print(f"  {animal}: early/late evoked window...", flush=True)
    ts, units, mask, max_t = load_data(animal)
    n_units = ts.shape[1]

    if mask is None:
        print(f"  {animal}: no mask, skipping")
        return

    starts_s, ends_s = get_pulse_edges(mask)
    if len(starts_s) == 0:
        return

    early_all = [[] for _ in range(n_units)]
    late_all = [[] for _ in range(n_units)]

    for pulse_idx in range(len(starts_s)):
        pulse_start = starts_s[pulse_idx]
        pulse_end = ends_s[pulse_idx]
        early_end = pulse_start + 0.5  # 500 ms
        late_start = pulse_end - 9.5   # 10s - 500ms onward
        # Actually, [500ms, 10s] means from 0.5s to 10s
        middle = pulse_start + 9.5
        early_end = pulse_start + 0.5

        for i in range(n_units):
            t = np.atleast_1d(ts[0, i][0]) / 1000.0
            early_mask = (t >= pulse_start) & (t < early_end)
            late_mask = (t >= middle) & (t < pulse_end)
            early_all[i].extend(t[early_mask].tolist())
            late_all[i].extend(t[late_mask].tolist())

    early_spikes = [np.array(x) for x in early_all]
    late_spikes = [np.array(x) for x in late_all]

    # Keep only cells with spikes
    early_ok = [i for i in range(n_units) if len(early_spikes[i]) > 1]
    late_ok = [i for i in range(n_units) if len(late_spikes[i]) > 1]

    early_rates, early_ccg = compute_ccg_and_rates(early_spikes)
    late_rates, late_ccg = compute_ccg_and_rates(late_spikes)

    # Rate comparison on cells with spikes in both
    both = np.array([i for i in range(n_units)
                     if len(early_spikes[i]) > 0 and len(late_spikes[i]) > 0])
    if len(both) > 2:
        from scipy.stats import ttest_rel
        _, p_val = ttest_rel(early_rates[both], late_rates[both])
    else:
        p_val = 1.0

    print(f"    Early: {len(early_ok)} cells, Late: {len(late_ok)} cells, "
          f"mean rates {early_rates[both].mean():.3f}/{late_rates[both].mean():.3f} Hz (p={p_val:.4f})",
          flush=True)

    return {
        "animal": animal,
        "n_pulses": len(starts_s),
        "n_units": n_units,
        "early_n_active": len(early_ok),
        "late_n_active": len(late_ok),
        "early_mean_rate": round(early_rates[both].mean(), 4) if len(both) > 0 else 0,
        "late_mean_rate": round(late_rates[both].mean(), 4) if len(both) > 0 else 0,
        "rate_p_value": round(p_val, 4),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--animal")
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()

    animals = ANIMALS if args.all else [args.animal]
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    results_51 = []
    results_52 = []

    for a in animals:
        try:
            r51 = run_adaptation(a)
            if r51:
                results_51.append(r51)
        except Exception as e:
            print(f"  {a}: adaptation ERROR: {e}", flush=True)

        try:
            r52 = run_early_late(a)
            if r52:
                results_52.append(r52)
        except Exception as e:
            print(f"  {a}: early/late ERROR: {e}", flush=True)

    # Write CSV
    if results_51:
        import csv
        with open(OUT_DIR / "adaptation_check.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=results_51[0].keys())
            w.writeheader()
            w.writerows(results_51)
        print(f"\nAdaptation check written to {OUT_DIR / 'adaptation_check.csv'}")

    if results_52:
        import csv
        with open(OUT_DIR / "early_late_evoked.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=results_52[0].keys())
            w.writeheader()
            w.writerows(results_52)
        print(f"Early/late evoked written to {OUT_DIR / 'early_late_evoked.csv'}")


if __name__ == "__main__":
    main()
