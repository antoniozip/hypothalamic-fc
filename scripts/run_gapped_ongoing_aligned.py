#!/usr/bin/env python3
"""Re-run GLMCC on `ongoing` data gapped with the REAL stimulus window geometry.

Improves on scripts/run_gapped_ongoing.py, which imposed a regular period fitted to each
animal's median window. Here the true pulse edges are read from `mask_stimulus` in the .mat
file (exactly 10.00 s pulses every 90.0 s) and shifted by --offset seconds so the windows fall
in the light-off intervals. The result is a control with identical observation geometry to
lightON -- same window count, duration and period -- but sampling genuinely light-off data.

`ongoing` is the complement of `lightON` (verified: <20 of ~10^6 ongoing spikes fall inside a
pulse), so a phase-shifted mask selects real ongoing data with no overlap.

Output: results/glmcc/adj_{animal}_ongoing_aligned.csv (+ .meta.json)
        results/aligned_ongoing_summary.csv
"""

import argparse
import csv
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from gap_spike_trains import gap_directory_windows, stimulus_windows  # noqa: E402

GLMCC_BIN = PROJECT_ROOT / "vendor" / "glmcc-c" / "glmcc"
RESULTS = PROJECT_ROOT / "results" / "glmcc"

ANIMAL_MAP = {
    "day1": "171019", "day2": "180131", "day3": "180302",
    "day4": "180419", "day5": "180420", "day6": "180423",
    "night1": "171207", "night2": "171208", "night3": "171213",
    "night4": "180110", "night5": "180111", "night6": "180221", "night7": "180228",
}


def run_one(animal: str, workroot: Path, offset_s: float) -> dict:
    old = ANIMAL_MAP[animal]
    ongoing_dir = PROJECT_ROOT / f"DATA{old}_ongoing"
    row = {"animal": animal, "old_id": old, "offset_s": offset_s, "status": "", "note": ""}

    mat = None
    for cand in (PROJECT_ROOT / f"DATA{old}_ongoing.mat", PROJECT_ROOT / f"DATA{old}_lightON.mat"):
        if cand.exists():
            mat = cand
            break
    if mat is None or not ongoing_dir.is_dir():
        row["status"] = "missing"
        row["note"] = f"no .mat with mask_stimulus, or missing {ongoing_dir.name}"
        return row

    try:
        starts, ends = stimulus_windows(mat)
    except KeyError as exc:
        row["status"] = "missing"
        row["note"] = str(exc)
        return row

    durations = ends - starts
    row.update({
        "n_pulses": len(starts),
        "pulse_s": round(float(np.median(durations)), 2),
        "period_s": round(float(np.median(np.diff(starts))), 1) if len(starts) > 1 else "",
        "mask_source": mat.name,
    })

    # Shift into the light-off interval, keeping window length identical.
    shifted_starts = starts + offset_s
    shifted_ends = shifted_starts + durations
    if len(starts) > 1:
        period = float(np.median(np.diff(starts)))
        if offset_s + float(np.median(durations)) > period:
            row["status"] = "error"
            row["note"] = f"offset {offset_s}s + pulse overruns period {period}s"
            return row

    work = workroot / animal
    if work.exists():
        shutil.rmtree(work)
    stats = gap_directory_windows(ongoing_dir, work, shifted_starts, shifted_ends)
    row.update({
        "spikes_kept": stats["spikes_kept"], "spikes_total": stats["spikes_total"],
        "duty_pct": round(stats["duty_cycle_pct"], 1),
        "observed_s": round(stats["observed_s"], 0),
        "rate_hz": round(stats["rate_hz_per_neuron"], 3),
    })

    n_cells = len(list(ongoing_dir.glob("cell*.txt")))
    started = time.time()
    proc = subprocess.run([str(GLMCC_BIN), ".", str(n_cells), "exp", "GLM"],
                          cwd=work, capture_output=True, text=True)
    row["elapsed_s"] = round(time.time() - started, 1)

    produced = sorted(work.glob("W_py_*.csv"))
    if proc.returncode != 0 or not produced:
        row["status"] = "error"
        row["note"] = (proc.stderr or "no W_py output").strip()[:200]
        return row

    W = np.loadtxt(produced[0], delimiter=",")
    out = RESULTS / f"adj_{animal}_ongoing_aligned.csv"
    np.savetxt(out, W, delimiter=",", fmt="%.6f")

    n = W.shape[0]
    off = ~np.eye(n, dtype=bool)
    nz = W[off][W[off] != 0]
    row.update({
        "status": "ok", "n_units": n, "nonzero": len(nz),
        "neg_pct": round(100 * float((nz < 0).mean()), 1) if len(nz) else "",
        "median_w": round(float(np.median(nz)), 4) if len(nz) else "",
    })

    out.with_suffix(".meta.json").write_text(json.dumps({
        "encoding": "glmcc_weight_raw",
        "condition": "ongoing_aligned",
        "mask_source": mat.name,
        "pulse_s": float(np.median(durations)),
        "offset_s": offset_s,
        "n_windows": len(starts),
        "duty_cycle_pct": stats["duty_cycle_pct"],
        "spikes_kept": stats["spikes_kept"], "spikes_total": stats["spikes_total"],
        "glmcc": "vendor/glmcc-c/glmcc exp GLM",
        "written_at": datetime.now(timezone.utc).isoformat(),
    }, indent=2))
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--animals", nargs="+", default=list(ANIMAL_MAP))
    parser.add_argument("--offset", type=float, default=45.0,
                        help="seconds to shift each pulse window into the light-off interval")
    parser.add_argument("--workroot", type=Path, default=Path("/tmp/claude-1000/aligned_glmcc"))
    args = parser.parse_args()

    args.workroot.mkdir(parents=True, exist_ok=True)
    rows = []
    for animal in args.animals:
        print(f"\n{'=' * 66}\n{animal}\n{'=' * 66}", flush=True)
        row = run_one(animal, args.workroot, args.offset)
        rows.append(row)
        print(f"[{animal}] {row['status']} pulses={row.get('n_pulses','-')} "
              f"duty={row.get('duty_pct','-')}% neg={row.get('neg_pct','-')}% "
              f"in {row.get('elapsed_s','-')}s {row.get('note','')}", flush=True)

    fields = ["animal", "old_id", "status", "n_units", "mask_source", "n_pulses", "pulse_s",
              "period_s", "offset_s", "duty_pct", "observed_s", "spikes_kept", "spikes_total",
              "rate_hz", "nonzero", "neg_pct", "median_w", "elapsed_s", "note"]
    out = PROJECT_ROOT / "results" / "aligned_ongoing_summary.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})
    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"\n{ok}/{len(rows)} animals completed. Summary: {out}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
