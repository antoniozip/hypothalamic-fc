#!/usr/bin/env python3
"""Re-run GLMCC on `ongoing` data gapped to match each animal's own lightON duty cycle.

The lightON epoch is a ~10% duty-cycle concatenation of stimulus windows, while ongoing is
near-continuous. GLMCC fits its baseline assuming continuous recording, so the two conditions
are not comparable (confirmed on night3: gapping ongoing flips it from 98.6% negative weights
to 85.2% positive). Imposing the same gap structure on ongoing puts both conditions under the
same bias, making the condition contrast interpretable.

Output: results/glmcc/adj_{animal}_ongoing_gapped.csv (+ .meta.json provenance)
        results/gapped_ongoing_summary.csv

Usage: python scripts/run_gapped_ongoing.py [--animals day1 night3] [--jobs 1]
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
from gap_spike_trains import gap_directory, measure_windows  # noqa: E402

GLMCC_BIN = PROJECT_ROOT / "vendor" / "glmcc-c" / "glmcc"
RESULTS = PROJECT_ROOT / "results" / "glmcc"

ANIMAL_MAP = {
    "day1": "171019", "day2": "180131", "day3": "180302",
    "day4": "180419", "day5": "180420", "day6": "180423",
    "night1": "171207", "night2": "171208", "night3": "171213",
    "night4": "180110", "night5": "180111", "night6": "180221", "night7": "180228",
}


def run_one(animal: str, workroot: Path) -> dict:
    old = ANIMAL_MAP[animal]
    light_dir = PROJECT_ROOT / f"DATA{old}"
    ongoing_dir = PROJECT_ROOT / f"DATA{old}_ongoing"
    row = {"animal": animal, "old_id": old, "status": "", "note": ""}

    if not light_dir.is_dir() or not ongoing_dir.is_dir():
        row["status"] = "missing"
        row["note"] = f"missing {light_dir.name} or {ongoing_dir.name}"
        return row

    # Each animal's own stimulus structure, not a global constant.
    on_s, period_s, duty_light = measure_windows(light_dir)
    row.update({"lightON_window_s": round(on_s, 2),
                "lightON_period_s": round(period_s, 1),
                "lightON_duty_pct": round(duty_light, 1)})

    work = workroot / animal
    if work.exists():
        shutil.rmtree(work)
    stats = gap_directory(ongoing_dir, work, on_s, period_s)
    row.update({
        "n_windows": stats["n_windows"],
        "spikes_kept": stats["spikes_kept"],
        "spikes_total": stats["spikes_total"],
        "gapped_duty_pct": round(stats["duty_cycle_pct"], 1),
        "gapped_rate_hz": round(stats["rate_hz_per_neuron"], 3),
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
    out = RESULTS / f"adj_{animal}_ongoing_gapped.csv"
    np.savetxt(out, W, delimiter=",", fmt="%.6f")

    n = W.shape[0]
    off = ~np.eye(n, dtype=bool)
    nz = W[off][W[off] != 0]
    row.update({
        "status": "ok", "n_units": n, "nonzero": len(nz),
        "neg_pct": round(100 * (nz < 0).mean(), 1) if len(nz) else "",
        "median_w": round(float(np.median(nz)), 4) if len(nz) else "",
    })

    out.with_suffix(".meta.json").write_text(json.dumps({
        "encoding": "glmcc_weight_raw",
        "condition": "ongoing_gapped",
        "gap_source": f"DATA{old} lightON window structure",
        "window_s": on_s, "period_s": period_s,
        "duty_cycle_pct": stats["duty_cycle_pct"],
        "spikes_kept": stats["spikes_kept"], "spikes_total": stats["spikes_total"],
        "glmcc": "vendor/glmcc-c/glmcc exp GLM",
        "written_at": datetime.now(timezone.utc).isoformat(),
    }, indent=2))
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--animals", nargs="+", default=list(ANIMAL_MAP))
    parser.add_argument("--workroot", type=Path,
                        default=Path("/tmp/claude-1000/gapped_glmcc"))
    args = parser.parse_args()

    args.workroot.mkdir(parents=True, exist_ok=True)
    rows = []
    for animal in args.animals:
        print(f"\n{'=' * 66}\n{animal}\n{'=' * 66}", flush=True)
        row = run_one(animal, args.workroot)
        rows.append(row)
        print(f"[{animal}] {row['status']} "
              f"duty={row.get('gapped_duty_pct','-')}% neg={row.get('neg_pct','-')}% "
              f"in {row.get('elapsed_s','-')}s", flush=True)

    fields = ["animal", "old_id", "status", "n_units", "lightON_window_s", "lightON_period_s",
              "lightON_duty_pct", "n_windows", "gapped_duty_pct", "spikes_kept", "spikes_total",
              "gapped_rate_hz", "nonzero", "neg_pct", "median_w", "elapsed_s", "note"]
    out = PROJECT_ROOT / "results" / "gapped_ongoing_summary.csv"
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
