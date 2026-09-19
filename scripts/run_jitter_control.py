#!/usr/bin/env python3
"""Jitter control: is the lightON signal fine-timing coupling or slow common drive?

The duty-cycle work established that the lightON/ongoing sign asymmetry is mostly an artifact
of GLMCC's baseline, and that a duty-cycle-matched light-off control lands at ~28% negative
against lightON's ~0.3%. That residual could be genuine synaptic coupling, or it could be
stimulus-locked rate co-modulation -- many neurons driven together by a 10 s light pulse produce
CCG peaks with no synapse between them.

Jittering each spike by +/-25 ms destroys 1-5 ms synaptic structure while leaving the 10 s
stimulus envelope intact. Read the outcome as:

  lightON_jittered moves toward the light-off control  -> the residual was fine-timing coupling
  lightON_jittered stays near lightON                  -> the residual is slow common drive

Both lightON and the aligned light-off control are jittered, so the comparison is symmetric.

Output: results/glmcc/adj_{animal}_lightON_jittered.csv
        results/glmcc/adj_{animal}_ongoing_aligned_jittered.csv  (+ .meta.json each)
        results/jitter_control_summary.csv
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
from jitter_spike_trains import jitter_directory  # noqa: E402

GLMCC_BIN = PROJECT_ROOT / "vendor" / "glmcc-c" / "glmcc"
RESULTS = PROJECT_ROOT / "results" / "glmcc"

ANIMAL_MAP = {
    "day1": "171019", "day2": "180131", "day3": "180302",
    "day4": "180419", "day5": "180420", "day6": "180423",
    "night1": "171207", "night2": "171208", "night3": "171213",
    "night4": "180110", "night5": "180111", "night6": "180221", "night7": "180228",
}


def _glmcc(work: Path, n_cells: int) -> tuple[np.ndarray | None, float, str]:
    started = time.time()
    proc = subprocess.run([str(GLMCC_BIN), ".", str(n_cells), "exp", "GLM"],
                          cwd=work, capture_output=True, text=True)
    elapsed = round(time.time() - started, 1)
    produced = sorted(work.glob("W_py_*.csv"))
    if proc.returncode != 0 or not produced:
        return None, elapsed, (proc.stderr or "no W_py output").strip()[:200]
    return np.loadtxt(produced[0], delimiter=","), elapsed, ""


def _sign_stats(W: np.ndarray) -> dict:
    n = W.shape[0]
    off = ~np.eye(n, dtype=bool)
    nz = W[off][W[off] != 0]
    if not len(nz):
        return {"n_units": n, "nonzero": 0, "neg_pct": "", "median_w": ""}
    return {"n_units": n, "nonzero": len(nz),
            "neg_pct": round(100 * float((nz < 0).mean()), 1),
            "median_w": round(float(np.median(nz)), 4)}


def run_one(animal: str, workroot: Path, jitter_ms: float, offset_s: float, seed: int) -> list[dict]:
    old = ANIMAL_MAP[animal]
    light_dir = PROJECT_ROOT / f"DATA{old}"
    ongoing_dir = PROJECT_ROOT / f"DATA{old}_ongoing"
    rows = []

    mat = next((c for c in (PROJECT_ROOT / f"DATA{old}_ongoing.mat",
                            PROJECT_ROOT / f"DATA{old}_lightON.mat") if c.exists()), None)
    if mat is None or not light_dir.is_dir() or not ongoing_dir.is_dir():
        return [{"animal": animal, "variant": v, "status": "missing",
                 "note": "missing .mat or DATA dir"}
                for v in ("lightON_jittered", "ongoing_aligned_jittered")]

    starts, ends = stimulus_windows(mat)
    durations = ends - starts

    jobs = [
        # lightON: jitter in place, inside the true pulse windows.
        ("lightON_jittered", light_dir, starts, ends, None),
        # light-off control: gap ongoing to the shifted windows, then jitter identically.
        ("ongoing_aligned_jittered", ongoing_dir,
         starts + offset_s, starts + offset_s + durations, "gap_first"),
    ]

    for variant, src, w_start, w_end, mode in jobs:
        row = {"animal": animal, "old_id": old, "variant": variant,
               "jitter_ms": jitter_ms, "seed": seed, "status": "", "note": ""}
        work = workroot / f"{animal}_{variant}"
        if work.exists():
            shutil.rmtree(work)

        if mode == "gap_first":
            staging = workroot / f"{animal}_{variant}_staged"
            if staging.exists():
                shutil.rmtree(staging)
            gap_directory_windows(src, staging, w_start, w_end)
            src_for_jitter = staging
        else:
            src_for_jitter = src

        stats = jitter_directory(src_for_jitter, work, w_start, w_end, jitter_ms, seed)
        row.update({"spikes_kept": stats["spikes_kept"], "spikes_total": stats["spikes_total"],
                    "dropped_at_edges": stats["dropped_at_edges"],
                    "retention_pct": round(stats["retention_pct"], 2)})

        n_cells = len(list(src.glob("cell*.txt")))
        W, elapsed, err = _glmcc(work, n_cells)
        row["elapsed_s"] = elapsed
        if W is None:
            row.update({"status": "error", "note": err})
            rows.append(row)
            continue

        out = RESULTS / f"adj_{animal}_{variant}.csv"
        np.savetxt(out, W, delimiter=",", fmt="%.6f")
        row.update({"status": "ok", **_sign_stats(W)})
        out.with_suffix(".meta.json").write_text(json.dumps({
            "encoding": "glmcc_weight_raw", "condition": variant,
            "jitter_ms": jitter_ms, "seed": seed, "offset_s": offset_s if mode else 0,
            "mask_source": mat.name, "n_windows": len(w_start),
            "spikes_kept": stats["spikes_kept"], "spikes_total": stats["spikes_total"],
            "glmcc": "vendor/glmcc-c/glmcc exp GLM",
            "written_at": datetime.now(timezone.utc).isoformat(),
        }, indent=2))
        rows.append(row)

    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--animals", nargs="+", default=list(ANIMAL_MAP))
    parser.add_argument("--jitter-ms", type=float, default=25.0)
    parser.add_argument("--offset", type=float, default=45.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--workroot", type=Path, default=Path("/tmp/claude-1000/jitter_glmcc"))
    args = parser.parse_args()

    args.workroot.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for animal in args.animals:
        print(f"\n{'=' * 66}\n{animal}\n{'=' * 66}", flush=True)
        for row in run_one(animal, args.workroot, args.jitter_ms, args.offset, args.seed):
            all_rows.append(row)
            print(f"[{animal}/{row['variant']}] {row['status']} "
                  f"neg={row.get('neg_pct','-')}% retention={row.get('retention_pct','-')}% "
                  f"in {row.get('elapsed_s','-')}s {row.get('note','')}", flush=True)

    fields = ["animal", "old_id", "variant", "status", "n_units", "jitter_ms", "seed",
              "spikes_kept", "spikes_total", "dropped_at_edges", "retention_pct",
              "nonzero", "neg_pct", "median_w", "elapsed_s", "note"]
    out = PROJECT_ROOT / "results" / "jitter_control_summary.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in all_rows:
            w.writerow({k: r.get(k, "") for k in fields})
    ok = sum(1 for r in all_rows if r["status"] == "ok")
    print(f"\n{ok}/{len(all_rows)} runs completed. Summary: {out}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
