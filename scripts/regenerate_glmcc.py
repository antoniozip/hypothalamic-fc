#!/usr/bin/env python3
"""Regenerate GLMCC connectivity with correct units and window.

Both GLMCC implementations expect 0-based MILLISECOND spike times over a recording of
length T, and use a +/-50 unit correlogram window (WIN=50, DELTA=1). This dataset stores
spike times in SECONDS, spanning 10,000-19,800 s and starting near 6,000 s, so:

  fed as-is (seconds) -> every spike passes the filter, but the window is +/-50 SECONDS
  fed as milliseconds -> the hardcoded T_SEC=5400 cutoff discards 100% of spikes

Neither is usable. This script converts to 0-based ms and passes the true recording
duration to vendor/glmcc-c/glmcc_fixed, giving a genuine +/-50 ms window.

Output: results/glmcc_fixed/adj_{animal}_{condition}.csv (+ .meta.json)
        results/glmcc_regeneration_summary.csv
"""

import argparse, csv, json, shutil, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "vendor" / "glmcc-c" / "glmcc_fixed"
OUT = ROOT / "results" / "glmcc_fixed"
ANIMAL_MAP = {"day1":"171019","day2":"180131","day3":"180302","day4":"180419","day5":"180420",
              "day6":"180423","night1":"171207","night2":"171208","night3":"171213",
              "night4":"180110","night5":"180111","night6":"180221","night7":"180228"}


def prepare(src: Path, dst: Path) -> tuple[int, float, int]:
    """Write 0-based millisecond cell files; return (n_cells, duration_s, n_spikes)."""
    dst.mkdir(parents=True, exist_ok=True)
    cells = sorted(src.glob("cell*.txt"), key=lambda p: int(p.stem[4:]))
    trains = [np.atleast_1d(np.loadtxt(c)) for c in cells]
    nonempty = [t for t in trains if t.size]
    if not nonempty:
        return 0, 0.0, 0
    t0 = min(t.min() for t in nonempty)
    t1 = max(t.max() for t in nonempty)
    total = 0
    for path, t in zip(cells, trains):
        shifted = (t - t0) * 1000.0 if t.size else t   # seconds -> ms, 0-based
        total += shifted.size
        np.savetxt(dst / path.name, np.sort(shifted), fmt="%.3f")
    return len(cells), float(t1 - t0), total


def run(animal: str, condition: str, workroot: Path) -> dict:
    old = ANIMAL_MAP[animal]
    src = ROOT / (f"DATA{old}" if condition == "lightON" else f"DATA{old}_{condition}")
    row = {"animal": animal, "condition": condition, "status": "", "note": ""}
    if not src.is_dir():
        row.update(status="missing", note=f"no {src.name}")
        return row

    work = workroot / f"{animal}_{condition}"
    if work.exists():
        shutil.rmtree(work)
    n_cells, dur_s, n_spikes = prepare(src, work)
    if n_cells == 0:
        row.update(status="empty", note="no spikes")
        return row
    row.update(n_units=n_cells, duration_s=round(dur_s, 1), n_spikes=n_spikes)

    started = time.time()
    # Pad T so nothing is truncated by the [0, T*1000) filter.
    proc = subprocess.run([str(BIN), ".", str(n_cells), "exp", "GLM", f"{dur_s * 1.01:.1f}"],
                          cwd=work, capture_output=True, text=True)
    row["elapsed_s"] = round(time.time() - started, 1)
    produced = sorted(work.glob("W_py_*.csv"))
    if proc.returncode != 0 or not produced:
        row.update(status="error", note=(proc.stderr or "no output").strip()[:150])
        return row

    W = np.loadtxt(produced[0], delimiter=",")
    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / f"adj_{animal}_{condition}.csv"
    np.savetxt(dest, W, delimiter=",", fmt="%.6f")

    n = W.shape[0]
    off = ~np.eye(n, dtype=bool)
    nz = W[off][W[off] != 0]
    row.update(status="ok", nonzero=len(nz),
               density_pct=round(100 * len(nz) / off.sum(), 2),
               neg_pct=round(100 * float((nz < 0).mean()), 1) if len(nz) else "",
               median_w=round(float(np.median(nz)), 4) if len(nz) else "")
    dest.with_suffix(".meta.json").write_text(json.dumps({
        "encoding": "glmcc_weight_raw", "units": "0-based milliseconds",
        "window_ms": 50.0, "bin_ms": 1.0, "T_seconds": round(dur_s * 1.01, 1),
        "source_dir": src.name, "n_units": n_cells, "n_spikes": n_spikes,
        "binary": "vendor/glmcc-c/glmcc_fixed",
        "written_at": datetime.now(timezone.utc).isoformat()}, indent=2))
    return row


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--animals", nargs="+", default=list(ANIMAL_MAP))
    ap.add_argument("--conditions", nargs="+", default=["lightON", "ongoing"])
    ap.add_argument("--workroot", type=Path, default=Path("/tmp/claude-1000/glmcc_fixed"))
    args = ap.parse_args()
    args.workroot.mkdir(parents=True, exist_ok=True)

    rows = []
    for cond in args.conditions:
        for a in args.animals:
            r = run(a, cond, args.workroot)
            rows.append(r)
            print(f"[{a}_{cond}] {r['status']} density={r.get('density_pct','-')}% "
                  f"neg={r.get('neg_pct','-')}% n={r.get('n_units','-')} "
                  f"in {r.get('elapsed_s','-')}s {r.get('note','')}", flush=True)

    fields = ["animal","condition","status","n_units","duration_s","n_spikes","nonzero",
              "density_pct","neg_pct","median_w","elapsed_s","note"]
    out = ROOT / "results" / "glmcc_regeneration_summary.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
        for r in rows: w.writerow({k: r.get(k, "") for k in fields})
    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"\n{ok}/{len(rows)} completed. Summary: {out}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
