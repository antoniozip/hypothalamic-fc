#!/usr/bin/env python3
"""Jitter sensitivity sweep: at what timescale does the GLMCC edge set break?

`scripts/run_jitter_control.py` displaces every spike by +/-25 ms and found that 41,827 of
41,842 aligned-ongoing edges survive it (100.0%). A synaptic edge lives on 1-5 ms structure and
cannot survive a 25 ms smear, so either the edges are not synaptic or the control is not doing
what it claims. This sweep separates those by varying the displacement and watching where the
edge set starts to move.

Design
------
One arm only: aligned-ongoing (ongoing gapped to the real pulse geometry, phase-shifted into the
light-off intervals). Its geometry is fixed across every level, so the only thing that changes
between runs is how far spikes move. The unjittered run of that same arm is the baseline every
level is compared against.

Levels, in milliseconds, plus two reference points that make the curve readable:

  0        identity check -- no displacement. Must reproduce the baseline exactly. If it does
           not, the comparison machinery is broken and no other number here means anything.
  5..100   the requested sweep. 5 ms is at the upper edge of a monosynaptic lag; 100 ms is far
           beyond any synaptic timescale but still well inside a 10 s window.
  shuffle  floor. Every spike redrawn uniformly inside the window it came from, preserving
           per-cell per-window counts. All sub-window temporal structure is gone, so whatever
           survives this is carried by spike count and window geometry alone.

Spikes displaced outside their originating window are dropped, not clipped (clipping piles
spikes at edges and manufactures synchrony). That loss grows with the displacement, so
retention is recorded per level: an edge count that falls because spikes were lost is a
different finding from one that falls because coupling was destroyed.

Output: results/jitter_sensitivity.csv
        results/glmcc/sensitivity/adj_{animal}_jit{level}.csv
"""

import argparse
import csv
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from gap_spike_trains import gap_directory_windows, stimulus_windows  # noqa: E402
from jitter_spike_trains import jitter_directory  # noqa: E402

GLMCC_BIN = PROJECT_ROOT / "vendor" / "glmcc-c" / "glmcc"
OUT_DIR = PROJECT_ROOT / "results" / "glmcc" / "sensitivity"
SUMMARY = PROJECT_ROOT / "results" / "jitter_sensitivity.csv"

ANIMAL_MAP = {
    "day1": "171019", "day2": "180131", "day3": "180302",
    "day4": "180419", "day5": "180420", "day6": "180423",
    "night1": "171207", "night2": "171208", "night3": "171213",
    "night4": "180110", "night5": "180111", "night6": "180221", "night7": "180228",
}

FIELDS = [
    "animal", "level", "jitter_ms", "status", "n_units",
    "spikes_kept", "spikes_total", "retention_pct",
    "n_base", "n_jit", "overlap", "survival_pct", "jaccard",
    "sign_agree_pct", "weight_r", "neg_pct", "elapsed_s", "note",
]


def _glmcc(work: Path, n_cells: int):
    started = time.time()
    proc = subprocess.run([str(GLMCC_BIN), ".", str(n_cells), "exp", "GLM"],
                          cwd=work, capture_output=True, text=True)
    elapsed = round(time.time() - started, 1)
    produced = sorted(work.glob("W_py_*.csv"))
    if proc.returncode != 0 or not produced:
        return None, elapsed, (proc.stderr or "no W_py output").strip()[:180]
    return np.loadtxt(produced[0], delimiter=","), elapsed, ""


def shuffle_directory(src: Path, dst: Path, starts: np.ndarray, ends: np.ndarray,
                      seed: int = 42) -> dict:
    """Redraw every spike uniformly inside the window it came from.

    Preserves each cell's spike count per window and the window geometry; destroys every
    temporal relationship finer than a window. This is the floor the jitter curve should
    approach if the edges depend on timing at all.
    """
    dst.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    order = np.argsort(starts)
    starts, ends = starts[order], ends[order]

    kept = total = 0
    for path in sorted(src.glob("cell*.txt"), key=lambda p: int(p.stem[4:])):
        v = np.atleast_1d(np.loadtxt(path))
        total += v.size
        if v.size:
            idx = np.clip(np.searchsorted(starts, v, side="right") - 1, 0, len(starts) - 1)
            w = np.sort(rng.uniform(starts[idx], ends[idx]))
        else:
            w = v
        kept += w.size
        np.savetxt(dst / path.name, w, fmt="%.6f")
    return {"spikes_kept": kept, "spikes_total": total, "dropped_at_edges": 0,
            "retention_pct": 100.0 * kept / max(total, 1)}


def compare(base: np.ndarray, jit: np.ndarray) -> dict:
    """Edge-set and weight agreement between the baseline and a jittered run."""
    b, j = base != 0, jit != 0
    np.fill_diagonal(b, False)
    np.fill_diagonal(j, False)
    overlap = int((b & j).sum())
    union = int((b | j).sum())
    both = b & j
    sign_agree = (float((np.sign(base[both]) == np.sign(jit[both])).mean() * 100)
                  if overlap else float("nan"))
    if overlap > 2:
        r = float(np.corrcoef(base[both], jit[both])[0, 1])
    else:
        r = float("nan")
    nz = jit[j]
    return {
        "n_base": int(b.sum()), "n_jit": int(j.sum()), "overlap": overlap,
        "survival_pct": round(100.0 * overlap / max(int(b.sum()), 1), 2),
        "jaccard": round(overlap / max(union, 1), 4),
        "sign_agree_pct": round(sign_agree, 2) if overlap else "",
        "weight_r": round(r, 4) if overlap > 2 else "",
        "neg_pct": round(100.0 * float((nz < 0).mean()), 2) if nz.size else "",
    }


def run_animal(animal: str, levels: list, workroot: Path, offset_s: float, seed: int) -> list:
    old = ANIMAL_MAP[animal]
    ongoing_dir = PROJECT_ROOT / f"DATA{old}_ongoing"
    mat = next((c for c in (PROJECT_ROOT / f"DATA{old}_ongoing.mat",
                            PROJECT_ROOT / f"DATA{old}_lightON.mat") if c.exists()), None)
    if mat is None or not ongoing_dir.is_dir():
        return [{"animal": animal, "level": lv, "status": "missing",
                 "note": "no .mat or DATA dir"} for lv in levels]

    starts, ends = stimulus_windows(mat)
    durations = ends - starts
    w_start = starts + offset_s
    w_end = w_start + durations
    n_cells = len(list(ongoing_dir.glob("cell*.txt")))

    # Baseline: the aligned-ongoing arm with no displacement at all.
    staging = workroot / f"{animal}_staged"
    if staging.exists():
        shutil.rmtree(staging)
    gap_directory_windows(ongoing_dir, staging, w_start, w_end)
    base, base_elapsed, err = _glmcc(staging, n_cells)
    if base is None:
        return [{"animal": animal, "level": lv, "status": "error", "note": f"baseline: {err}"}
                for lv in levels]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    np.savetxt(OUT_DIR / f"adj_{animal}_baseline.csv", base, delimiter=",", fmt="%.6f")

    rows = []
    for lv in levels:
        row = {"animal": animal, "level": str(lv), "n_units": n_cells,
               "jitter_ms": "" if lv == "shuffle" else lv, "status": "", "note": ""}
        work = workroot / f"{animal}_jit{lv}"
        if work.exists():
            shutil.rmtree(work)

        if lv == "shuffle":
            stats = shuffle_directory(staging, work, w_start, w_end, seed)
        else:
            stats = jitter_directory(staging, work, w_start, w_end, float(lv), seed)
        row.update({k: stats[k] for k in ("spikes_kept", "spikes_total")})
        row["retention_pct"] = round(stats["retention_pct"], 3)

        W, elapsed, err = _glmcc(work, n_cells)
        row["elapsed_s"] = elapsed
        if W is None:
            row.update({"status": "error", "note": err})
            rows.append(row)
            continue
        np.savetxt(OUT_DIR / f"adj_{animal}_jit{lv}.csv", W, delimiter=",", fmt="%.6f")
        row.update({"status": "ok", **compare(base, W)})
        rows.append(row)
        print(f"  [{animal} jit={lv}] survival={row['survival_pct']}% "
              f"edges {row['n_base']}->{row['n_jit']} retention={row['retention_pct']}% "
              f"in {elapsed}s", flush=True)
        shutil.rmtree(work, ignore_errors=True)

    shutil.rmtree(staging, ignore_errors=True)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--levels", nargs="+", default=["0", "5", "10", "25", "50", "100", "shuffle"],
                    help="jitter half-widths in ms, plus the literal 'shuffle'")
    ap.add_argument("--animals", nargs="+", default=list(ANIMAL_MAP))
    ap.add_argument("--offset", type=float, default=45.0,
                    help="seconds to shift each pulse window into the light-off interval")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--workroot", type=Path,
                    default=Path("/tmp/claude-1000/jitter_sensitivity"))
    args = ap.parse_args()

    levels = [lv if lv == "shuffle" else int(lv) for lv in args.levels]
    args.workroot.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for animal in args.animals:
        print(f"\n=== {animal} ===", flush=True)
        all_rows.extend(run_animal(animal, levels, args.workroot, args.offset, args.seed))

    with SUMMARY.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)
    ok = sum(1 for r in all_rows if r.get("status") == "ok")
    print(f"\n{ok}/{len(all_rows)} runs ok. Summary: {SUMMARY}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
