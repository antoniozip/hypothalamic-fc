#!/usr/bin/env python3
"""Rate-matched lightON vs light-off contrast.

The problem
-----------
GLMCC reports a coupling when |J| exceeds `J_min = sqrt(16.3 / tau / cc_0)`, and `cc_0`, the
fitted correlogram height, scales with the pair's coincidence count and so with
`n_i * n_j / T_observed`. More spikes therefore lower the threshold and admit more edges,
independently of any change in circuitry. Across the 26 main matrices, Spearman rho between
spikes-per-unit and density is **+0.955** (p = 3.7e-14), and `ongoing` carries 10.2-13.7x more
spikes per unit than `lightON`. Any `lightON` vs `ongoing` comparison is therefore confounded
with detection power at its root.

The contrast built here
-----------------------
Two arms that differ only in whether the light was on:

  lightON    DATA{old}, the evoked epoch, inside its true stimulus windows
  lightOFF   DATA{old}_ongoing, gapped to those same windows shifted --offset seconds
             into the light-off interval

Same window count, same window durations, same period: observation geometry is identical by
construction, so `T_observed` matches. On top of that, each unit is thinned to the smaller of its
two spike counts, so `n_i` matches per unit as well. Both terms of `n_i * n_j / T_observed` are
then equal between arms and `J_min` is on the same footing.

Thinning is uniform random sampling without replacement. It preserves the *shape* of the
cross-correlogram -- a coincidence survives with probability p_i * p_j regardless of lag -- while
lowering the counts, which is exactly the intent: equalise detection power without inventing or
destroying temporal structure.

Arms are written as `validated_adj_{animal}_{lightON,ongoing}.csv` under `results/glmcc_rm/` so
the existing R pipeline can consume them unchanged. **"ongoing" in that directory means the
rate-matched light-off control, not the raw ongoing epoch.**

Output: results/glmcc_rm/validated_adj_{animal}_{lightON,ongoing}.csv (+ .meta.json)
        results/rate_matched_summary.csv
"""

import argparse
import csv
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from gap_spike_trains import gap_directory_windows, stimulus_windows  # noqa: E402
from glmcc_units import run_glmcc  # noqa: E402

OUT = PROJECT_ROOT / "results" / "glmcc_rm"
SUMMARY = PROJECT_ROOT / "results" / "rate_matched_summary.csv"

ANIMAL_MAP = {
    "day1": "171019", "day2": "180131", "day3": "180302",
    "day4": "180419", "day5": "180420", "day6": "180423",
    "night1": "171207", "night2": "171208", "night3": "171213",
    "night4": "180110", "night5": "180111", "night6": "180221", "night7": "180228",
}

FIELDS = ["animal", "condition", "status", "n_units", "n_windows", "observed_s",
          "spikes_before", "spikes_after", "thinned_pct", "nonzero", "density_pct",
          "neg_pct", "median_w", "elapsed_s", "note"]


def counts(d: Path) -> np.ndarray:
    return np.array([np.atleast_1d(np.loadtxt(p)).size
                     for p in sorted(d.glob("cell*.txt"), key=lambda q: int(q.stem[4:]))])


def thin_to(src: Path, dst: Path, target: np.ndarray, seed: int) -> int:
    """Uniformly subsample each cell down to target[i] spikes, without replacement."""
    dst.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    kept = 0
    for i, p in enumerate(sorted(src.glob("cell*.txt"), key=lambda q: int(q.stem[4:]))):
        v = np.atleast_1d(np.loadtxt(p))
        t = int(target[i])
        if v.size > t:
            v = np.sort(rng.choice(v, size=t, replace=False))
        kept += v.size
        np.savetxt(dst / p.name, v, fmt="%.6f")
    return kept


def run_animal(animal: str, workroot: Path, offset_s: float, seed: int) -> list:
    old = ANIMAL_MAP[animal]
    light_dir = PROJECT_ROOT / f"DATA{old}"
    ongoing_dir = PROJECT_ROOT / f"DATA{old}_ongoing"
    mat = next((c for c in (PROJECT_ROOT / f"DATA{old}_ongoing.mat",
                            PROJECT_ROOT / f"DATA{old}_lightON.mat") if c.exists()), None)
    if mat is None or not light_dir.is_dir() or not ongoing_dir.is_dir():
        return [{"animal": animal, "condition": c, "status": "missing",
                 "note": "missing .mat or DATA dir"} for c in ("lightON", "ongoing")]

    starts, ends = stimulus_windows(mat)
    durations = ends - starts
    work = workroot / animal
    shutil.rmtree(work, ignore_errors=True)

    # Both arms restricted to identical observation geometry.
    on_raw, off_raw = work / "on_raw", work / "off_raw"
    s_on = gap_directory_windows(light_dir, on_raw, starts, ends)
    s_off = gap_directory_windows(ongoing_dir, off_raw, starts + offset_s,
                                  starts + offset_s + durations)

    c_on, c_off = counts(on_raw), counts(off_raw)
    if c_on.size != c_off.size:
        return [{"animal": animal, "condition": c, "status": "error",
                 "note": f"unit count differs: {c_on.size} vs {c_off.size}"}
                for c in ("lightON", "ongoing")]

    target = np.minimum(c_on, c_off)          # per-unit common floor
    on_m, off_m = work / "on", work / "off"
    k_on = thin_to(on_raw, on_m, target, seed)
    k_off = thin_to(off_raw, off_m, target, seed + 1)

    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for cond, src, before, after, stats in (("lightON", on_m, int(c_on.sum()), k_on, s_on),
                                            ("ongoing", off_m, int(c_off.sum()), k_off, s_off)):
        row = {"animal": animal, "condition": cond, "n_units": int(c_on.size),
               "n_windows": len(starts), "observed_s": round(stats["observed_s"], 1),
               "spikes_before": before, "spikes_after": after,
               "thinned_pct": round(100 * (1 - after / max(before, 1)), 2),
               "status": "", "note": ""}
        W, elapsed, err = run_glmcc(src, int(c_on.size))
        row["elapsed_s"] = elapsed
        if W is None:
            row.update(status="error", note=err)
            rows.append(row)
            continue
        out = OUT / f"validated_adj_{animal}_{cond}.csv"
        np.savetxt(out, W, delimiter=",", fmt="%.6f")
        n = W.shape[0]
        off_diag = ~np.eye(n, dtype=bool)
        nz = W[off_diag][W[off_diag] != 0]
        row.update({
            "status": "ok", "nonzero": len(nz),
            "density_pct": round(100 * len(nz) / max(n * (n - 1), 1), 2),
            "neg_pct": round(100 * float((nz < 0).mean()), 1) if len(nz) else "",
            "median_w": round(float(np.median(nz)), 4) if len(nz) else "",
        })
        # The R pipeline requires a sidecar naming the encoding.
        out.with_suffix(".meta.json").write_text(json.dumps({
            "encoding": "glmcc_weight_masked",
            "note": ("rate- and geometry-matched arms; 'ongoing' here is the light-off "
                     "control, not the raw ongoing epoch"),
            "condition": cond, "offset_s": offset_s, "seed": seed,
            "n_windows": len(starts), "spikes_before": before, "spikes_after": after,
            "matching": "per-unit thinning to min(lightON, lightOFF) count",
            "written_at": datetime.now(timezone.utc).isoformat(),
        }, indent=2))
        rows.append(row)
        print(f"  [{animal}/{cond}] {after:,} spikes (thinned {row['thinned_pct']}%) "
              f"-> {row['nonzero']} edges, {row['density_pct']}% dense, "
              f"{row['neg_pct']}% neg, {elapsed}s", flush=True)

    shutil.rmtree(work, ignore_errors=True)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--animals", nargs="+", default=list(ANIMAL_MAP))
    ap.add_argument("--offset", type=float, default=45.0)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--workroot", type=Path, default=Path("/tmp/claude-1000/rate_matched"))
    args = ap.parse_args()
    args.workroot.mkdir(parents=True, exist_ok=True)

    rows = []
    for a in args.animals:
        print(f"\n=== {a} ===", flush=True)
        rows.extend(run_animal(a, args.workroot, args.offset, args.seed))

    with SUMMARY.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    ok = sum(1 for r in rows if r.get("status") == "ok")
    print(f"\n{ok}/{len(rows)} arms ok. Summary: {SUMMARY}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
