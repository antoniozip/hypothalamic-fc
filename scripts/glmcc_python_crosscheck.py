#!/usr/bin/env python3
"""Run the vendored Python GLMCC and the C port on one identical fixture.

The positive control showed the C pipeline recovers 15% of known couplings at best and
returns 100% negative weights against a 50/50 excitatory/inhibitory ground truth. That
either reflects a porting bug or something intrinsic to GLMCC on this kind of data.
Running Kobayashi's own Python implementation on the same spikes separates the two.

Fixture is sized so Python (~50x slower) is feasible and so the duration stays under
Est_Data.py's hardcoded T = 5400 s, avoiding truncation on that side.
"""

import shutil, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from glmcc_positive_control import simulate, score  # noqa: E402

C_BIN = ROOT / "vendor" / "glmcc-c" / "glmcc_fixed"
PY_EST = ROOT / "vendor" / "GLMCC" / "Est_Data.py"

N, RATE, DUR = 12, 2.0, 5000.0
N_EXC, N_INH = 8, 8


def write_cells(trains, d: Path):
    d.mkdir(parents=True, exist_ok=True)
    for i, t in enumerate(trains):   # 0-indexed, as both implementations expect
        np.savetxt(d / f"cell{i}.txt", t * 1000.0, fmt="%.3f")


def main() -> int:
    rng = np.random.default_rng(42)
    trains, truth = simulate(N, RATE, DUR, N_EXC, N_INH, 0.30, 2.0, 0.5, rng)
    spk = int(np.mean([t.size for t in trains]))
    print(f"fixture: {N} neurons, {spk:,} spikes/neuron, {len(truth)} true connections "
          f"({N_EXC} E / {N_INH} I) among {N*(N-1)} ordered pairs\n")
    print(f"{'implementation':16s} {'density%':>9s} {'recall%':>8s} {'precis%':>8s} "
          f"{'sign%':>6s} {'det.neg%':>9s} {'secs':>7s}")

    work = Path("/tmp/claude-1000/glmcc_xcheck")
    if work.exists():
        shutil.rmtree(work)

    # --- C ---
    cdir = work / "c"
    write_cells(trains, cdir)
    t0 = time.time()
    p = subprocess.run([str(C_BIN), ".", str(N), "exp", "GLM", f"{DUR*1.01:.1f}"],
                       cwd=cdir, capture_output=True, text=True)
    c_secs = time.time() - t0
    got = sorted(cdir.glob("W_py_*.csv"))
    if p.returncode == 0 and got:
        s = score(np.loadtxt(got[0], delimiter=","), truth, N)
        print(f"{'C (glmcc_fixed)':16s} {s['density_pct']:9.2f} {s['recall_pct']:8.1f} "
              f"{s['precision_pct']:8.1f} {str(s['sign_correct_pct']):>6s} "
              f"{str(s['detected_neg_pct']):>9s} {c_secs:7.1f}")
    else:
        print(f"{'C':16s} FAILED: {(p.stderr or '')[:70]}")

    # --- Python (Kobayashi reference) ---
    pdir = work / "py"
    write_cells(trains, pdir)
    t0 = time.time()
    p = subprocess.run([sys.executable, str(PY_EST), ".", str(N), "exp", "GLM"],
                       cwd=pdir, capture_output=True, text=True)
    py_secs = time.time() - t0
    got = sorted(pdir.glob("W_py_*.csv"))
    if p.returncode == 0 and got:
        W = np.loadtxt(got[0], delimiter=",")
        s = score(W, truth, N)
        print(f"{'Python (vendor)':16s} {s['density_pct']:9.2f} {s['recall_pct']:8.1f} "
              f"{s['precision_pct']:8.1f} {str(s['sign_correct_pct']):>6s} "
              f"{str(s['detected_neg_pct']):>9s} {py_secs:7.1f}")
    else:
        tail = ((p.stderr or "") + (p.stdout or "")).strip().splitlines()
        print(f"{'Python':16s} FAILED after {py_secs:.1f}s: {tail[-1][:90] if tail else '?'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
