"""Parity test: the C port of GLMCC must reproduce Kobayashi's Python reference.

The C port is only usable if it computes the same connectivity matrix as
`vendor/GLMCC/Est_Data.py`, which is the published implementation. This runs both
on one identical fixture and compares the weight matrices element-wise.

It also checks the sign of the recovered couplings against injected ground truth,
because the failure this test was written for -- the port omitted the two coupling
sufficient statistics from the design vector, so its log posterior had no term
rewarding J > 0 -- showed up as every weight coming out inhibitory.

Run: python3 tests/test_glmcc_c_parity.py
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
C_BIN = PROJECT_ROOT / "vendor" / "glmcc-c" / "glmcc_fixed"
PY_DIR = PROJECT_ROOT / "vendor" / "GLMCC"
PY_EST = PY_DIR / "Est_Data.py"

# Est_Data.py hard-codes T = 5400 s and names its output W_py_5400.csv. Keeping the
# fixture shorter than that, and passing the same T to the C port, makes both sides
# apply an identical spike filter and write an identically named file.
T_SEC = 5400.0

N_NEURONS = 5
RATE_HZ = 6.0
DURATION_S = 2000.0
DELAY_MS = 3.0
P_TRANSMIT = 0.40
# The inhibitory pair needs a stronger coupling than the excitatory one to clear GLMCC's
# detection threshold; below this the connection is missed by BOTH implementations, which
# leaves the negative branch of calc_PSP untested.
P_SUPPRESS = 0.80
SUPPRESS_MS = 4.0

# Weights agree to ~1e-5 relative; the two implementations differ only in
# floating-point accumulation order and in scipy's vs. this port's Ei().
ATOL = 2e-3
RTOL = 1e-3


def build_fixture(rng):
    """Poisson trains with two injected couplings: 0 -> 1 excitatory, 2 -> 3 inhibitory."""
    trains = [
        np.sort(rng.uniform(0, DURATION_S, rng.poisson(RATE_HZ * DURATION_S)))
        for _ in range(N_NEURONS)
    ]
    truth = {(0, 1): 1, (2, 3): -1}

    # excitatory: copy presynaptic spikes forward by the synaptic delay
    src = trains[0]
    fire = rng.random(src.size) < P_TRANSMIT
    lat = (DELAY_MS + rng.normal(0, 0.4, src.size)) / 1000.0
    added = src[fire] + lat[fire]
    trains[1] = np.sort(
        np.concatenate([trains[1], added[(added >= 0) & (added < DURATION_S)]])
    )

    # inhibitory: delete postsynaptic spikes in a short window after each presynaptic spike
    src, tgt = trains[2], trains[3]
    fire = rng.random(src.size) < P_SUPPRESS
    lo = src[fire] + DELAY_MS / 1000.0
    drop = np.zeros(tgt.size, bool)
    for a, b in zip(np.searchsorted(tgt, lo), lo + SUPPRESS_MS / 1000.0):
        while a < tgt.size and tgt[a] <= b:
            drop[a] = True
            a += 1
    trains[3] = tgt[~drop]

    return trains, truth


def write_cells(trains, d: Path) -> None:
    d.mkdir(parents=True, exist_ok=True)
    for i, t in enumerate(trains):
        np.savetxt(d / f"cell{i}.txt", t * 1000.0, fmt="%.4f")  # 0-indexed, milliseconds


def run_c(trains, d: Path, mode: str, method: str) -> np.ndarray:
    write_cells(trains, d)
    proc = subprocess.run(
        [str(C_BIN), ".", str(N_NEURONS), mode, method, f"{T_SEC:.1f}"],
        cwd=d,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"C port failed: {proc.stderr[-400:]}")
    return np.loadtxt(d / f"W_py_{T_SEC:.0f}.csv", delimiter=",")


def run_python(trains, d: Path, mode: str, method: str) -> np.ndarray:
    write_cells(trains, d)
    env = {"PYTHONPATH": str(PY_DIR), "PATH": "/usr/bin:/bin"}
    proc = subprocess.run(
        [sys.executable, str(PY_EST), ".", str(N_NEURONS), mode, method],
        cwd=d,
        capture_output=True,
        text=True,
        env=env,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Python reference failed: {(proc.stderr + proc.stdout)[-400:]}")
    return np.loadtxt(d / f"W_py_{T_SEC:.0f}.csv", delimiter=",")


# All four combinations Est_Data.py accepts. "exp" scans delays 1-4 ms and picks the best;
# "sim" fixes the delay at 3 ms. "LR" replaces the J_min threshold with a likelihood-ratio
# test, which is a separate code path in both implementations.
MODES = [("exp", "GLM"), ("sim", "GLM"), ("exp", "LR"), ("sim", "LR")]


def compare_one(trains, truth, mode: str, method: str) -> bool:
    """Run both implementations under one mode and compare."""
    work = Path(tempfile.mkdtemp(prefix="glmcc_parity_"))
    try:
        W_c = run_c(trains, work / "c", mode, method)
        W_py = run_python(trains, work / "py", mode, method)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    ok = True
    print(f"\n-- {mode}/{method} " + "-" * 52)

    diff = np.abs(W_c - W_py)
    worst = np.unravel_index(np.argmax(diff), diff.shape)
    if np.allclose(W_c, W_py, atol=ATOL, rtol=RTOL):
        print(f"PASS  weight matrices agree (max |diff| = {diff.max():.2e})")
    else:
        ok = False
        print(
            f"FAIL  weight matrices differ: max |diff| = {diff.max():.4g} at {worst}, "
            f"C = {W_c[worst]:.6f}, Python = {W_py[worst]:.6f}"
        )

    # Both sides must agree on which entries are non-zero (the detection decision).
    if np.array_equal(W_c != 0, W_py != 0):
        print(f"PASS  identical support ({int((W_c != 0).sum())} non-zero entries)")
    else:
        ok = False
        n_only_c = int(((W_c != 0) & (W_py == 0)).sum())
        n_only_py = int(((W_py != 0) & (W_c == 0)).sum())
        print(f"FAIL  support differs: {n_only_c} only in C, {n_only_py} only in Python")

    # Guards the original defect: a port that drives every J to the negative clamp
    # still produces a plausible-looking matrix, but it is uniformly inhibitory.
    nonzero = W_c[W_c != 0]
    if nonzero.size and (nonzero < 0).all():
        ok = False
        print(f"FAIL  all {nonzero.size} detected C weights are negative (sign inversion)")
    else:
        print(f"PASS  C weights are not uniformly negative ({int((nonzero > 0).sum())} positive)")

    for (pre, post), sign in sorted(truth.items()):
        got = W_c[pre, post]
        mark = "PASS" if np.sign(got) == sign else ("MISS" if got == 0 else "FAIL")
        if mark == "FAIL":
            ok = False
        kind = "excitatory" if sign > 0 else "inhibitory"
        print(f"{mark}  ground truth {pre}->{post} {kind}: C = {got:+.4f}, Python = {W_py[pre, post]:+.4f}")

    return ok


def test_c_matches_python_reference() -> bool:
    """The C port reproduces the reference weight matrix in every mode."""
    if not C_BIN.exists():
        print(f"FAIL  C binary not built: {C_BIN}")
        print("      build with: gcc -O3 -fopenmp -o glmcc_fixed glmcc_fixed.c -lm")
        return False

    rng = np.random.default_rng(20260919)
    trains, truth = build_fixture(rng)
    print(
        f"fixture: {N_NEURONS} neurons, "
        f"{int(np.mean([t.size for t in trains])):,} spikes/neuron, "
        f"{DURATION_S:.0f} s, couplings {sorted(truth)}"
    )

    return all([compare_one(trains, truth, m, meth) for m, meth in MODES])


if __name__ == "__main__":
    raise SystemExit(0 if test_c_matches_python_reference() else 1)
