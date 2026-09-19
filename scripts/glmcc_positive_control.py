#!/usr/bin/env python3
"""Positive control: can GLMCC recover known couplings at this dataset's operating point?

The regenerated connectivity has edge count tracking spike count (rho=+0.94) and comes out
99.8-100% inhibitory, so it is unclear whether the estimator is measuring circuitry at all.
This injects KNOWN couplings into independent Poisson trains and asks three questions:

  recall     what fraction of true connections are recovered
  precision  what fraction of reported connections are real
  sign       of the recovered ones, how many get excitatory/inhibitory right

Scenarios match the real recording regimes, so the answer is specific to this data:

  ongoing-like    1.4 Hz over 15,000 s  (~21,000 spikes/neuron)
  lightON-like    0.65 Hz over 1,000 s  (~650 spikes/neuron)  <- data-starved
  generous        2.0 Hz over 15,000 s  (best case)

Ground truth uses the standard injection model: an excitatory connection copies a presynaptic
spike to the postsynaptic train after a short latency with probability p; an inhibitory one
deletes postsynaptic spikes in the latency window with probability p.
"""

import argparse, csv, json, shutil, subprocess, tempfile, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "vendor" / "glmcc-c" / "glmcc_fixed"


def simulate(n_neurons, rate_hz, duration_s, n_exc, n_inh, p_transmit, delay_ms,
             jitter_ms, rng):
    """Independent Poisson trains plus injected couplings. Times in seconds."""
    trains = [np.sort(rng.uniform(0, duration_s, rng.poisson(rate_hz * duration_s)))
              for _ in range(n_neurons)]

    pairs = [(i, j) for i in range(n_neurons) for j in range(n_neurons) if i != j]
    chosen = rng.choice(len(pairs), size=n_exc + n_inh, replace=False)
    truth = {}
    for k, idx in enumerate(chosen):
        truth[pairs[idx]] = 1 if k < n_exc else -1

    for (pre, post), sign in truth.items():
        src = trains[pre]
        if src.size == 0:
            continue
        fire = rng.random(src.size) < p_transmit
        lat = (delay_ms + rng.normal(0, jitter_ms, src.size)) / 1000.0
        if sign > 0:
            added = src[fire] + lat[fire]
            added = added[(added >= 0) & (added < duration_s)]
            trains[post] = np.sort(np.concatenate([trains[post], added]))
        else:
            tgt = trains[post]
            if tgt.size == 0:
                continue
            lo = src[fire] + delay_ms / 1000.0
            hi = lo + 3.0 / 1000.0          # 3 ms suppression window
            idxs = np.searchsorted(tgt, lo)
            drop = np.zeros(tgt.size, bool)
            for a, b in zip(idxs, hi):
                while a < tgt.size and tgt[a] <= b:
                    drop[a] = True
                    a += 1
            trains[post] = tgt[~drop]
    return trains, truth


def run_glmcc(trains, duration_s, workdir):
    workdir.mkdir(parents=True, exist_ok=True)
    for i, t in enumerate(trains, start=1):
        np.savetxt(workdir / f"cell{i}.txt", t * 1000.0, fmt="%.3f")   # 0-based ms
    proc = subprocess.run([str(BIN), ".", str(len(trains)), "exp", "GLM",
                           f"{duration_s * 1.01:.1f}"],
                          cwd=workdir, capture_output=True, text=True)
    produced = sorted(workdir.glob("W_py_*.csv"))
    if proc.returncode != 0 or not produced:
        return None
    return np.loadtxt(produced[0], delimiter=",")


def score(W, truth, n):
    detected = {(i, j): W[i, j] for i in range(n) for j in range(n)
                if i != j and W[i, j] != 0}
    tp = [p for p in detected if p in truth]
    fp = [p for p in detected if p not in truth]
    correct_sign = sum(1 for p in tp if np.sign(detected[p]) == truth[p])
    n_pairs = n * (n - 1)
    return {
        "n_true": len(truth), "n_detected": len(detected),
        "density_pct": round(100 * len(detected) / n_pairs, 2),
        "recall_pct": round(100 * len(tp) / max(len(truth), 1), 1),
        "precision_pct": round(100 * len(tp) / max(len(detected), 1), 1),
        "sign_correct_pct": round(100 * correct_sign / max(len(tp), 1), 1) if tp else "",
        "detected_neg_pct": round(100 * np.mean([v < 0 for v in detected.values()]), 1)
                            if detected else "",
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--neurons", type=int, default=30)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    scenarios = [
        ("generous",     2.00, 15000.0),
        ("ongoing-like", 1.40, 15000.0),
        ("lightON-like", 0.65,  1000.0),
    ]
    n = args.neurons
    rows = []
    print(f"{n} neurons, 20 excitatory + 20 inhibitory injected connections "
          f"({100*40/(n*(n-1)):.1f}% of pairs), p_transmit=0.3, latency 2+/-0.5 ms\n")
    print(f"{'scenario':14s} {'rate':>5s} {'dur_s':>7s} {'spk/nrn':>8s} {'density%':>9s} "
          f"{'recall%':>8s} {'precis%':>8s} {'sign%':>6s} {'det.neg%':>9s}")

    with tempfile.TemporaryDirectory(dir="/tmp/claude-1000") as tmp:
        for name, rate, dur in scenarios:
            trains, truth = simulate(n, rate, dur, 20, 20, 0.30, 2.0, 0.5, rng)
            spk = int(np.mean([t.size for t in trains]))
            t0 = time.time()
            W = run_glmcc(trains, dur, Path(tmp) / name)
            if W is None:
                print(f"{name:14s} GLMCC failed"); continue
            s = score(W, truth, n)
            s.update(scenario=name, rate_hz=rate, duration_s=dur,
                     spikes_per_neuron=spk, elapsed_s=round(time.time() - t0, 1))
            rows.append(s)
            print(f"{name:14s} {rate:5.2f} {dur:7.0f} {spk:8,d} {s['density_pct']:9.2f} "
                  f"{s['recall_pct']:8.1f} {s['precision_pct']:8.1f} "
                  f"{str(s['sign_correct_pct']):>6s} {str(s['detected_neg_pct']):>9s}")

    out = ROOT / "results" / "glmcc_positive_control.csv"
    if rows:
        with out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
        print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
