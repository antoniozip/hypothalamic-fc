#!/usr/bin/env python3
"""Compare GLMCC weight sign distributions across the three conditions.

lightON        - native ~10% duty cycle (stimulus windows)
ongoing        - native ~89% duty cycle (near-continuous)
ongoing_gapped - ongoing re-run after imposing that animal's own lightON gap structure

If the sign asymmetry is an artifact of the duty cycle, ongoing_gapped should move away from
ongoing and toward lightON.
"""

import statistics
from pathlib import Path

import numpy as np

RESULTS = Path(__file__).resolve().parent.parent / "results" / "glmcc"
ANIMALS = ["day1", "day2", "day3", "day4", "day5", "day6",
           "night1", "night2", "night3", "night4", "night5", "night6", "night7"]
CONDITIONS = [("lightON", "adj_{a}_lightON.csv"),
              ("ongoing", "adj_{a}_ongoing.csv"),
              ("ongoing_gapped", "adj_{a}_ongoing_gapped.csv"),
              ("ongoing_aligned", "adj_{a}_ongoing_aligned.csv"),
              ("lightON_jit", "adj_{a}_lightON_jittered.csv"),
              ("ongoing_al_jit", "adj_{a}_ongoing_aligned_jittered.csv")]


def summarise(path: Path) -> dict | None:
    if not path.exists():
        return None
    m = np.loadtxt(path, delimiter=",", ndmin=2)
    n = m.shape[0]
    off = ~np.eye(n, dtype=bool)
    nz = m[off][m[off] != 0]
    if not len(nz):
        return {"n": n, "nonzero": 0, "neg": float("nan"), "median": float("nan")}
    return {"n": n, "nonzero": len(nz),
            "neg": 100 * float((nz < 0).mean()), "median": float(np.median(nz))}


def main() -> int:
    print(f"{'animal':8s} " + " ".join(f"{c:>22s}" for c, _ in CONDITIONS))
    print(f"{'':8s} " + " ".join(f"{'neg%   median  n':>22s}" for _ in CONDITIONS))
    agg = {c: {"neg": [], "med": []} for c, _ in CONDITIONS}

    for a in ANIMALS:
        cells = []
        for cond, pattern in CONDITIONS:
            s = summarise(RESULTS / pattern.format(a=a))
            if s is None or s["nonzero"] == 0:
                cells.append(f"{'-':>22s}")
                continue
            agg[cond]["neg"].append(s["neg"])
            agg[cond]["med"].append(s["median"])
            cells.append(f"{s['neg']:6.1f} {s['median']:8.3f} {s['nonzero']:6d}")
        print(f"{a:8s} " + " ".join(cells))

    print()
    print(f"{'condition':16s} {'median neg%':>12s} {'range':>16s} {'median weight':>14s}  n")
    for cond, _ in CONDITIONS:
        negs, meds = agg[cond]["neg"], agg[cond]["med"]
        if not negs:
            continue
        print(f"{cond:16s} {statistics.median(negs):11.1f}% "
              f"{min(negs):6.1f}-{max(negs):5.1f}% {statistics.median(meds):14.3f}  {len(negs)}")

    if agg['ongoing']['neg'] and agg['ongoing_gapped']['neg']:
        drops = [a - b for a, b in zip(agg['ongoing']['neg'], agg['ongoing_gapped']['neg'])]
        print(f"\nPer-animal drop in negative fraction, ongoing -> ongoing_gapped:")
        print(f"  median {statistics.median(drops):.1f} points, "
              f"range {min(drops):.1f} to {max(drops):.1f}, all {len(drops)} animals "
              f"{'decrease' if all(d > 0 for d in drops) else 'MIXED'}")

    # Does using the real stimulus geometry instead of a fitted period change the answer?
    g, al = agg['ongoing_gapped']['neg'], agg['ongoing_aligned']['neg']
    if g and al and len(g) == len(al):
        diffs = [b - a for a, b in zip(g, al)]
        absd = [abs(d) for d in diffs]
        print(f"\nFitted-period mask vs real-stimulus-aligned mask (per animal, neg% points):")
        print(f"  mean signed diff {statistics.mean(diffs):+.2f}, "
              f"median |diff| {statistics.median(absd):.2f}, max |diff| {max(absd):.2f}")
        try:
            from scipy.stats import wilcoxon
            w = wilcoxon(g, al)
            print(f"  Wilcoxon signed-rank: W={w.statistic:.1f} p={w.pvalue:.3f}")
        except Exception as exc:  # scipy optional
            print(f"  (wilcoxon unavailable: {exc})")
        span = statistics.median(agg['ongoing']['neg']) - statistics.median(agg['lightON']['neg'])
        print(f"  for scale, the ongoing->lightON gap being explained is {span:.1f} points")

    # Jitter control: does destroying <25 ms structure collapse lightON onto the light-off control?
    lt, ltj = agg['lightON']['neg'], agg['lightON_jit']['neg']
    al, alj = agg['ongoing_aligned']['neg'], agg['ongoing_al_jit']['neg']
    if lt and ltj and al and alj and len({len(lt), len(ltj), len(al), len(alj)}) == 1:
        print("\nJitter control (+/-25 ms displacement, destroys structure finer than that):")
        d_light = [b - a for a, b in zip(lt, ltj)]
        d_ctrl = [b - a for a, b in zip(al, alj)]
        print(f"  lightON        {statistics.median(lt):5.1f}% -> {statistics.median(ltj):5.1f}% "
              f"(median shift {statistics.median(d_light):+5.1f} pts)")
        print(f"  light-off ctrl {statistics.median(al):5.1f}% -> {statistics.median(alj):5.1f}% "
              f"(median shift {statistics.median(d_ctrl):+5.1f} pts)")
        gap_before = [a - b for a, b in zip(al, lt)]
        gap_after = [a - b for a, b in zip(alj, ltj)]
        print(f"  control-minus-lightON gap: {statistics.median(gap_before):+5.1f} pts before "
              f"-> {statistics.median(gap_after):+5.1f} pts after jitter")
        closed = 100 * (1 - abs(statistics.median(gap_after)) / abs(statistics.median(gap_before)))
        print(f"  => jitter closes {closed:.0f}% of the lightON/control separation")
        try:
            from scipy.stats import wilcoxon
            print(f"  Wilcoxon lightON vs lightON_jittered: p={wilcoxon(lt, ltj).pvalue:.4f}")
            print(f"  Wilcoxon control vs control_jittered: p={wilcoxon(al, alj).pvalue:.4f}")
            print(f"  Wilcoxon lightON_jit vs control_jit:  p={wilcoxon(ltj, alj).pvalue:.4f}")
        except Exception as exc:
            print(f"  (wilcoxon unavailable: {exc})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
