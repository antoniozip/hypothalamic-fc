#!/usr/bin/env python3
"""Estimate the effective high-pass corner per recording, and harmonise the two blocks.

Context
-------
The nominal acquisition filter was 600-6000 Hz for every session. The waveforms do not agree
with that. Median peak-to-peak duration is 0.500 ms in the seven recordings from 2018-01-31
onward -- consistent with a 600 Hz high-pass -- but 0.775 ms (up to 1.05 ms) in the six earlier
ones, which is what a substantially *lower* corner produces. Trough depth and positive overshoot
shift across the same boundary and in the same direction (see
`reports/ei_differential/waveform-diagnostic.md`).

Since the pre-sorting data no longer exists, the raw traces cannot be re-filtered. What can be
done is to apply a compensating high-pass to the stored mean waveforms, so that both blocks sit
at a comparable effective passband before classification.

Method
------
1. For each recording, sweep a compensating 2nd-order Butterworth high-pass over a range of
   corners, applied zero-phase (filtfilt, so the trough does not move), and find the corner that
   brings that recording's median PP to the late-block target.
2. Apply it, recompute FW and PP, and report how much the between-recording spread shrinks.

Caveats, which are serious
--------------------------
Filtering an already-filtered *mean* waveform is not equivalent to filtering the raw trace: the
original filter's phase response is already baked in, and averaging across spikes has smoothed
it. The record is also short -- 81 samples at 20 kHz is 4.05 ms, only ~2.4 periods of a 600 Hz
component -- so edge effects matter; reflection padding is used to limit them. This is a
harmonisation for the purpose of making a threshold comparable, not a reconstruction of what a
correctly filtered recording would have looked like.

Output: data/processed/waveform_features_harmonised.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
from scipy.signal import butter, filtfilt

ROOT = Path(__file__).resolve().parent.parent
FS = 20000.0
PAD = 120                      # reflection padding, samples
TARGET_PP_MS = 0.500           # late-block median, the block matching the nominal 600 Hz
AM = {"day1": "171019", "night1": "171207", "night2": "171208", "night3": "171213",
      "night4": "180110", "night5": "180111", "day2": "180131", "night6": "180221",
      "night7": "180228", "day3": "180302", "day4": "180419", "day5": "180420",
      "day6": "180423"}
EARLY = {"day1", "night1", "night2", "night3", "night4", "night5"}


def feats(wf):
    """(FW ms, PP ms) measured on the negative trough; NaN if unmeasurable."""
    wf = np.asarray(wf, float)
    tr = int(np.argmin(wf))
    base = float(np.median(wf[:max(tr // 2, 3)])) if tr >= 6 else float(np.median(wf))
    depth = base - wf[tr]
    if depth <= 0 or tr >= wf.size - 1:
        return np.nan, np.nan
    pk = tr + int(np.argmax(wf[tr:]))
    pp = (pk - tr) / FS * 1000 if pk > tr else np.nan
    half = base - depth / 2.0
    lo = hi = None
    for i in range(tr, 0, -1):
        if (wf[i] - half) * (wf[i - 1] - half) <= 0:
            d = wf[i - 1] - wf[i]
            lo = (i - 1) + (wf[i - 1] - half) / d if d else float(i)
            break
    for i in range(tr, wf.size - 1):
        if (wf[i] - half) * (wf[i + 1] - half) <= 0:
            d = wf[i + 1] - wf[i]
            hi = i + (half - wf[i]) / d if d else float(i)
            break
    fw = (hi - lo) / FS * 1000 if (lo is not None and hi is not None and hi > lo) else np.nan
    return fw, pp


def highpass(wf, fc):
    """Zero-phase 2nd-order Butterworth high-pass with reflection padding."""
    if fc <= 0:
        return wf
    b, a = butter(2, fc / (FS / 2), btype="high")
    p = np.r_[wf[PAD:0:-1], wf, wf[-2:-PAD - 2:-1]]
    return filtfilt(b, a, p)[len(wf[PAD:0:-1]): len(wf[PAD:0:-1]) + len(wf)]


def main() -> int:
    raw = {}
    for a, o in AM.items():
        U = scipy.io.loadmat(ROOT / f"DATA{o}.mat", squeeze_me=True,
                             struct_as_record=False)["Units"]
        raw[a] = [np.asarray(u.wf, float) for u in U]

    grid = np.r_[0, np.arange(50, 1550, 25)]
    rows, chosen = [], {}
    print(f"{'animal':8s} {'block':6s} {'PP before':>10s} {'corner Hz':>10s} {'PP after':>9s}")
    for a in AM:
        pp0 = np.nanmedian([feats(w)[1] for w in raw[a]])
        best, bestd = 0.0, abs(pp0 - TARGET_PP_MS)
        for fc in grid[1:]:
            m = np.nanmedian([feats(highpass(w, fc))[1] for w in raw[a]])
            if np.isfinite(m) and abs(m - TARGET_PP_MS) < bestd:
                best, bestd = float(fc), abs(m - TARGET_PP_MS)
        chosen[a] = best
        pp1 = np.nanmedian([feats(highpass(w, best))[1] for w in raw[a]]) if best else pp0
        print(f"{a:8s} {'early' if a in EARLY else 'late':6s} {pp0:10.3f} "
              f"{best:10.0f} {pp1:9.3f}")
        for i, w in enumerate(raw[a]):
            fw, pp = feats(highpass(w, best) if best else w)
            rows.append({"animal": a, "neuron_id": i + 1,
                         "block": "early" if a in EARLY else "late",
                         "corner_hz": best, "fw_ms": fw, "pp_ms": pp,
                         "usable": bool(np.isfinite(fw) and np.isfinite(pp))})

    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "data" / "processed" / "waveform_features_harmonised.csv", index=False)

    g = df[df.usable]
    med = g.groupby("animal").pp_ms.median()
    before = pd.Series({a: np.nanmedian([feats(w)[1] for w in raw[a]]) for a in AM})
    print(f"\nbetween-recording median PP spread:")
    print(f"  before  {before.min():.3f}-{before.max():.3f} ms  ({before.max()/before.min():.2f}x)")
    print(f"  after   {med.min():.3f}-{med.max():.3f} ms  ({med.max()/med.min():.2f}x)")
    print(f"\ncompensating corners: early {[int(chosen[a]) for a in AM if a in EARLY]}, "
          f"late {[int(chosen[a]) for a in AM if a not in EARLY]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
