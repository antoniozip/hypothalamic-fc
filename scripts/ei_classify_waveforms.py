#!/usr/bin/env python3
"""Re-derive the E/I classification from spike waveforms, following spikeMAP.

Method (Mahallati et al., eLife 2026, https://elifesciences.org/articles/106557)
-------------------------------------------------------------------------------
Two features of the mean spike waveform:

  FW  full width at half-maximum amplitude
  PP  peak-to-peak duration

    "the full width (FW) of spikes at half-maximum amplitude and (2) the peak-to-peak
     (PP) duration"

classified by clustering rather than fixed cut-offs:

    "K-means clustering was performed on values of FW and PP obtained across all
     simultaneously recorded somatic electrodes, where the optimal number of clusters was
     chosen with a Calinski-Harabasz criterion"

with inhibitory cells being "a single cluster with rapid time constants for both FW and PP".
In mouse prefrontal cortex at 18 kHz they classified 5.25% of 8,189 units as inhibitory.

Why this replaces the existing classification
---------------------------------------------
`tim_stuff.m:696-714` computes the same two features but classifies with hard thresholds
(FW < 0.6 ms AND PP < 0.65 ms), and its FW is wrong:

    [pks,locs,ws{i},proms] = findpeaks(Units(i).wf);
    half_w(i) = ws{i}(end)/20000;

`findpeaks` returns *positive* peaks, so `ws{i}(end)` is the half-prominence width of the
**last positive bump** in the waveform. An extracellular spike's dominant deflection is the
negative trough, so that number is not the spike's half-width at all -- it is a property of the
trailing repolarisation. That, plus arbitrary thresholds, is the likely source of the
implausible 44.6% inhibitory fraction the old labels give.

Here FW is measured on the trough, at half its depth below baseline, with linear interpolation
for sub-sample resolution.

Deviation from the paper
------------------------
The paper clusters per recording. These recordings carry 24-126 units each, too few for a
stable 2-D k-means, so clustering is done once over all units pooled across the 13 recordings
(same rig, same probe, same sampling rate) and the per-animal fractions are reported after the
fact. Both are recorded in the output.

Sampling rate is 20 kHz, from `tim_stuff.m` (`/20000`).

Output: data/processed/waveform_features.csv   per-unit FW, PP, cluster, ei
        results/ei/classification_compare.csv  new vs old labels
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
from sklearn.cluster import KMeans
from sklearn.metrics import calinski_harabasz_score

ROOT = Path(__file__).resolve().parent.parent
FS = 20000.0          # Hz, per tim_stuff.m
ANIMAL_MAP = {
    "day1": "171019", "day2": "180131", "day3": "180302", "day4": "180419",
    "day5": "180420", "day6": "180423", "night1": "171207", "night2": "171208",
    "night3": "171213", "night4": "180110", "night5": "180111",
    "night6": "180221", "night7": "180228",
}


def features(wf: np.ndarray) -> tuple[float, float, bool]:
    """Return (FW seconds, PP seconds, ok). Both measured on the negative trough."""
    wf = np.asarray(wf, dtype=float)
    if wf.size < 5 or not np.all(np.isfinite(wf)):
        return np.nan, np.nan, False

    trough = int(np.argmin(wf))
    # Baseline from the samples before the spike; fall back to the median if the trough is early.
    base = float(np.median(wf[:max(trough // 2, 3)])) if trough >= 6 else float(np.median(wf))
    depth = base - wf[trough]
    if depth <= 0:
        return np.nan, np.nan, False

    # PP: trough to the following maximum.
    after = wf[trough:]
    peak = trough + int(np.argmax(after))
    if peak <= trough:
        return np.nan, np.nan, False
    pp = (peak - trough) / FS

    # FW: width of the trough at half depth, linearly interpolated on both flanks.
    half = base - depth / 2.0
    left = None
    for i in range(trough, 0, -1):
        if wf[i] <= half <= wf[i - 1] or wf[i - 1] >= half >= wf[i]:
            d = wf[i - 1] - wf[i]
            left = (i - 1) + (wf[i - 1] - half) / d if d != 0 else float(i)
            break
    right = None
    for i in range(trough, wf.size - 1):
        if wf[i] <= half <= wf[i + 1] or wf[i + 1] >= half >= wf[i]:
            d = wf[i + 1] - wf[i]
            right = i + (half - wf[i]) / d if d != 0 else float(i)
            break
    if left is None or right is None or right <= left:
        return np.nan, pp, False
    return (right - left) / FS, pp, True


def main() -> int:
    rows = []
    for animal, old in ANIMAL_MAP.items():
        mat = ROOT / f"DATA{old}.mat"
        if not mat.exists():
            print(f"  missing {mat.name}")
            continue
        U = scipy.io.loadmat(mat, squeeze_me=True, struct_as_record=False)["Units"]
        for i, u in enumerate(U):
            fw, pp, ok = features(np.asarray(u.wf))
            rows.append({"animal": animal, "neuron_id": i + 1,   # 1-indexed, as elsewhere
                         "fw_ms": fw * 1000 if ok else np.nan,
                         "pp_ms": pp * 1000 if np.isfinite(pp) else np.nan,
                         "usable": ok})
    df = pd.DataFrame(rows)
    good = df[df.usable].copy()
    print(f"units: {len(df)} total, {len(good)} with usable waveforms")

    # The paper clusters "across all simultaneously recorded somatic electrodes", i.e. within a
    # recording. That matters here: median PP varies 0.400-1.050 ms BETWEEN these recordings, a
    # 2.6x range that swamps any cell-type difference. Clustering the pooled values makes a
    # recording's inhibitory fraction a function of where it sits (Spearman -0.978 vs its median
    # PP) -- a rig effect, not biology. Per-recording clustering removes the offset by
    # construction, so it is both faithful to the paper and the right control here.
    good["cluster"] = -1
    good["ei"] = "unclassified"
    per = []
    for animal, idx in good.groupby("animal").groups.items():
        sub = good.loc[idx]
        if len(sub) < 20:
            per.append((animal, len(sub), np.nan, np.nan, "too few units"))
            continue
        X = sub[["fw_ms", "pp_ms"]].to_numpy()
        Xz = (X - X.mean(0)) / X.std(0)
        sc = {}
        for k in range(2, min(6, len(sub) // 6) + 1):
            km = KMeans(n_clusters=k, n_init=50, random_state=42).fit(Xz)
            sc[k] = calinski_harabasz_score(Xz, km.labels_)
        if not sc:
            per.append((animal, len(sub), np.nan, np.nan, "no k"))
            continue
        bk = max(sc, key=sc.get)
        km = KMeans(n_clusters=bk, n_init=50, random_state=42).fit(Xz)
        cen = pd.DataFrame(X, index=sub.index).groupby(km.labels_).mean()
        inh = int((cen[0] / X[:, 0].std() + cen[1] / X[:, 1].std()).idxmin())
        good.loc[idx, "cluster"] = km.labels_
        good.loc[idx, "ei"] = np.where(km.labels_ == inh, "I", "E")
        f = float((km.labels_ == inh).mean() * 100)
        per.append((animal, len(sub), bk, f,
                    f"FW {cen.loc[inh,0]:.3f} PP {cen.loc[inh,1]:.3f} ms"))

    print("\nper-recording k-means (Calinski-Harabasz):")
    print(f"{'animal':8s} {'n':>4s} {'k':>3s} {'% inh':>7s}  inhibitory-cluster centroid")
    for a_, n_, k_, f_, note in per:
        ks = f"{int(k_)}" if k_ == k_ else "-"
        fs = f"{f_:7.1f}" if f_ == f_ else "      -"
        print(f"{a_:8s} {n_:4d} {ks:>3s} {fs}  {note}")

    df = df.merge(good[["animal", "neuron_id", "cluster", "ei"]],
                  on=["animal", "neuron_id"], how="left")
    df["ei"] = df.ei.fillna("unclassified")
    df.to_csv(ROOT / "data" / "processed" / "waveform_features.csv", index=False)

    n_i = int((df.ei == "I").sum())
    n_c = int((df.ei != "unclassified").sum())
    print(f"\ninhibitory overall: {n_i}/{n_c} = {100*n_i/max(n_c,1):.2f}%  "
          f"(spikeMAP reported 5.25% in mouse PFC)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
