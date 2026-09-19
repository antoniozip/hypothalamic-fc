#!/usr/bin/env python3
"""Differential analysis by putative excitatory / inhibitory class.

Classification
--------------
`{animal}_put_inh.txt` lists putative inhibitory units, produced in `tim_stuff.m:714` by a
narrow-spiking waveform criterion:

    putative_inh = intersect(find(half_w < 0.0006), find(trough_to_peak < 0.00065))

i.e. half-peak width < 0.6 ms AND trough-to-peak < 0.65 ms. MATLAB `find()` is 1-indexed, which
matches `neuron_id` in `neurons.csv` and in the metrics files. The criterion is derived from
spike waveforms and is therefore **independent of GLMCC**, so using it to interrogate the
connectivity is not circular.

**Only 7 of 13 animals are classifiable.** day2/3/4/5/6 and night6/7 label 38.7-51.1% of their
units inhibitory; day1, night1-5 label 0-3.2%. That gap (3.2% -> 38.7%) is not a biological
range, it is the difference between the classification having been run and not. Units in the
six unrun animals are marked `unclassified`, **not** excitatory -- the existing
`src/py/ei_classification.py` defaults them to "E", which silently turns missing data into
221 false excitatory labels.

The confound
------------
Fast-spiking interneurons fire faster, and GLMCC's detection threshold falls with spike count
(`J_min ~ 1/sqrt(cc_0)`, `cc_0 ~ n_i*n_j/T`). Any E vs I difference in connectivity could
therefore be a firing-rate difference in disguise. Spike counts are reported per class, and the
downstream LME carries log10(spike count) as a covariate.

Analyses
--------
1. firing rate by class -- sizes the confound
2. validation: does the sign of a unit's outgoing couplings track its waveform class?
3. node metrics by class
4. edge-type composition: E->E, E->I, I->E, I->I

Matrices: `results/glmcc_rm/` (the rate- and geometry-matched arms), because the raw contrast is
a detection-power artifact (`reports/rate_matched/`).

Output: data/processed/neurons_ei_strict.csv
        results/ei/{rates,outgoing_sign,node_metrics,edge_types}.csv
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
OUT = ROOT / "results" / "ei"

ANIMAL_MAP = {
    "day1": "171019", "day2": "180131", "day3": "180302", "day4": "180419",
    "day5": "180420", "day6": "180423", "night1": "171207", "night2": "171208",
    "night3": "171213", "night4": "180110", "night5": "180111",
    "night6": "180221", "night7": "180228",
}
MIN_INH_PCT = 10.0   # below this the classification was not run; see module docstring


def load_put_inh(old: str) -> set[int]:
    p = ROOT / f"{old}_put_inh.txt"
    if not p.exists():
        return set()
    return {int(x) for x in re.split(r"[,\s]+", p.read_text().strip()) if x.strip().isdigit()}


def build_table() -> pd.DataFrame:
    """neurons.csv + ei in {E, I, unclassified} + classifiable flag. neuron_id is 1-indexed."""
    n = pd.read_csv(ROOT / "data" / "processed" / "neurons.csv")
    ei, classifiable = [], []
    frac = {}
    for a, old in ANIMAL_MAP.items():
        inh = load_put_inh(old)
        units = int((n.animal == a).sum())
        frac[a] = 100.0 * len(inh) / max(units, 1)
    for _, r in n.iterrows():
        a = str(r["animal"])
        ok = frac.get(a, 0.0) >= MIN_INH_PCT
        classifiable.append(ok)
        if not ok:
            ei.append("unclassified")
        else:
            ei.append("I" if int(r["neuron_id"]) in load_put_inh(ANIMAL_MAP[a]) else "E")
    n["ei"] = ei
    n["classifiable"] = classifiable
    n["pct_inh_animal"] = n.animal.map(frac).round(1)
    return n


def spike_counts(animal: str, condition: str) -> np.ndarray:
    """Spikes per unit in the arm actually analysed, indexed 0..n-1 (neuron_id - 1)."""
    old = ANIMAL_MAP[animal]
    d = ROOT / (f"DATA{old}" if condition == "lightON" else f"DATA{old}_ongoing")
    cells = sorted(d.glob("cell*.txt"), key=lambda p: int(p.stem[4:]))
    return np.array([np.atleast_1d(np.loadtxt(c)).size for c in cells])


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    tab = build_table()
    tab.to_csv(ROOT / "data" / "processed" / "neurons_ei_strict.csv", index=False)
    cls = tab[tab.classifiable]
    animals = sorted(cls.animal.unique())
    print(f"classifiable: {len(animals)} animals, {len(cls)} units "
          f"({(cls.ei=='I').sum()} I / {(cls.ei=='E').sum()} E)\n")

    rate_rows, sign_rows, node_rows, edge_rows = [], [], [], []

    for a in animals:
        sub = cls[cls.animal == a].sort_values("neuron_id")
        lab = sub.set_index("neuron_id").ei.to_dict()       # 1-indexed
        for cond in ("lightON", "ongoing"):
            adj_p = ROOT / "results" / "glmcc_rm" / f"validated_adj_{a}_{cond}.csv"
            met_p = ROOT / "results" / "glmcc_rm" / f"metrics_{a}_{cond}.csv"
            if not adj_p.exists():
                continue
            W = np.loadtxt(adj_p, delimiter=",")
            n = W.shape[0]
            counts = spike_counts(a, cond)
            if counts.size != n:
                print(f"  skip {a}/{cond}: {counts.size} cell files vs {n}x{n} matrix")
                continue
            klass = np.array([lab.get(i + 1, "unclassified") for i in range(n)])

            # 1. firing rate by class
            for k in ("E", "I"):
                m = klass == k
                if m.any():
                    rate_rows.append({"animal": a, "condition": cond, "ei": k,
                                      "n_units": int(m.sum()),
                                      "median_spikes": float(np.median(counts[m])),
                                      "mean_spikes": float(counts[m].mean())})

            # 2. validation: sign of a unit's OUTGOING couplings vs its waveform class
            for i in range(n):
                if klass[i] == "unclassified":
                    continue
                out = W[i, :][np.arange(n) != i]
                nz = out[out != 0]
                if nz.size:
                    sign_rows.append({"animal": a, "condition": cond, "neuron_id": i + 1,
                                      "ei": klass[i], "n_out": int(nz.size),
                                      "frac_neg": float((nz < 0).mean()),
                                      "median_w": float(np.median(nz)),
                                      "spikes": int(counts[i])})

            # 3. node metrics by class
            if met_p.exists():
                mt = pd.read_csv(met_p)
                mt["ei"] = mt.neuron_id.map(lab)
                mt = mt[mt.ei.isin(["E", "I"])].copy()
                mt["spikes"] = mt.neuron_id.map(lambda q: int(counts[q - 1]))
                node_rows.append(mt)

            # 4. edge-type composition
            for pre in ("E", "I"):
                for post in ("E", "I"):
                    pm, qm = klass == pre, klass == post
                    if not (pm.any() and qm.any()):
                        continue
                    blk = W[np.ix_(pm, qm)]
                    if pre == post:               # drop the diagonal of a same-class block
                        blk = blk[~np.eye(blk.shape[0], dtype=bool)]
                    nz = blk[blk != 0]
                    poss = blk.size
                    edge_rows.append({
                        "animal": a, "condition": cond, "pre": pre, "post": post,
                        "possible": int(poss), "edges": int(nz.size),
                        "density_pct": round(100 * nz.size / max(poss, 1), 3),
                        "frac_neg": round(float((nz < 0).mean()), 4) if nz.size else np.nan,
                        "median_w": round(float(np.median(nz)), 4) if nz.size else np.nan,
                    })

    pd.DataFrame(rate_rows).to_csv(OUT / "rates.csv", index=False)
    pd.DataFrame(sign_rows).to_csv(OUT / "outgoing_sign.csv", index=False)
    pd.concat(node_rows).to_csv(OUT / "node_metrics.csv", index=False)
    pd.DataFrame(edge_rows).to_csv(OUT / "edge_types.csv", index=False)
    print(f"wrote {OUT}/{{rates,outgoing_sign,node_metrics,edge_types}}.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
