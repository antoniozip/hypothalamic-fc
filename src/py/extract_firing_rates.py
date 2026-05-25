#!/usr/bin/env python3
"""Extract firing rates from .mat files and generate comparison plots.

Reads the original .mat spike time data, computes per-neuron firing rates
for ongoing and lightON conditions, and generates a summary figure.

Output: figures/main/figS5_firing_rates.png
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
FIGURES_DIR = PROJECT_ROOT / "figures" / "main"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

ANIMALS = [
    "day1", "night1", "night2", "night3", "night4", "night5",
    "day2", "night6", "night7", "day3", "day4", "day5", "day6",
]
DAYNIGHT = {
    "day1": "D", "night1": "N", "night2": "N", "night3": "N",
    "night4": "N", "night5": "N", "day2": "D", "night6": "N",
    "night7": "N", "day3": "D", "day4": "D", "day5": "D",
    "day6": "D",
}

TS_KEY_MAP = {
    "lightON": "evoked_ts",
    "ongoing": "ongoing_ts",
    "evoked": "evoked_ts",
}


def load_spike_times(animal: str, condition: str) -> list[np.ndarray]:
    """Load spike times (seconds) from .mat file."""
    mat_path = PROJECT_ROOT / f"DATA{animal}_{condition}.mat"
    if not mat_path.exists():
        if condition == "ongoing":
            alts = sorted(PROJECT_ROOT.glob(f"DATA{animal}*ongoing*.mat"))
            if alts:
                mat_path = alts[0]
    if not mat_path.exists():
        return None

    m = loadmat(str(mat_path), simplify_cells=False)
    ts_key = TS_KEY_MAP.get(condition, "evoked_ts")
    if ts_key not in m:
        ts_keys = [k for k in m if k.endswith("_ts")]
        if not ts_keys:
            return None
        ts_key = ts_keys[0]

    ts = m[ts_key]
    n_units = min(ts.shape[1], 64)  # max 64 channels
    times_sec = []
    for i in range(n_units):
        try:
            t = np.atleast_1d(ts[0, i][0]).astype(float)
            # .mat values are in ms — convert to seconds
            if len(t) > 0 and t.max() > 100:
                t = t / 1000.0
            times_sec.append(t)
        except (IndexError, ValueError):
            times_sec.append(np.array([]))
    return times_sec


def compute_firing_rate(spike_times: np.ndarray, duration_s: float | None = None) -> float:
    """Compute firing rate in Hz."""
    if len(spike_times) < 1:
        return 0.0
    if duration_s is None:
        duration_s = spike_times[-1] - spike_times[0]
        if duration_s <= 0:
            duration_s = 1.0
    return len(spike_times) / duration_s


def main():
    # Collect firing rates
    fr_data = []  # list of dicts

    for animal in ANIMALS:
        times_l = load_spike_times(animal, "lightON")
        times_o = load_spike_times(animal, "ongoing")

        if times_l is None:
            print(f"  SKIP {animal} — no lightON .mat file")
            continue

        n_neurons = len(times_l)
        fr_light = [compute_firing_rate(t) for t in times_l]
        fr_ongoing = [compute_firing_rate(t) for t in times_o] if times_o else None

        for i in range(n_neurons):
            entry = {
                "animal": animal,
                "day_night": DAYNIGHT.get(animal, "?"),
                "neuron_id": i + 1,
                "fr_lightON": fr_light[i],
                "fr_ongoing": fr_ongoing[i] if fr_ongoing and i < len(fr_ongoing) else None,
            }
            # Compute change if both available
            if entry["fr_ongoing"] is not None and entry["fr_ongoing"] > 0:
                entry["fr_change_pct"] = (
                    (entry["fr_lightON"] - entry["fr_ongoing"]) / entry["fr_ongoing"] * 100
                )
            else:
                entry["fr_change_pct"] = None
            fr_data.append(entry)

    if not fr_data:
        print("No data loaded.")
        return

    # Convert to structured arrays for plotting
    animals_list = [d["animal"] for d in fr_data]
    fr_l = np.array([d["fr_lightON"] for d in fr_data])
    fr_o = np.array([d["fr_ongoing"] for d in fr_data if d["fr_ongoing"] is not None])
    dn = np.array([d["day_night"] for d in fr_data])
    fr_change = np.array([d["fr_change_pct"] for d in fr_data if d["fr_change_pct"] is not None])

    print(f"  Loaded {len(fr_data)} neurons from {len(set(animals_list))} animals")
    print(f"  LightON FR: {np.mean(fr_l):.3f} ± {np.std(fr_l):.3f} Hz")
    if len(fr_o) > 0:
        print(f"  Ongoing FR: {np.mean(fr_o):.3f} ± {np.std(fr_o):.3f} Hz")
    if len(fr_change) > 0:
        print(f"  Mean change: {np.mean(fr_change):.1f}%")

    # ─── Figure ───────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Firing Rate Analysis — LightON vs Ongoing",
                 fontsize=15, fontweight="bold")

    # Panel A: Per-animal mean FR (bar chart)
    ax = axes[0, 0]
    animal_order = sorted(set(animals_list))
    fr_light_by_animal = []
    fr_ongoing_by_animal = []
    for a in animal_order:
        mask = np.array(animals_list) == a
        fr_light_by_animal.append(np.mean(fr_l[mask]))
        ongoing_mask = [i for i, d in enumerate(fr_data)
                        if d["animal"] == a and d["fr_ongoing"] is not None]
        if ongoing_mask:
            fr_ongoing_by_animal.append(np.mean(
                [fr_data[i]["fr_ongoing"] for i in ongoing_mask]))
        else:
            fr_ongoing_by_animal.append(0)

    x = np.arange(len(animal_order))
    w = 0.35
    colors = ["#2166ac" if DAYNIGHT.get(a) == "D" else "#b2182b" for a in animal_order]
    ax.bar(x - w/2, fr_light_by_animal, w, label="LightON",
           color="#f4a582", alpha=0.85)
    ax.bar(x + w/2, fr_ongoing_by_animal, w, label="Ongoing",
           color="#92c5de", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(animal_order, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mean Firing Rate (Hz)")
    ax.set_title("A: Per-animal mean firing rate")
    ax.legend(fontsize=8)

    # Panel B: Scatter — LightON vs Ongoing per neuron
    ax = axes[0, 1]
    paired = [(d["fr_ongoing"], d["fr_lightON"]) for d in fr_data
              if d["fr_ongoing"] is not None and d["fr_lightON"] > 0]
    if paired:
        x_vals, y_vals = zip(*paired)
        ax.scatter(x_vals, y_vals, alpha=0.3, s=8, c="#444444")
        max_val = max(max(x_vals), max(y_vals)) * 1.1
        ax.plot([0, max_val], [0, max_val], "k--", alpha=0.3, linewidth=0.8)
        ax.set_xlim(0, max_val)
        ax.set_ylim(0, max_val)
    ax.set_xlabel("Ongoing FR (Hz)")
    ax.set_ylabel("LightON FR (Hz)")
    ax.set_title(f"B: Per-neuron FR (n={len(paired)} pairs)")

    # Panel C: FR change histogram
    ax = axes[1, 0]
    if len(fr_change) > 0:
        ax.hist(fr_change, bins=40, color="#6baed6", edgecolor="white", alpha=0.85)
        ax.axvline(0, color="red", linestyle="--", linewidth=1)
        ax.axvline(np.mean(fr_change), color="#08519c", linestyle="-", linewidth=1.5,
                   label=f"Mean = {np.mean(fr_change):.1f}%")
        ax.legend(fontsize=9)
    ax.set_xlabel("FR Change (LightON − Ongoing) %")
    ax.set_ylabel("Neuron Count")
    ax.set_title(f"C: Firing rate change distribution (n={len(fr_change)} neurons)")

    # Panel D: FR change by animal
    ax = axes[1, 1]
    fr_change_by_animal = {}
    for d in fr_data:
        if d["fr_change_pct"] is not None:
            a = d["animal"]
            if a not in fr_change_by_animal:
                fr_change_by_animal[a] = []
            fr_change_by_animal[a].append(d["fr_change_pct"])

    animals_with_data = sorted(fr_change_by_animal.keys())
    means = [np.mean(fr_change_by_animal[a]) for a in animals_with_data]
    stds = [np.std(fr_change_by_animal[a]) for a in animals_with_data]
    colors_dn = ["#2166ac" if DAYNIGHT.get(a) == "D" else "#b2182b"
                 for a in animals_with_data]
    xs = np.arange(len(animals_with_data))
    ax.bar(xs, means, color=colors_dn, alpha=0.85, yerr=stds, capsize=3)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.7)
    ax.set_xticks(xs)
    ax.set_xticklabels(animals_with_data, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mean FR Change (%)")
    ax.set_title("D: FR change by animal (mean ± SD)")

    plt.tight_layout()
    out_path = FIGURES_DIR / "figS5_firing_rates.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
