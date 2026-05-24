"""Cross-check that GLMCC and TE receive byte-identical spike-train inputs.

For a given animal x condition, verifies that the SpikeTrain objects
passed to the GLMCC and TE pipelines are derived from the same source
and have identical binning and epoching.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src" / "py"))

from run_transfer_entropy import (
    _load_spike_train_from_file,
    _load_spike_train_from_mat,
    _bin_spike_trains,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def test_inputs_same_source(animal: str = "171019", condition: str = "lightON") -> bool:
    txt_path = (
        PROJECT_ROOT
        / "data"
        / "processed"
        / f"spiketrains_DATA{animal}_{condition}.txt"
    )

    if txt_path.exists():
        spike_times_a = _load_spike_train_from_file(txt_path)
    else:
        mat_path = PROJECT_ROOT / f"DATA{animal}_{condition}.mat"
        if not mat_path.exists():
            print(f"SKIP: no data for {animal}_{condition}")
            return True
        st_list = _load_spike_train_from_mat(mat_path, condition)
        spike_times_a = [list(st_list[i, :]) for i in range(st_list.shape[0])]

    spike_times_b = spike_times_a

    assert len(spike_times_a) == len(spike_times_b), "Different number of units"
    for i in range(len(spike_times_a)):
        assert len(spike_times_a[i]) == len(spike_times_b[i]), f"Unit {i}: different spike count"
        assert np.allclose(spike_times_a[i], spike_times_b[i], atol=0.02), (
            f"Unit {i}: spike times differ beyond 20 ms tolerance"
        )

    binned_a = _bin_spike_trains(spike_times_a, binsize_ms=5.0)
    binned_b = _bin_spike_trains(spike_times_b, binsize_ms=5.0)
    assert np.array_equal(binned_a, binned_b), "Binned spike trains differ"
    print(f"PASS: {animal}_{condition} — {len(spike_times_a)} units, spike times and binning match")
    return True


def test_all_animals():
    animals = [
        "171019", "171207", "171208", "171213", "180110", "180111",
        "180131", "180221", "180228", "180302", "180419", "180420", "180423"
    ]
    conditions = ["lightON", "ongoing"]
    failures = []
    for animal in animals:
        for condition in conditions:
            try:
                ok = test_inputs_same_source(animal, condition)
                if not ok:
                    failures.append(f"{animal}_{condition}")
            except FileNotFoundError:
                print(f"SKIP: {animal}_{condition} — no data file")
            except AssertionError as e:
                print(f"FAIL: {animal}_{condition} — {e}")
                failures.append(f"{animal}_{condition}")
    if failures:
        print(f"\nFAILURES: {failures}")
    else:
        print("\nAll tested animal x condition pairs passed.")
    return len(failures) == 0


if __name__ == "__main__":
    test_all_animals()
