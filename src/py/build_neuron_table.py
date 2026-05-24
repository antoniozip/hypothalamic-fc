#!/usr/bin/env python3
"""Build per-neuron region mapping by matching .mat file electrode coords to region xlsx.

Output: data/processed/neurons.csv with columns:
  animal, electrode_id, neuron_id, region6, region9, ap, ml, dv
"""

from pathlib import Path
import csv

import numpy as np
from scipy.io import loadmat
import openpyxl


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


def _parse_unit_list(units_str) -> list[int]:
    """Parse unit entry which may be comma-separated string or single int."""
    if units_str is None:
        return []
    if isinstance(units_str, (int, float)):
        return [int(units_str)]
    if isinstance(units_str, str):
        return [int(x.strip()) for x in units_str.split(",") if x.strip()]
    return []


def build_neuron_table() -> None:
    wb = openpyxl.load_workbook(
        str(PROJECT_ROOT / "regions_match_revision_March2024.xlsx")
    )
    ws = wb["Sheet1"]

    electrode_rows = []
    for row in ws.iter_rows(min_row=2, max_row=ws.max_row, values_only=True):
        rec_id = row[0]
        if rec_id is None or not str(rec_id).isdigit():
            continue
        units_str = row[1]
        if units_str is None:
            continue
        ap = float(row[2]) if row[2] is not None else None
        ml = float(row[3]) if row[3] is not None else None
        dv = float(row[4]) if row[4] is not None else None
        region = row[5] if row[5] else None
        macro_region = row[8] if row[8] else None
        further = row[9] if row[9] else None
        if region and macro_region:
            # Normalize names
            if macro_region.lower().endswith("thalamus"):
                region6 = "VM-thalamus"
            elif macro_region.lower().startswith("paraventricular"):
                region6 = "PVH"
            elif macro_region.lower().startswith("dorsomedial"):
                region6 = "DMH"
            elif macro_region.lower().startswith("arcuate"):
                region6 = "ARH"
            elif macro_region.lower().startswith("posterior"):
                region6 = "PH"
            elif macro_region.lower().startswith("zona"):
                region6 = "ZI"
            elif macro_region.lower().startswith("ventromedial"):
                region6 = "VMH"
            elif macro_region.lower().startswith("anterior"):
                region6 = "AH"
            elif macro_region.lower().endswith("complex"):
                region6 = "MC"
            else:
                region6 = macro_region
        else:
            region6 = "unknown"
        electrode_rows.append({
            "animal": str(rec_id),
            "units": units_str,
            "ap": ap, "ml": ml, "dv": dv,
            "region": region,
            "region6": region6,
            "region9": further if further else region6,
        })

    neurons = []
    for rec in electrode_rows:
        animal = rec["animal"]
        unit_ids = _parse_unit_list(rec["units"])
        for uid in unit_ids:
            neurons.append({
                "animal": animal,
                "electrode_id": uid,
                "neuron_id": uid,
                "region6": rec["region6"],
                "region9": rec["region9"],
                "ap": rec["ap"],
                "ml": rec["ml"],
                "dv": rec["dv"],
            })

    with open(PROJECT_ROOT / "data" / "processed" / "neurons.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "animal", "electrode_id", "neuron_id", "region6", "region9", "ap", "ml", "dv"
        ])
        writer.writeheader()
        writer.writerows(neurons)
    print(f"Wrote {len(neurons)} neuron rows to data/processed/neurons.csv")


if __name__ == "__main__":
    build_neuron_table()
