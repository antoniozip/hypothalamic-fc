#!/usr/bin/env python3
"""Re-validate every animal x condition through one validator.

Replaces the split pipeline in which lightON came from validate_ccg.py (.mat,
milliseconds, +/-50 ms window, 100 surrogates, weighted output) while ongoing came
from the untracked scripts/fast_validate_ongoing.py (cell*.txt, seconds, so a
+/-50 s window, 50 surrogates, binary output). Every pair now uses the same
source, units, window, surrogate count and encoding.

Output: results/glmcc/validated_adj_{animal}_{condition}.csv holding GLMCC
coupling weights masked by the CCG shuffle-ISI test, plus a run summary at
results/revalidation_summary.csv.

Usage: python scripts/revalidate_all.py [--n-surrogates 100] [--conditions lightON ongoing]
"""

import argparse
import csv
import sys
import time
import traceback
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src" / "py"))

from validate_edges import validate_edges, RankMismatchError  # noqa: E402

ANIMALS = [
    "day1", "day2", "day3", "day4", "day5", "day6",
    "night1", "night2", "night3", "night4", "night5", "night6", "night7",
]

SUMMARY_FIELDS = [
    "animal", "condition", "status", "n_units", "n_possible_edges",
    "n_significant", "edge_survival_rate", "n_surrogates", "encoding",
    "elapsed_s", "note",
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-surrogates", type=int, default=100)
    parser.add_argument("--fdr-q", type=float, default=0.05)
    parser.add_argument("--conditions", nargs="+", default=["lightON", "ongoing"])
    parser.add_argument("--animals", nargs="+", default=ANIMALS)
    args = parser.parse_args()

    rows = []
    for condition in args.conditions:
        for animal in args.animals:
            started = time.time()
            label = f"{animal}_{condition}"
            print(f"\n{'=' * 70}\n{label}\n{'=' * 70}", flush=True)
            row = {
                "animal": animal, "condition": condition, "n_surrogates": args.n_surrogates,
                "encoding": "glmcc_weight_masked", "note": "",
            }
            try:
                _, stats = validate_edges(
                    animal, condition, estimator="glmcc", method="ccg",
                    fdr_q=args.fdr_q, n_surrogates=args.n_surrogates,
                )
                row.update({
                    "status": "ok",
                    "n_units": stats["n_units"],
                    "n_possible_edges": stats["n_possible_edges"],
                    "n_significant": stats["n_significant"],
                    "edge_survival_rate": stats["edge_survival_rate"],
                })
            except RankMismatchError as exc:
                # A matrix that does not describe this unit set must not be
                # masked and shipped; it would silently mislabel every region.
                row.update({"status": "skipped", "note": str(exc)})
                print(f"SKIPPED: {exc}", flush=True)
            except FileNotFoundError as exc:
                row.update({"status": "missing", "note": str(exc)})
                print(f"MISSING: {exc}", flush=True)
            except Exception as exc:  # noqa: BLE001 - one bad pair must not stop the run
                row.update({"status": "error", "note": f"{type(exc).__name__}: {exc}"})
                traceback.print_exc()

            row["elapsed_s"] = round(time.time() - started, 1)
            rows.append(row)
            print(f"[{label}] {row['status']} in {row['elapsed_s']}s", flush=True)

    tag = "_".join(args.conditions)
    out_path = PROJECT_ROOT / "results" / f"revalidation_summary_{tag}.csv"
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in SUMMARY_FIELDS})

    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"\n{'=' * 70}\n{ok}/{len(rows)} pairs validated. Summary: {out_path}")
    for row in rows:
        if row["status"] != "ok":
            print(f"  {row['status'].upper():8s} {row['animal']}_{row['condition']}: {row['note']}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
