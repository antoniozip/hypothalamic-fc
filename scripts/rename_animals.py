#!/usr/bin/env python3
"""Rename animals from date-based IDs to dayN/nightN labels.

Mapping:
  Day:   171019→day1, 180131→day2, 180302→day3, 180419→day4, 180420→day5, 180423→day6
  Night: 171207→night1, 171208→night2, 171213→night3, 180110→night4, 180111→night5,
         180221→night6, 180228→night7

Actions:
  1. Rename all CSV/PNG files containing animal IDs
  2. Update animal column values inside all CSV files
  3. Update neurons.csv and neurons_ei.csv
  4. Update hardcoded animal lists in scripts (*.R, *.py)
"""

import os
import re
from pathlib import Path
import pandas as pd
import shutil

PROJECT_ROOT = Path("/media/antonio/data/tim_brown")

# ─── Mapping ──────────────────────────────────────────────────────
OLD_IDS = [
    "171019", "171207", "171208", "171213", "180110", "180111",
    "180131", "180221", "180228", "180302", "180419", "180420", "180423",
]

DAYNIGHT = {
    "171019": "D", "171207": "N", "171208": "N", "171213": "N",
    "180110": "N", "180111": "N", "180131": "D", "180221": "N",
    "180228": "N", "180302": "D", "180419": "D", "180420": "D",
    "180423": "D",
}

# Build mapping: old → new
day_counter = 0
night_counter = 0
ANIMAL_MAP = {}
for old_id in OLD_IDS:
    if DAYNIGHT[old_id] == "D":
        day_counter += 1
        ANIMAL_MAP[old_id] = f"day{day_counter}"
    else:
        night_counter += 1
        ANIMAL_MAP[old_id] = f"night{night_counter}"

# Reverse map for sorting
NEW_ORDER_DAY = [v for k, v in ANIMAL_MAP.items() if DAYNIGHT[k] == "D"]
NEW_ORDER_NIGHT = [v for k, v in ANIMAL_MAP.items() if DAYNIGHT[k] == "N"]
NEW_ORDER = NEW_ORDER_DAY + NEW_ORDER_NIGHT

print("Animal mapping:")
for old, new in ANIMAL_MAP.items():
    print(f"  {old} → {new} ({DAYNIGHT[old]})")


def rename_path_component(path_str: str) -> str:
    """Replace all old IDs in a path string with new IDs."""
    result = path_str
    # Sort by length descending to avoid partial matches (e.g., 180131 before 18013)
    for old_id in sorted(ANIMAL_MAP.keys(), key=len, reverse=True):
        if old_id in result:
            result = result.replace(old_id, ANIMAL_MAP[old_id])
    return result


def update_animal_column(csv_path: Path):
    """Update animal column values in a CSV file."""
    try:
        df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    except Exception:
        return False

    modified = False
    for col in df.columns:
        if "animal" in col.lower():
            for old_id, new_id in ANIMAL_MAP.items():
                mask = df[col].str.strip().str.fullmatch(old_id, na=False)
                if mask.any():
                    df.loc[mask, col] = new_id
                    modified = True

    if modified:
        df.to_csv(csv_path, index=False)
        return True
    return False


# ═══════════════════════════════════════════════════════════════════
# STEP 1: Rename CSV files
# ═══════════════════════════════════════════════════════════════════
print("\n=== Step 1: Renaming CSV files ===")

csv_renames = []
for root, dirs, files in os.walk(PROJECT_ROOT / "results"):
    for fname in files:
        if not fname.endswith(".csv"):
            continue
        old_path = Path(root) / fname
        new_name = rename_path_component(fname)
        if new_name != fname:
            new_path = Path(root) / new_name
            csv_renames.append((old_path, new_path))

# Also data/processed/
for root, dirs, files in os.walk(PROJECT_ROOT / "data" / "processed"):
    for fname in files:
        if not fname.endswith(".csv"):
            continue
        old_path = Path(root) / fname
        new_name = rename_path_component(fname)
        if new_name != fname:
            csv_renames.append((old_path, new_path))

for old_p, new_p in csv_renames:
    if new_p.exists():
        print(f"  SKIP (exists): {new_p.name}")
        continue
    shutil.move(str(old_p), str(new_p))
    print(f"  {old_p.name} → {new_p.name}")

print(f"  Renamed {len(csv_renames)} files")


# ═══════════════════════════════════════════════════════════════════
# STEP 2: Update animal columns inside all CSV files
# ═══════════════════════════════════════════════════════════════════
print("\n=== Step 2: Updating animal columns in CSVs ===")

updated_count = 0
for root, dirs, files in os.walk(PROJECT_ROOT / "results"):
    for fname in files:
        if not fname.endswith(".csv"):
            continue
        if update_animal_column(Path(root) / fname):
            updated_count += 1

# Also neurons.csv and neurons_ei.csv
for fname in ["neurons.csv", "neurons_ei.csv"]:
    p = PROJECT_ROOT / "data" / "processed" / fname
    if p.exists() and update_animal_column(p):
        updated_count += 1
        print(f"  Updated: {fname}")

print(f"  Updated {updated_count} CSV files")


# ═══════════════════════════════════════════════════════════════════
# STEP 3: Rename figure files
# ═══════════════════════════════════════════════════════════════════
print("\n=== Step 3: Renaming figure files ===")

png_renames = []
for root, dirs, files in os.walk(PROJECT_ROOT / "figures"):
    for fname in files:
        if not fname.endswith(".png"):
            continue
        old_path = Path(root) / fname
        new_name = rename_path_component(fname)
        if new_name != fname:
            new_path = Path(root) / new_name
            png_renames.append((old_path, new_path))

for old_p, new_p in png_renames:
    if new_p.exists():
        print(f"  SKIP (exists): {new_p.name}")
        continue
    shutil.move(str(old_p), str(new_p))
    print(f"  {old_p.name} → {new_p.name}")

print(f"  Renamed {len(png_renames)} files")


# ═══════════════════════════════════════════════════════════════════
# STEP 4: Update hardcoded animal lists in scripts
# ═══════════════════════════════════════════════════════════════════
print("\n=== Step 4: Updating script animal lists ===")

# Find scripts with hardcoded animal lists
script_dirs = [
    PROJECT_ROOT / "src" / "py",
    PROJECT_ROOT / "src" / "r",
]

for sdir in script_dirs:
    if not sdir.exists():
        continue
    for fpath in sdir.iterdir():
        if fpath.suffix not in (".py", ".R", ".sh"):
            continue
        content = fpath.read_text()
        modified = False
        for old_id, new_id in ANIMAL_MAP.items():
            if old_id in content:
                content = content.replace(f'"{old_id}"', f'"{new_id}"')
                content = content.replace(f"'{old_id}'", f"'{new_id}'")
                # Also replace in unquoted contexts (e.g., c(171019, ...))
                # But be careful not to replace partial matches in numbers
                # Use word boundary
                content = re.sub(rf'\b{old_id}\b', new_id, content)
                modified = True
        if modified:
            fpath.write_text(content)
            print(f"  Updated: {fpath.name}")

# Also update Snakefile
snakefile = PROJECT_ROOT / "Snakefile"
if snakefile.exists():
    content = snakefile.read_text()
    modified = False
    for old_id, new_id in ANIMAL_MAP.items():
        if old_id in content:
            content = content.replace(f'"{old_id}"', f'"{new_id}"')
            modified = True
    if modified:
        snakefile.write_text(content)
        print(f"  Updated: Snakefile")

# Also update the daynight mapping in scripts
daynight_old_str = str(DAYNIGHT)
# Build new daynight mapping
NEW_DAYNIGHT = {ANIMAL_MAP[k]: v for k, v in DAYNIGHT.items()}

print(f"\nNew day/night mapping:")
for k, v in NEW_DAYNIGHT.items():
    print(f"  {k}: {v}")

# Save mapping for reference
mapping_df = pd.DataFrame([
    {"old_id": k, "new_id": v, "day_night": DAYNIGHT[k], "neurons": ""}
    for k, v in ANIMAL_MAP.items()
])
mapping_df.to_csv(PROJECT_ROOT / "data" / "processed" / "animal_mapping.csv", index=False)
print(f"\n  Saved animal_mapping.csv")

print("\n=== Rename complete ===")
print("NOTE: You should re-run key scripts (make_figures.R, tier1_analyses.R, etc.)")
print("to regenerate figures with the new labels embedded in titles/axes.")
