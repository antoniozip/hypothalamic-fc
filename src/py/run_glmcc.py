#!/usr/bin/env python3
"""Wrapper around GLMCC (Kobayashi et al. 2019) for connectivity inference.

Reads spike trains from data/processed/spiketrains_{animal}_{condition}.txt,
calls the GLMCC estimator, and writes the adjacency matrix.

Environment:
    GLMCC_HOME — path to GLMCC installation directory (default: auto-detect)
    GLMCC_CMD  — override the GLMCC command (default: python Est_Data.py)

Usage:
    python run_glmcc.py --animal day1 --condition lightON
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results" / "glmcc"


def _find_glmcc() -> Path | None:
    """Locate GLMCC installation. Checks in order:
    1. GLMCC_HOME environment variable
    2. ~/Documents/GLMCC
    3. ../GLMCC relative to project root
    """
    candidates = [
        os.environ.get("GLMCC_HOME"),
        os.path.expanduser("~/Documents/GLMCC"),
        str(PROJECT_ROOT.parent / "GLMCC"),
        os.path.expanduser("~/GLMCC"),
    ]
    for c in candidates:
        if c:
            p = Path(c)
            if p.exists():
                return p
    return None


def _count_units(spike_path: Path) -> int:
    """Count max columns (neurons) in spike-train file efficiently."""
    max_cols = 0
    with open(spike_path) as f:
        for line in f:
            if line.strip():
                n = len(line.strip().split())
                if n > max_cols:
                    max_cols = n
    return max_cols


def _write_placeholder(out_path: Path, n_units: int, reason: str,
                       metadata: dict | None = None) -> Path:
    """Write a zero-filled adjacency + metadata JSON explaining the failure."""
    placeholder = np.zeros((n_units, n_units))
    np.savetxt(str(out_path), placeholder, delimiter=",")

    # Write metadata alongside the CSV
    meta_path = out_path.with_suffix(".meta.json")
    meta = {
        "status": "placeholder",
        "reason": reason,
        "n_units": n_units,
        "command": metadata.get("command") if metadata else None,
        "stderr": metadata.get("stderr") if metadata else None,
    }
    meta_path.write_text(json.dumps(meta, indent=2))

    return out_path


def run_glmcc(animal: str, condition: str, glmcc_path: str | None = None) -> Path:
    # --- Resolve spike train file ---
    spike_path = (
        PROJECT_ROOT / "data" / "processed"
        / f"spiketrains_DATA{animal}_{condition}.txt"
    )
    if not spike_path.exists():
        # Try aliases
        alt_patterns = [
            f"spiketrains_DATA{animal}_evoked_LightON.txt",
            f"spiketrains_DATA{animal}*{condition}*.txt",
        ]
        for pat in alt_patterns:
            hits = sorted((PROJECT_ROOT / "data" / "processed").glob(pat))
            if hits:
                spike_path = hits[0]
                break
    if not spike_path.exists():
        raise FileNotFoundError(
            f"No spike-train file for animal={animal}, condition={condition}"
        )

    n_units = _count_units(spike_path)

    # --- Locate GLMCC ---
    glmcc_home = glmcc_path or _find_glmcc()
    if glmcc_home is None:
        reason = (
            "GLMCC not found. Set GLMCC_HOME environment variable "
            "or install GLMCC from https://github.com/naokimas/GLMCC"
        )
        print(f"WARNING: {reason}", file=sys.stderr)
        out = RESULTS_DIR / f"adj_{animal}_{condition}.csv"
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        return _write_placeholder(out, n_units, reason)

    glmcc_home = Path(glmcc_home)
    est_data = glmcc_home / "Est_Data.py"
    if not est_data.exists():
        reason = f"Est_Data.py not found in {glmcc_home}"
        print(f"WARNING: {reason}", file=sys.stderr)
        out = RESULTS_DIR / f"adj_{animal}_{condition}.csv"
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        return _write_placeholder(out, n_units, reason)

    # --- Run GLMCC ---
    cmd = [sys.executable, str(est_data), str(spike_path), "exp", "GLM"]
    cmd_str = " ".join(cmd)

    print(f"Running: {cmd_str}", flush=True)
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour max for large datasets
            cwd=str(RESULTS_DIR),
        )
    except subprocess.TimeoutExpired:
        reason = "GLMCC timed out (1 hour limit)"
        print(f"WARNING: {reason}", file=sys.stderr)
        out = RESULTS_DIR / f"adj_{animal}_{condition}.csv"
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        return _write_placeholder(out, n_units, reason, {"command": cmd_str})
    except Exception as e:
        reason = f"GLMCC subprocess failed: {e}"
        print(f"WARNING: {reason}", file=sys.stderr)
        out = RESULTS_DIR / f"adj_{animal}_{condition}.csv"
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        return _write_placeholder(out, n_units, reason, {"command": cmd_str})

    # --- Collect output ---
    out_path = RESULTS_DIR / f"adj_{animal}_{condition}.csv"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # GLMCC writes result_<animal>.csv to its working directory
    default_output = RESULTS_DIR / f"result_{animal}.csv"
    if default_output.exists():
        import shutil
        shutil.move(str(default_output), str(out_path))
        print(f"GLMCC result moved to {out_path}")
    elif result.returncode != 0:
        reason = f"GLMCC exited with code {result.returncode}"
        stderr_tail = result.stderr.strip()[-500:] if result.stderr else "(none)"
        print(f"WARNING: {reason}", file=sys.stderr)
        print(f"  stderr: {stderr_tail}", file=sys.stderr)
        _write_placeholder(out_path, n_units, reason, {
            "command": cmd_str,
            "stderr": stderr_tail,
        })
    else:
        # Check for result files in other locations
        alt_outputs = sorted(Path.cwd().glob(f"result*{animal}*.csv"))
        if alt_outputs:
            import shutil
            shutil.move(str(alt_outputs[0]), str(out_path))
            print(f"GLMCC result moved from {alt_outputs[0]} to {out_path}")
        else:
            reason = "GLMCC completed but no output CSV found"
            print(f"WARNING: {reason}", file=sys.stderr)
            _write_placeholder(out_path, n_units, reason, {
                "command": cmd_str,
                "stdout": result.stdout.strip()[-500:] if result.stdout else "(none)",
            })

    return out_path


def main():
    parser = argparse.ArgumentParser(description="Run GLMCC connectivity inference")
    parser.add_argument("--animal", required=True)
    parser.add_argument("--condition", required=True,
                        choices=["ongoing", "lightON", "ongoing_bis"])
    parser.add_argument("--glmcc", default=None,
                        help="Path to GLMCC installation directory (overrides GLMCC_HOME)")
    args = parser.parse_args()

    out = run_glmcc(args.animal, args.condition, args.glmcc)
    print(f"Output: {out}")


if __name__ == "__main__":
    main()
