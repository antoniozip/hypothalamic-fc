#!/bin/bash
# DO NOT USE -- WRONG TIME BASE.
#
# This calls the binary with four arguments on DATA*/cell*.txt, which are in SECONDS.
# GLMCC reads WIN=50, DELTA=1 and tau=4 in the file's own unit, so it builds a
# +/-50 SECOND correlogram at 1 s resolution instead of +/-50 ms at 1 ms -- not a
# synaptic measurement. Anything it produced is invalid.
#
# Use `python3 scripts/regenerate_glmcc.py` instead: it converts to 0-based
# milliseconds and passes the true recording duration. See scripts/glmcc_units.py.
echo "run_glmcc_c_missing.sh is disabled: wrong time base. Use scripts/regenerate_glmcc.py" >&2
exit 1
#
# Original script retained below for reference only.
# Run C GLMCC for all missing ongoing conditions sequentially.
# day3 already done. Estimated total: ~13 hours for remaining 5.
set -euo pipefail
cd /media/antonio/data/tim_brown

GLMCC=vendor/glmcc-c/glmcc
mkdir -p results/glmcc

declare -A NEURONS=(
    [day6]=56
    [day4]=61
    [day5]=71
    [night7]=79
    [day1]=126
)

ORDER="day6 day4 day5 night7 day1"

for animal in $ORDER; do
    n=${NEURONS[$animal]}
    out="results/glmcc/adj_${animal}_ongoing.csv"
    
    # Skip if already has non-zero data
    if [ -f "$out" ]; then
        nz=$(python3 -c "import numpy as np; print(np.count_nonzero(np.loadtxt('$out',delimiter=',')))" 2>/dev/null)
        if [ "$nz" -gt 0 ] 2>/dev/null; then
            echo "SKIP $animal/ongoing ($nz non-zero edges)"
            continue
        fi
    fi
    
    echo "=== $animal/ongoing ($n neurons) ==="
    echo "Start: $(date)"
    
    $GLMCC "DATA${animal}_ongoing" $n exp GLM
    mv W_py_5400.csv "$out"
    
    nz=$(python3 -c "import numpy as np; print(np.count_nonzero(np.loadtxt('$out',delimiter=',')))" 2>/dev/null)
    echo "Done: $nz non-zero edges"
    echo "End: $(date)"
    echo ""
done

echo "=== ALL DONE: $(date) ==="
