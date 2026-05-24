#!/bin/bash
# GLMCC container entrypoint.
# Maps arguments to Est_Data.py call.
#
# Usage:
#   docker run --rm -v $(pwd)/data:/data -v $(pwd)/results:/results \
#       glmcc --animal 171019 --condition lightON
#
# Or via environment:
#   GLMCC_HOME=/opt/GLMCC python src/py/run_glmcc.py --animal 171019 --condition lightON

set -euo pipefail

ANIMAL=""
CONDITION=""
SPIKE_PATH=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --animal) ANIMAL="$2"; shift 2 ;;
        --condition) CONDITION="$2"; shift 2 ;;
        --spike-path) SPIKE_PATH="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

if [[ -z "$ANIMAL" || -z "$CONDITION" ]]; then
    echo "Usage: run_glmcc --animal <id> --condition <cond> [--spike-path <file>]"
    exit 1
fi

# Default spike train path
SPIKE_PATH="${SPIKE_PATH:-/data/processed/spiketrains_DATA${ANIMAL}_${CONDITION}.txt}"

if [[ ! -f "$SPIKE_PATH" ]]; then
    echo "ERROR: Spike train file not found: $SPIKE_PATH"
    exit 1
fi

echo "Running GLMCC on $SPIKE_PATH (animal=$ANIMAL, condition=$CONDITION)"

cd /results/glmcc
python /opt/GLMCC/Est_Data.py "$SPIKE_PATH" exp GLM

# Move result to expected location
if [[ -f "result_${ANIMAL}.csv" ]]; then
    mv "result_${ANIMAL}.csv" "adj_${ANIMAL}_${CONDITION}.csv"
    echo "Output: /results/glmcc/adj_${ANIMAL}_${CONDITION}.csv"
else
    echo "WARNING: GLMCC completed but no result_${ANIMAL}.csv found"
    ls -la /results/glmcc/
fi
