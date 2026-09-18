#!/usr/bin/env bash
# Regenerate graph metrics and first/second-level statistics from the
# re-validated adjacency matrices.
#
# Run after scripts/revalidate_all.py. Any animal x condition whose adjacency
# could not be produced is skipped here too, so a stale matrix in the old
# encoding never reaches the models.
set -uo pipefail
cd "$(dirname "$0")/.."

ESTIMATOR="${ESTIMATOR:-glmcc}"
CONDITIONS=(lightON ongoing)
ANIMALS=(day1 day2 day3 day4 day5 day6 night1 night2 night3 night4 night5 night6 night7)

ok=0; skipped=0; failed=0
for condition in "${CONDITIONS[@]}"; do
  for animal in "${ANIMALS[@]}"; do
    adj="results/${ESTIMATOR}/validated_adj_${animal}_${condition}.csv"
    if [[ ! -f "$adj" ]]; then
      echo "SKIP  ${animal}_${condition}: no validated adjacency"
      skipped=$((skipped+1)); continue
    fi
    # Existence is not enough: a pair the re-validation could not produce may
    # still have a file left over from the old split pipeline, in the old
    # encoding. The sidecar written beside each matrix states what it holds.
    meta="results/${ESTIMATOR}/validated_adj_${animal}_${condition}.meta.json"
    if [[ ! -f "$meta" ]] || ! grep -q "weight_masked" "$meta"; then
      echo "SKIP  ${animal}_${condition}: no current-encoding provenance sidecar"
      skipped=$((skipped+1)); continue
    fi
    if Rscript src/r/compute_graph_metrics.R --animal "$animal" --condition "$condition" \
         --estimator "$ESTIMATOR" >/dev/null 2>"logs/metrics_${animal}_${condition}.err"; then
      echo "OK    ${animal}_${condition}"
      ok=$((ok+1))
    else
      echo "FAIL  ${animal}_${condition}: $(tail -2 "logs/metrics_${animal}_${condition}.err" | tr '\n' ' ')"
      failed=$((failed+1))
    fi
    Rscript src/r/aggregate_density.R --animal "$animal" --condition "$condition" \
      --estimator "$ESTIMATOR" >/dev/null 2>&1 || true
  done
done

echo
echo "metrics: ${ok} ok, ${skipped} skipped, ${failed} failed"
echo
echo "=== first-level LME ==="
Rscript src/r/stats_first_level.R --estimator "$ESTIMATOR" 2>&1 | tail -25
echo
echo "=== second-level LME ==="
Rscript src/r/stats_second_level.R --estimator "$ESTIMATOR" 2>&1 | tail -15
