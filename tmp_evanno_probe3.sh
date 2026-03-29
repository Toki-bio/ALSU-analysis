#!/usr/bin/env bash
set -euo pipefail
BASE=/staging/ALSU-analysis/admixture_analysis/evanno_runs

echo "=== BASE ==="
echo "$BASE"
[ -d "$BASE" ] || { echo "MISSING_BASE"; exit 0; }

for k in 2 3 4 5 6 7 8; do
  total_dirs=$(find "$BASE" -maxdepth 1 -type d -name "K${k}_rep*" | wc -l)
  logs_present=$(find "$BASE" -maxdepth 2 -type f -path "$BASE/K${k}_rep*/admixture_K${k}_rep*.log" | wc -l)
  q_present=$(find "$BASE" -maxdepth 2 -type f -path "$BASE/K${k}_rep*/UZB_for_admixture.${k}.Q" | wc -l)
  p_present=$(find "$BASE" -maxdepth 2 -type f -path "$BASE/K${k}_rep*/UZB_for_admixture.${k}.P" | wc -l)
  echo "K=$k dirs=$total_dirs logs=$logs_present Q=$q_present P=$p_present"
done

echo "=== MISSING LOG DIRS ==="
for k in 2 3 4 5 6 7 8; do
  for d in "$BASE"/K${k}_rep*; do
    [ -d "$d" ] || continue
    rep=${d##*rep}
    log="$d/admixture_K${k}_rep${rep}.log"
    [ -f "$log" ] || echo "missing_log K${k}_rep${rep}"
  done
done

echo "=== DONE ==="
