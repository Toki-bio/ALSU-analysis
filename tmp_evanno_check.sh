#!/usr/bin/env bash
set -euo pipefail

cd /staging/ALSU-analysis/admixture_analysis/evanno_runs

echo "=== CONNECT ==="
hostname
date

echo "=== PROCESS ==="
pgrep -af 'evanno_k78|admixture --cv UZB_for_admixture.bed [78]' || true

echo "=== COUNTS ==="
for K in 7 8; do
  q=$(find . -maxdepth 2 -type f -path "./K${K}_rep*/UZB_for_admixture.${K}.Q" | wc -l)
  p=$(find . -maxdepth 2 -type f -path "./K${K}_rep*/UZB_for_admixture.${K}.P" | wc -l)
  logs=$(find . -maxdepth 2 -type f -path "./K${K}_rep*/*.log" | wc -l)
  echo "K${K}: Q=${q} P=${p} logs=${logs}"
done

echo "=== LOG TAIL ==="
tail -n 30 evanno_k78_final.log 2>/dev/null || tail -n 30 evanno_k78.nohup 2>/dev/null || true
