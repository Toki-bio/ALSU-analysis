#!/bin/bash
cd /staging/ALSU-analysis/admixture_analysis/evanno_runs
echo "=== Current status ==="
for K in 2 3 4 5 6 7 8; do
  q=$(find . -maxdepth 2 -type f -path "./K${K}_rep*/UZB_for_admixture.${K}.Q" | wc -l)
  echo "K${K}: Q=${q}/10"
done
tail -n 30 evanno_k78_final.log 2>/dev/null || true
