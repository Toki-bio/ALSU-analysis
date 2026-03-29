#!/usr/bin/env bash
set -euo pipefail
cd /staging/ALSU-analysis/admixture_analysis/evanno_runs

echo "=== NOW ==="
date

echo "=== ACTIVE ADMIXTURE PROCESSES ==="
pgrep -af 'admixture --cv UZB_for_admixture.bed [2-8]' || true

echo "=== PID DETAILS ==="
for pid in $(pgrep -f 'admixture --cv UZB_for_admixture.bed [2-8]' || true); do
  ps -p "$pid" -o pid,etime,%cpu,%mem,stat,cmd
  echo "---"
done

echo "=== K7 REP4 FILE STATUS ==="
ls -l --time-style=long-iso K7_rep4/admixture_K7_rep4.log 2>/dev/null || true
ls -l --time-style=long-iso K7_rep4/UZB_for_admixture.7.Q 2>/dev/null || true
ls -l --time-style=long-iso K7_rep4/UZB_for_admixture.7.P 2>/dev/null || true

echo "=== K8 REP1 FILE STATUS ==="
ls -l --time-style=long-iso K8_rep1/admixture_K8_rep1.log 2>/dev/null || true
ls -l --time-style=long-iso K8_rep1/UZB_for_admixture.8.Q 2>/dev/null || true
ls -l --time-style=long-iso K8_rep1/UZB_for_admixture.8.P 2>/dev/null || true

echo "=== K7 REP4 LOG (last 30) ==="
tail -n 30 K7_rep4/admixture_K7_rep4.log 2>/dev/null || true

echo "=== QP MASTER LOG TAIL ==="
tail -n 40 evanno_qp_nohup.out 2>/dev/null || true

echo "=== COUNTS K2..K8 ==="
for K in 2 3 4 5 6 7 8; do
  q=$(find . -maxdepth 2 -type f -path "./K${K}_rep*/UZB_for_admixture.${K}.Q" | wc -l)
  p=$(find . -maxdepth 2 -type f -path "./K${K}_rep*/UZB_for_admixture.${K}.P" | wc -l)
  logs=$(find . -maxdepth 2 -type f -path "./K${K}_rep*/admixture_K${K}_rep*.log" | wc -l)
  echo "K${K}: logs=${logs} Q=${q} P=${p}"
done
