#!/usr/bin/env bash
set -euo pipefail

BASE=/staging/ALSU-analysis/admixture_analysis
RUNS=$BASE/evanno_runs
ADMIXTURE=/staging/conda/envs/bioinfo/bin/admixture
LOG=$RUNS/evanno_k78_final.log

mkdir -p "$RUNS"

echo "[$(date)] Evanno K7+K8 final run started" | tee "$LOG"

for K in 7 8; do
  for REP in 1 2 3 4 5 6 7 8 9 10; do
    SEED=$((K * 1000 + REP))
    DIR="$RUNS/K${K}_rep${REP}"
    
    mkdir -p "$DIR"
    cd "$DIR"
    ln -sf "$BASE/UZB_for_admixture.bed" .
    ln -sf "$BASE/UZB_for_admixture.bim" .
    ln -sf "$BASE/UZB_for_admixture.fam" .
    
    Q="UZB_for_admixture.${K}.Q"
    P="UZB_for_admixture.${K}.P"
    
    if [[ -s "$Q" && -s "$P" ]]; then
      echo "[$(date)] K=$K rep=$REP SKIP (exists)" | tee -a "$LOG"
      continue
    fi
    
    echo "[$(date)] K=$K rep=$REP seed=$SEED START" | tee -a "$LOG"
    "$ADMIXTURE" --cv UZB_for_admixture.bed "$K" -j16 -s "$SEED" > admixture_K${K}_rep${REP}.log 2>&1
    echo "[$(date)] K=$K rep=$REP DONE" | tee -a "$LOG"
  done
done

{
  echo "=== FINAL PER-K COMPLETION ==="
  for K in 2 3 4 5 6 7 8; do
    dirs=$(find "$RUNS" -maxdepth 1 -type d -name "K${K}_rep*" | wc -l)
    logs=$(find "$RUNS" -maxdepth 2 -type f -path "*/K${K}_rep*/*.log" | wc -l)
    q=$(find "$RUNS" -maxdepth 2 -type f -path "*/K${K}_rep*/UZB_for_admixture.${K}.Q" | wc -l)
    p=$(find "$RUNS" -maxdepth 2 -type f -path "*/K${K}_rep*/UZB_for_admixture.${K}.P" | wc -l)
    echo "K${K}: dirs=${dirs} logs=${logs} Q=${q} P=${p}"
  done
  echo "=== done $(date) ==="
} | tee -a "$LOG"

echo "[$(date)] Evanno K7+K8 final run COMPLETE" | tee -a "$LOG"
