#!/usr/bin/env bash
set -euo pipefail

BASE=/staging/ALSU-analysis/admixture_analysis
RUNS=$BASE/evanno_runs
ADMIXTURE=/staging/conda/envs/bioinfo/bin/admixture
MASTER_LOG=$RUNS/evanno_qp_complete.log
SUMMARY_TSV=$RUNS/evanno_qp_completion_summary.tsv
STATUS_TXT=$RUNS/evanno_qp_completion_status.txt

mkdir -p "$RUNS"

echo "[$(date)] Evanno Q/P completion started" | tee "$MASTER_LOG"
printf 'K\trep\tseed\taction\n' > "$SUMMARY_TSV"

for K in 2 3 4 5 6 7 8; do
  for REP in 1 2 3 4 5 6 7 8 9 10; do
    SEED=$((K * 1000 + REP))
    DIR="$RUNS/K${K}_rep${REP}"
    REP_LOG="admixture_K${K}_rep${REP}.log"

    mkdir -p "$DIR"
    cd "$DIR"
    ln -sf "$BASE/UZB_for_admixture.bed" .
    ln -sf "$BASE/UZB_for_admixture.bim" .
    ln -sf "$BASE/UZB_for_admixture.fam" .

    Q="UZB_for_admixture.${K}.Q"
    P="UZB_for_admixture.${K}.P"

    if [[ -s "$Q" && -s "$P" ]]; then
      printf '%s\t%s\t%s\tskip_qp_exists\n' "$K" "$REP" "$SEED" >> "$SUMMARY_TSV"
      continue
    fi

    echo "[$(date)] RUN K=$K rep=$REP seed=$SEED (missing Q/P)" | tee -a "$MASTER_LOG"
    "$ADMIXTURE" --cv UZB_for_admixture.bed "$K" -j16 -s "$SEED" > "$REP_LOG" 2>&1
    echo "[$(date)] DONE K=$K rep=$REP" | tee -a "$MASTER_LOG"
    printf '%s\t%s\t%s\tran\n' "$K" "$REP" "$SEED" >> "$SUMMARY_TSV"
  done
done

{
  echo "=== Per-K completion ==="
  for K in 2 3 4 5 6 7 8; do
    dirs=$(find "$RUNS" -maxdepth 1 -type d -name "K${K}_rep*" | wc -l)
    logs=$(find "$RUNS" -maxdepth 2 -type f -path "*/K${K}_rep*/*.log" | wc -l)
    q=$(find "$RUNS" -maxdepth 2 -type f -path "*/K${K}_rep*/UZB_for_admixture.${K}.Q" | wc -l)
    p=$(find "$RUNS" -maxdepth 2 -type f -path "*/K${K}_rep*/UZB_for_admixture.${K}.P" | wc -l)
    echo "K${K}: dirs=${dirs} logs=${logs} Q=${q} P=${p}"
  done
  echo "=== done $(date) ==="
} | tee "$STATUS_TXT"

echo "[$(date)] Evanno Q/P completion finished" | tee -a "$MASTER_LOG"
