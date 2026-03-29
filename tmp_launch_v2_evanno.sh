#!/usr/bin/env bash
set -euo pipefail

echo "=== KILL STALE V1 EVANNO ==="
# Kill the runner and any child admixture processes
pkill -f 'tmp_evanno_complete_qp' 2>/dev/null && echo "Killed runner" || echo "Runner not found"
pkill -f 'admixture --cv UZB_for_admixture.bed' 2>/dev/null && echo "Killed admixture" || echo "No admixture process"
sleep 2
pgrep -af admixture || echo "All admixture processes stopped"

echo ""
echo "=== SETUP V2 EVANNO: GLOBAL ==="
GLOBAL_BASE=/home/copilot/v2/global_admixture
GLOBAL_EVANNO=/home/copilot/v2/global_admixture/evanno
ADMIXTURE=/staging/conda/envs/bioinfo/bin/admixture

mkdir -p "$GLOBAL_EVANNO"

echo "Input: global_v2_admix — $(wc -l < $GLOBAL_BASE/global_v2_admix.fam) samples, $(wc -l < $GLOBAL_BASE/global_v2_admix.bim) SNPs"

echo ""
echo "=== SETUP V2 EVANNO: UZB-ONLY ==="
UZB_BASE=/home/copilot/v2/admixture
UZB_EVANNO=/home/copilot/v2/admixture/evanno

mkdir -p "$UZB_EVANNO"

echo "Input: UZB_v2_admix — $(wc -l < $UZB_BASE/UZB_v2_admix.fam) samples, $(wc -l < $UZB_BASE/UZB_v2_admix.bim) SNPs"

echo ""
echo "=== LAUNCH GLOBAL EVANNO (K2-K8 x 10 reps) ==="

cat > /home/copilot/v2/run_evanno_global.sh << 'EOFG'
#!/usr/bin/env bash
set -euo pipefail
ADMIXTURE=/staging/conda/envs/bioinfo/bin/admixture
BASE=/home/copilot/v2/global_admixture
RUNS=/home/copilot/v2/global_admixture/evanno
LOG=$RUNS/evanno_global.log

echo "[$(date)] V2 Global Evanno started" | tee "$LOG"
echo "Input: $(wc -l < $BASE/global_v2_admix.fam) samples, $(wc -l < $BASE/global_v2_admix.bim) SNPs" | tee -a "$LOG"

for K in 2 3 4 5 6 7 8; do
  for REP in 1 2 3 4 5 6 7 8 9 10; do
    SEED=$((K * 1000 + REP))
    DIR="$RUNS/K${K}_rep${REP}"
    mkdir -p "$DIR"
    cd "$DIR"
    ln -sf "$BASE/global_v2_admix.bed" .
    ln -sf "$BASE/global_v2_admix.bim" .
    ln -sf "$BASE/global_v2_admix.fam" .

    Q="global_v2_admix.${K}.Q"
    P="global_v2_admix.${K}.P"

    if [[ -s "$Q" && -s "$P" ]]; then
      echo "[$(date)] K=$K rep=$REP SKIP (exists)" | tee -a "$LOG"
      continue
    fi

    echo "[$(date)] K=$K rep=$REP seed=$SEED START" | tee -a "$LOG"
    "$ADMIXTURE" --cv global_v2_admix.bed "$K" -j16 -s "$SEED" > "admixture_K${K}_rep${REP}.log" 2>&1
    echo "[$(date)] K=$K rep=$REP DONE" | tee -a "$LOG"
  done
done

echo ""
echo "=== FINAL COUNTS ===" | tee -a "$LOG"
for K in 2 3 4 5 6 7 8; do
  q=$(find "$RUNS" -maxdepth 2 -type f -path "*/K${K}_rep*/global_v2_admix.${K}.Q" | wc -l)
  p=$(find "$RUNS" -maxdepth 2 -type f -path "*/K${K}_rep*/global_v2_admix.${K}.P" | wc -l)
  echo "K${K}: Q=${q} P=${p}" | tee -a "$LOG"
done
echo "[$(date)] V2 Global Evanno COMPLETE" | tee -a "$LOG"
EOFG
chmod +x /home/copilot/v2/run_evanno_global.sh

echo ""
echo "=== LAUNCH UZB-ONLY EVANNO (K2-K8 x 10 reps) ==="

cat > /home/copilot/v2/run_evanno_uzb.sh << 'EOFU'
#!/usr/bin/env bash
set -euo pipefail
ADMIXTURE=/staging/conda/envs/bioinfo/bin/admixture
BASE=/home/copilot/v2/admixture
RUNS=/home/copilot/v2/admixture/evanno
LOG=$RUNS/evanno_uzb.log

echo "[$(date)] V2 UZB-only Evanno started" | tee "$LOG"
echo "Input: $(wc -l < $BASE/UZB_v2_admix.fam) samples, $(wc -l < $BASE/UZB_v2_admix.bim) SNPs" | tee -a "$LOG"

for K in 2 3 4 5 6 7 8; do
  for REP in 1 2 3 4 5 6 7 8 9 10; do
    SEED=$((K * 1000 + REP))
    DIR="$RUNS/K${K}_rep${REP}"
    mkdir -p "$DIR"
    cd "$DIR"
    ln -sf "$BASE/UZB_v2_admix.bed" .
    ln -sf "$BASE/UZB_v2_admix.bim" .
    ln -sf "$BASE/UZB_v2_admix.fam" .

    Q="UZB_v2_admix.${K}.Q"
    P="UZB_v2_admix.${K}.P"

    if [[ -s "$Q" && -s "$P" ]]; then
      echo "[$(date)] K=$K rep=$REP SKIP (exists)" | tee -a "$LOG"
      continue
    fi

    echo "[$(date)] K=$K rep=$REP seed=$SEED START" | tee -a "$LOG"
    "$ADMIXTURE" --cv UZB_v2_admix.bed "$K" -j16 -s "$SEED" > "admixture_K${K}_rep${REP}.log" 2>&1
    echo "[$(date)] K=$K rep=$REP DONE" | tee -a "$LOG"
  done
done

echo ""
echo "=== FINAL COUNTS ===" | tee -a "$LOG"
for K in 2 3 4 5 6 7 8; do
  q=$(find "$RUNS" -maxdepth 2 -type f -path "*/K${K}_rep*/UZB_v2_admix.${K}.Q" | wc -l)
  p=$(find "$RUNS" -maxdepth 2 -type f -path "*/K${K}_rep*/UZB_v2_admix.${K}.P" | wc -l)
  echo "K${K}: Q=${q} P=${p}" | tee -a "$LOG"
done
echo "[$(date)] V2 UZB-only Evanno COMPLETE" | tee -a "$LOG"
EOFU
chmod +x /home/copilot/v2/run_evanno_uzb.sh

# Launch global first (higher priority), then UZB sequentially after
nohup bash -c 'bash /home/copilot/v2/run_evanno_global.sh; bash /home/copilot/v2/run_evanno_uzb.sh' > /home/copilot/v2/evanno_v2_nohup.out 2>&1 < /dev/null &
echo "LAUNCHED PID: $!"

sleep 3
echo ""
echo "=== VERIFY LAUNCH ==="
pgrep -af 'run_evanno|admixture' || echo "Nothing running yet"
tail -n 5 /home/copilot/v2/evanno_v2_nohup.out 2>/dev/null || true
