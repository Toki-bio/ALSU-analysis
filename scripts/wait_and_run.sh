#!/bin/bash
# Waits for step6 to complete, then launches steps 7-11
# Run: nohup bash /tmp/wait_and_run.sh > /tmp/wait_and_run.log 2>&1 &
set -euo pipefail

POSTDIR=/staging/ALSU-analysis/spring2026/post_imputation
STEP6_LOG=/tmp/step6_out.txt
STEP7_11=/tmp/run_steps_7_to_11.sh

echo "[$(date)] Waiting for step6 to complete..."
echo "[$(date)] Looking for UZB_imputed_HQ_clean.fam in ${POSTDIR}/"

while true; do
    if [ -f "${POSTDIR}/UZB_imputed_HQ_clean.fam" ]; then
        # Double-check the step6 log says DONE
        if grep -q "=== DONE ===" "$STEP6_LOG" 2>/dev/null; then
            echo "[$(date)] Step 6 complete! Launching steps 7-11..."
            break
        fi
    fi
    sleep 300  # Check every 5 minutes
done

# Wait 30s for any file I/O to flush
sleep 30

echo "[$(date)] Step 6 output:"
tail -10 "$STEP6_LOG"
echo ""

# Launch steps 7-11
bash "$STEP7_11" > /tmp/step7_11_out.txt 2>&1
echo "[$(date)] Steps 7-11 pipeline finished. Check /tmp/step7_11_out.txt"
