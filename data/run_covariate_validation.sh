#!/bin/bash
# ─────────────────────────────────────────────────────────────────────
# run_covariate_validation.sh
# Waits for ADMIXTURE to finish, then runs covariate validation
# ─────────────────────────────────────────────────────────────────────
set -euo pipefail

BASE="/staging/ALSU-analysis/admixture_analysis"
SCRIPT="$BASE/validate_admixture_covariates.py"
LOG="$BASE/validation_runner.log"

echo "╔══════════════════════════════════════════════════════════════╗" | tee "$LOG"
echo "║  ADMIXTURE Covariate Validation Scheduler                   ║" | tee -a "$LOG"
echo "║  Started: $(date)                           ║" | tee -a "$LOG"
echo "╚══════════════════════════════════════════════════════════════╝" | tee -a "$LOG"

# ── 1. Wait for Uzbek-only ADMIXTURE (K=2) ──
UZB_Q2="$BASE/UZB_for_admixture.2.Q"
if [ -f "$UZB_Q2" ]; then
    echo "[$(date +%H:%M)] Uzbek K=2 Q file already exists." | tee -a "$LOG"
else
    echo "[$(date +%H:%M)] Waiting for Uzbek ADMIXTURE K=2..." | tee -a "$LOG"
    while [ ! -f "$UZB_Q2" ]; do
        sleep 60
    done
    echo "[$(date +%H:%M)] Uzbek K=2 Q file found." | tee -a "$LOG"
fi

# ── 2. Run Uzbek-only validation ──
echo "[$(date +%H:%M)] Running Uzbek-only covariate validation..." | tee -a "$LOG"
python3 "$SCRIPT" 2>&1 | tee -a "$LOG"
echo "[$(date +%H:%M)] Uzbek validation complete." | tee -a "$LOG"

# ── 3. Check for Global ADMIXTURE ──
GLOBAL_Q2="$BASE/global_admixture/admix_results/global_for_admixture.2.Q"
GLOBAL_FAM="$BASE/global_admixture/global_for_admixture.fam"

if [ -f "$GLOBAL_FAM" ]; then
    echo "[$(date +%H:%M)] Global ADMIXTURE FAM file exists." | tee -a "$LOG"

    # Wait for at least K=2 to finish
    if [ -f "$GLOBAL_Q2" ]; then
        echo "[$(date +%H:%M)] Global K=2 Q file already exists." | tee -a "$LOG"
    else
        echo "[$(date +%H:%M)] Waiting for global ADMIXTURE K=2..." | tee -a "$LOG"
        # Wait up to 48 hours (global takes much longer)
        WAIT_COUNT=0
        MAX_WAIT=2880  # 48 hours at 60s intervals
        while [ ! -f "$GLOBAL_Q2" ] && [ $WAIT_COUNT -lt $MAX_WAIT ]; do
            sleep 60
            WAIT_COUNT=$((WAIT_COUNT + 1))
            if [ $((WAIT_COUNT % 60)) -eq 0 ]; then
                echo "[$(date +%H:%M)] Still waiting for global K=2... (${WAIT_COUNT} min)" | tee -a "$LOG"
            fi
        done

        if [ ! -f "$GLOBAL_Q2" ]; then
            echo "[$(date +%H:%M)] WARNING: Timed out waiting for global ADMIXTURE." | tee -a "$LOG"
            echo "[$(date +%H:%M)] Run manually later: python3 $SCRIPT --global" | tee -a "$LOG"
            exit 0
        fi
    fi

    # Run global validation
    echo "[$(date +%H:%M)] Running global covariate validation..." | tee -a "$LOG"
    python3 "$SCRIPT" --global 2>&1 | tee -a "$LOG"
    echo "[$(date +%H:%M)] Global validation complete." | tee -a "$LOG"
else
    echo "[$(date +%H:%M)] Global ADMIXTURE not started yet. Skipping global validation." | tee -a "$LOG"
    echo "[$(date +%H:%M)] Run later: python3 $SCRIPT --global" | tee -a "$LOG"
fi

echo "" | tee -a "$LOG"
echo "╔══════════════════════════════════════════════════════════════╗" | tee -a "$LOG"
echo "║  All validation tasks complete!                             ║" | tee -a "$LOG"
echo "║  Finished: $(date)                          ║" | tee -a "$LOG"
echo "╚══════════════════════════════════════════════════════════════╝" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "Output files:" | tee -a "$LOG"
echo "  Uzbek:  $BASE/validation/" | tee -a "$LOG"
echo "  Global: $BASE/global_admixture/validation/" | tee -a "$LOG"
