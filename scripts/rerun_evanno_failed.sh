#!/bin/bash
set -e
cd /staging/ALSU-analysis/admixture_analysis

ADMIXTURE=/staging/conda/envs/bioinfo/bin/admixture
LOG=evanno_runs/evanno_rerun.log
echo "Started rerun: $(date)" > $LOG
echo "Using binary: $ADMIXTURE" >> $LOG

# K=7 reps 2-10 (rep1 already exists)
for REP in 2 3 4 5 6 7 8 9 10; do
  SEED=$((7 * 1000 + REP))
  OUTDIR="evanno_runs/K7_rep${REP}"
  mkdir -p "$OUTDIR"
  cd "$OUTDIR"
  for ext in bed bim fam; do
    ln -sf /staging/ALSU-analysis/admixture_analysis/UZB_for_admixture.$ext .
  done
  echo "$(date): K=7 rep=$REP seed=$SEED" >> /staging/ALSU-analysis/admixture_analysis/$LOG
  $ADMIXTURE --cv UZB_for_admixture.bed 7 -j16 -s $SEED >> /staging/ALSU-analysis/admixture_analysis/$LOG 2>&1
  echo "$(date): K=7 rep=$REP DONE" >> /staging/ALSU-analysis/admixture_analysis/$LOG
  cd /staging/ALSU-analysis/admixture_analysis
done

# K=8 all 10 reps
for REP in 1 2 3 4 5 6 7 8 9 10; do
  SEED=$((8 * 1000 + REP))
  OUTDIR="evanno_runs/K8_rep${REP}"
  mkdir -p "$OUTDIR"
  cd "$OUTDIR"
  for ext in bed bim fam; do
    ln -sf /staging/ALSU-analysis/admixture_analysis/UZB_for_admixture.$ext .
  done
  echo "$(date): K=8 rep=$REP seed=$SEED" >> /staging/ALSU-analysis/admixture_analysis/$LOG
  $ADMIXTURE --cv UZB_for_admixture.bed 8 -j16 -s $SEED >> /staging/ALSU-analysis/admixture_analysis/$LOG 2>&1
  echo "$(date): K=8 rep=$REP DONE" >> /staging/ALSU-analysis/admixture_analysis/$LOG
  cd /staging/ALSU-analysis/admixture_analysis
done

echo "ALL RERUNS COMPLETE: $(date)" >> $LOG
echo "Done! Check $LOG for results."
