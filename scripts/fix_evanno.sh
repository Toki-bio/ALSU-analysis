#!/bin/bash
export PATH=/usr/local/bin:/usr/bin:/bin:/staging/conda/envs/bioinfo/bin
cd /staging/ALSU-analysis/admixture_analysis
BED=UZB_for_admixture
THREADS=16
LOG=evanno_runs/evanno_fix.log
echo "$(date): Fix evanno started" >> $LOG
for K in 7 8; do
  for REP in $(seq 1 10); do
    DIR=evanno_runs/K${K}_rep${REP}
    mkdir -p ${DIR}
    for ext in bed bim fam; do
      [ ! -e ${DIR}/${BED}.${ext} ] && ln -sf /staging/ALSU-analysis/admixture_analysis/${BED}.${ext} ${DIR}/${BED}.${ext}
    done
    if [ -f ${DIR}/${BED}.${K}.Q ]; then
      echo "SKIP K=${K} rep=${REP}" >> $LOG
      continue
    fi
    SEED=$((K*1000+REP))
    echo "$(date): RUN K=${K} rep=${REP} seed=${SEED}" >> $LOG
    cd ${DIR}
    admixture --cv -j${THREADS} -s ${SEED} ${BED}.bed ${K} > admixture_K${K}_rep${REP}.log 2>&1
    RC=$?
    echo "$(date): K=${K} rep=${REP} exit=${RC}" >> /staging/ALSU-analysis/admixture_analysis/$LOG
    cd /staging/ALSU-analysis/admixture_analysis
  done
done
echo "$(date): ALL FIX RUNS DONE" >> $LOG
