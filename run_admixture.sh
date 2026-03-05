#!/bin/bash
cd /staging/ALSU-analysis/admixture_analysis
echo "Started: $(date)" > admixture_all.log
for K in 2 3 4 5 6 7 8; do
  echo "=== Running K=$K ===" >> admixture_all.log
  admixture --cv UZB_for_admixture.bed $K -j32 >> admixture_all.log 2>&1
  echo "=== K=$K DONE ===" >> admixture_all.log
done
echo "ALL_COMPLETE: $(date)" >> admixture_all.log
