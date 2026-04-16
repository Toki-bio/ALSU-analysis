#!/bin/bash
set -euo pipefail
ADMDIR=/staging/ALSU-analysis/spring2026/admixture
POSTDIR=/staging/ALSU-analysis/spring2026/post_imputation
SPRING=/staging/ALSU-analysis/spring2026

echo "===GLOBAL_CV==="
for k in 2 3 4 5 6 7 8; do
  echo -n "K=$k: "
  grep -i 'CV\|cross-validation' "$ADMDIR/admixture_K${k}.log" 2>/dev/null || echo "N/A"
done

echo "===UZB_CV==="
for k in 2 3 4 5 6 7 8; do
  echo -n "K=$k: "
  grep -i 'CV\|cross-validation' "$ADMDIR/uzb_admixture_K${k}.log" 2>/dev/null || echo "N/A"
done

echo "===GLOBAL_LOG_TAILS==="
for k in 2 3 4 5 6 7 8; do
  echo "--- Global K=$k ---"
  tail -10 "$ADMDIR/admixture_K${k}.log" 2>/dev/null || echo "N/A"
done

echo "===UZB_LOG_TAILS==="
for k in 2 3; do
  echo "--- UZB K=$k ---"
  tail -10 "$ADMDIR/uzb_admixture_K${k}.log" 2>/dev/null || echo "N/A"
done

echo "===STEP7_LOCAL_PCA==="
ls "$POSTDIR"/*local* "$POSTDIR"/*pca* "$POSTDIR"/*eigenval* "$POSTDIR"/*eigenvec* 2>/dev/null || echo "no local PCA files in postdir"
ls "$SPRING"/pca/ 2>/dev/null && echo "pca dir found" || echo "no pca dir"
ls "$SPRING"/local_pca/ 2>/dev/null && echo "local_pca dir found" || echo "no local_pca dir"
find "$SPRING" -maxdepth 3 -name "*.eigenval" -o -name "*.eigenvec" -o -name "*pca*log" 2>/dev/null | head -20

echo "===STEP8_GLOBAL_PCA==="
find "$SPRING" -maxdepth 3 -name "*global*pca*" -o -name "*merged*" -o -name "*1000g*" 2>/dev/null | head -20

echo "===STEP9_FST==="
ls "$SPRING"/fst/ 2>/dev/null || echo "no fst dir"
find "$SPRING" -maxdepth 3 -name "*fst*" -o -name "*Fst*" 2>/dev/null | head -20

echo "===STEP10_PBS==="
ls "$SPRING"/pbs/ 2>/dev/null || echo "no pbs dir"
find "$SPRING" -maxdepth 3 -name "*pbs*" 2>/dev/null | head -20

echo "===STEP11_ADMIX_DETAILS==="
cat "$ADMDIR/global_for_admixture.log" 2>/dev/null
echo "---"
cat "$ADMDIR/uzb_for_admixture.log" 2>/dev/null
echo "---"
cat "$ADMDIR/global_ld.log" 2>/dev/null
echo "---"
cat "$ADMDIR/uzb_ld.log" 2>/dev/null
echo "---"
cat "$ADMDIR/global_qc.log" 2>/dev/null
echo "---"
cat "$ADMDIR/uzb_qc.log" 2>/dev/null

echo "===STEP6_HET==="
cat "$POSTDIR/het_outliers.txt" 2>/dev/null || echo "no het outliers file"
echo "---"
head -5 "$POSTDIR/UZB_imputed_HQ.het" 2>/dev/null || echo "no het file"
echo "---"
wc -l "$POSTDIR/UZB_imputed_HQ.het" 2>/dev/null || echo "no het file"

echo "===DISK==="
df -h /staging | tail -1
echo "===DONE==="
