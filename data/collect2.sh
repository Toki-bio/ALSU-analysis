#!/bin/bash
set -euo pipefail
SPRING=/staging/ALSU-analysis/spring2026

echo "===STEP7_PCA_EIGENVAL==="
cat "$SPRING/pca/UZB_final_pca.eigenval"

echo "===STEP7_PCA_LOG==="
cat "$SPRING/pca/UZB_final_pca.log"

echo "===STEP7_QC_LOG==="
cat "$SPRING/pca/UZB_qc.log"

echo "===STEP7_UNIQUE_LOG==="
cat "$SPRING/pca/UZB_unique.log"

echo "===STEP7_PRUNE_LOG==="
cat "$SPRING/pca/UZB_pruned.log"

echo "===STEP7_LD_LOG==="
cat "$SPRING/pca/UZB_ldpruned.log"

echo "===STEP7_EIGENVEC_HEAD==="
head -5 "$SPRING/pca/UZB_final_pca.eigenvec"
echo "..."
wc -l "$SPRING/pca/UZB_final_pca.eigenvec"

echo "===STEP8_GLOBAL_PCA==="
cat "$SPRING/global_pca/GLOBAL_PCA.eigenval" 2>/dev/null || echo "no eigenval"
echo "---"
head -5 "$SPRING/global_pca/GLOBAL_PCA.eigenvec" 2>/dev/null || echo "no eigenvec"
echo "..."
wc -l "$SPRING/global_pca/GLOBAL_PCA.eigenvec" 2>/dev/null || echo "none"
echo "---MERGED_LOG---"
cat "$SPRING/global_pca/UZB_1kG_merged.log" 2>/dev/null | head -30
echo "---PCA_FILES---"
ls -la "$SPRING/global_pca/"

echo "===STEP9_FST_RESULTS==="
for pair in UZB_EUR UZB_EAS UZB_SAS UZB_AFR EUR_EAS; do
  echo "--- $pair ---"
  grep -i 'weighted\|mean' "$SPRING/fst/fst_${pair}.log" 2>/dev/null || echo "N/A"
done

echo "===STEP9_POP_COUNTS==="
for pop in UZB EUR EAS SAS AFR; do
  echo -n "$pop: "
  wc -l < "$SPRING/fst/keep_${pop}.txt" 2>/dev/null || echo "N/A"
done

echo "===STEP10_PBS_STATS==="
cat "$SPRING/pbs/pbs_stats.json" 2>/dev/null || echo "no stats"
echo ""
echo "===STEP10_PBS_CANDIDATES==="
cat "$SPRING/pbs/pbs_candidates.json" 2>/dev/null | head -50 || echo "no candidates"
echo ""
echo "===STEP10_PBS_HEAD==="
head -20 "$SPRING/pbs/pbs_all.tsv" 2>/dev/null || echo "no tsv"
wc -l "$SPRING/pbs/pbs_all.tsv" 2>/dev/null

echo "===STEP8_PCA_LOG==="
find "$SPRING/global_pca" -name "*.log" -exec echo "--- {} ---" \; -exec cat {} \; 2>/dev/null

echo "===DONE==="
