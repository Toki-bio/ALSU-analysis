#!/bin/bash
# Get the detailed 1000G sample→population mapping  
# Check in the 1000G reference directory and staging area
echo "=== LISTING_KG_DIR ==="
ls /staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/1000g_grch38/
echo "=== FIND_PANEL ==="
find /staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/ -name "*sample*" -o -name "*pop*" -o -name "*panel*" 2>/dev/null
echo "=== CHECK_PANEL ==="
for f in /staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/1000g_grch38/*.panel \
         /staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/1000g_grch38/pop_mapping* \
         /staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/1000g_grch38/sample* \
         /staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/*.panel \
         /staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/*.txt; do
  if [ -f "$f" ]; then
    echo "--- $f ---"
    head -3 "$f"
    echo "..."
    wc -l "$f"
  fi
done

# Also build mapping: for each ref IID in global_v2_admix.fam, get superpop from pca pop_mapping
# Then output CSV: IID,superpop,position  
echo "=== BUILD_MAPPING ==="
head -2548 ~/v2/global_admixture/global_v2_admix.fam | awk '{print $2}' | \
  awk 'NR==FNR{pop[$1]=$2;next} {print NR-1,$1,($1 in pop ? pop[$1] : "UNKNOWN")}' ~/v2/pca/pop_mapping.txt - | \
  tail -5

echo "=== SUPERPOP_BREAKDOWN ==="
# For each superpop, show the position ranges
head -2548 ~/v2/global_admixture/global_v2_admix.fam | awk '{print $2}' | \
  awk 'NR==FNR{pop[$1]=$2;next} {p=($1 in pop ? pop[$1] : "UNK"); printf "%d\t%s\n", NR-1, p}' ~/v2/pca/pop_mapping.txt - | \
  awk '{pop=$2; if(!(pop in first)) first[pop]=NR; last[pop]=NR; count[pop]++} END {for(p in count) print p, count[p], first[p]"-"last[p]}' | \
  sort -k2 -rn

echo "=== END ==="
