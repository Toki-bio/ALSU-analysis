#!/bin/bash
# Get detailed pop mapping for the 2548 ref samples in global ADMIXTURE
# The FAM has all 3595 samples; crop to first 2548 (the ref ones)
# Cross-reference with PCA pop_mapping

echo "=== REF_FAM_IDS ==="
head -2548 ~/v2/global_admixture/global_v2_admix.fam | awk '{print $2}' > /tmp/ref_ids.txt
wc -l /tmp/ref_ids.txt

echo "=== MATCH_POP ==="
# Join ref IDs with PCA pop_mapping (sample_id TAB pop)
# pop_mapping format: IID TAB SUPERPOP
awk 'NR==FNR{ids[$1]=1;next} $1 in ids{print $1,$2}' /tmp/ref_ids.txt ~/v2/pca/pop_mapping.txt | \
  awk '{print $2}' | sort | uniq -c | sort -rn

echo "=== MATCH_COUNT ==="
awk 'NR==FNR{ids[$1]=1;next} $1 in ids{print}' /tmp/ref_ids.txt ~/v2/pca/pop_mapping.txt | wc -l

# Now get the actual 1000G population (sub-pop level) from the KG_reference_final FAM
echo "=== KG_FAM_HEAD ==="
head -5 /staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/1000g_grch38/KG_reference_final.fam
echo "=== KG_FAM_POPS ==="
awk '{print $1}' /staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/1000g_grch38/KG_reference_final.fam | sort | uniq -c | sort -rn
echo "=== KG_FAM_TOTAL ==="
wc -l /staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/1000g_grch38/KG_reference_final.fam

echo "=== END ==="
