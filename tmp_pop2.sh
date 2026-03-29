#!/bin/bash
echo "=== PCA_POPMAPPING ==="
cat ~/v2/pca/pop_mapping.txt | head -20
echo "=== PCA_POPMAPPING_COUNTS ==="
cut -f2 ~/v2/pca/pop_mapping.txt | sort | uniq -c | sort -rn
echo "=== PCA_POPMAPPING_TOTAL ==="
wc -l ~/v2/pca/pop_mapping.txt
# Also check reanalysis log for pop info
echo "=== REANALYSIS_PHASES ==="
grep -n 'Phase\|population\|pop\|1000G\|merge\|reference\|ref_samples' ~/v2/reanalysis.log | head -30
echo "=== REANALYSIS_PHASEC ==="
grep -A5 'Phase C\|global' ~/v2/reanalysis.log | head -40
echo "=== END ==="
