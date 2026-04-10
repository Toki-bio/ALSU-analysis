#!/bin/bash
cd /staging/ALSU-analysis/winter2025/PLINK_301125_0312/

echo "=== EXTRACT ACTUAL 99 SAMPLES (CORRECT F_MISS > 0.20) ==="
awk 'NR>1 && $6+0 > 0.20 {print $1"\t"$2}' ConvSK_raw_miss.imiss | sort > remove_correct_99.txt
echo "Correct 99 samples:"
wc -l remove_correct_99.txt

echo ""
echo "=== EXTRACT THE 92 SAMPLES FROM CURRENT remove_miss20.txt ==="
sort remove_miss20.txt > remove_current_92_sorted.txt
wc -l remove_current_92_sorted.txt

echo ""
echo "=== FIND THE 7 MISSING SAMPLES ==="
echo "Samples in correct list but NOT in current 92-sample list:"
comm -23 remove_correct_99.txt remove_current_92_sorted.txt
echo ""
echo "These 7 are the samples that should have been removed but weren't"
