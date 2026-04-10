#!/bin/bash
cd /staging/ALSU-analysis/winter2025/PLINK_301125_0312/

echo "=== THE 7 MISSED SAMPLES (samples with F_MISS > 0.20 not in original 92-sample list) ==="
echo ""

# Create temp files in /tmp which should be writable
awk 'NR>1 && $6+0 > 0.20 {print $1"\t"$2}' ConvSK_raw_miss.imiss | sort > /tmp/correct_99.txt 2>/dev/null || echo "Error: cannot create temp file"
sort remove_miss20.txt > /tmp/current_92.txt 2>/dev/null || echo "Error: cannot sort remove_miss20.txt"

# If files were created, show the diff
if [ -f /tmp/correct_99.txt ] && [ -f /tmp/current_92.txt ]; then
  echo "Samples in CORRECT list (99) but NOT in buggy list (92):"
  echo "---"
  comm -23 /tmp/correct_99.txt /tmp/current_92.txt
  echo "---"
  echo ""
  
  # Count how many
  missed=$(comm -23 /tmp/correct_99.txt /tmp/current_92.txt | wc -l)
  echo "Total missed: $missed samples"
else
  # Alternative: just show first few from correct and see what's missing
  echo "Correct 99 (first 10):"
  awk 'NR>1 && $6+0 > 0.20 {print $1"\t"$2}' ConvSK_raw_miss.imiss | head -10
  echo ""
  echo "Buggy 92 (first 10):" 
  head -10 remove_miss20.txt
fi
