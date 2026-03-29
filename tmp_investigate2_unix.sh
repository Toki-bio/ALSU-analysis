#!/bin/bash
# Deep investigation of variant filtering - focus on winter2025 and SNPfreq directories

echo "=== 1. winter2025 directory structure ==="
find /staging/ALSU-analysis/winter2025 -maxdepth 3 -type d 2>/dev/null

echo ""
echo "=== 2. Key files in winter2025 (variant lists, freq tables, logs) ==="
find /staging/ALSU-analysis/winter2025 -maxdepth 4 -type f \( -name '*.txt' -o -name '*.tsv' -o -name '*.log' -o -name '*.sh' -o -name '*.py' \) 2>/dev/null | head -50

echo ""
echo "=== 3. SNPfreq directory ==="
ls -la /staging/ALSU-analysis/SNPfreq/ 2>/dev/null

echo ""
echo "=== 4. admixture_analysis directory ==="
ls -la /staging/ALSU-analysis/admixture_analysis/ 2>/dev/null

echo ""
echo "=== 5. Look for frequency table with R2 column ==="
find /staging/ALSU-analysis -maxdepth 6 -type f -name '*R2*' -o -name '*AF_DS*' -o -name '*freq*' -o -name '*HQ*' 2>/dev/null | head -30

echo ""
echo "=== 6. Look for imputation result VCFs ==="
find /staging/ALSU-analysis/winter2025 -maxdepth 6 -type f -name '*.vcf.gz' 2>/dev/null | head -20

echo ""
echo "=== 7. Look for filtered_clean directory ==="
find /staging/ALSU-analysis -maxdepth 6 -type d -name 'filtered_clean' 2>/dev/null
find /staging/ALSU-analysis -maxdepth 6 -type d -name 'unz' 2>/dev/null

echo ""
echo "=== 8. Check for the HQ variant ID file and count dots ==="
HQ_FILE=$(find /staging/ALSU-analysis -maxdepth 6 -name 'HQ_variant_ids.txt' 2>/dev/null | head -1)
if [ -n "$HQ_FILE" ]; then
  echo "Found: $HQ_FILE"
  echo "Total lines: $(wc -l < "$HQ_FILE")"
  echo "Lines with dot only: $(grep -c '^\.$' "$HQ_FILE")"
  echo "Lines with rsID: $(grep -c '^rs' "$HQ_FILE")"
  echo "First 10 lines:"
  head -10 "$HQ_FILE"
else
  echo "HQ_variant_ids.txt not found"
fi

echo ""
echo "=== 9. Look for the frequency table ==="
FREQ_FILE=$(find /staging/ALSU-analysis -maxdepth 6 -name 'UZB_all*R2*' 2>/dev/null | head -1)
if [ -n "$FREQ_FILE" ]; then
  echo "Found: $FREQ_FILE"
  echo "Size: $(ls -lh "$FREQ_FILE" | awk '{print $5}')"
  echo "Lines: $(wc -l < "$FREQ_FILE")"
  echo "Header:"
  head -1 "$FREQ_FILE"
  echo "First 5 data lines:"
  head -6 "$FREQ_FILE" | tail -5
else
  echo "Frequency table not found"
fi

echo ""
echo "=== 10. Look for UZB_all.AF_DS_R2.tsv (original unfiltered freq table) ==="
ORIG_FILE=$(find /staging/ALSU-analysis -maxdepth 6 -name 'UZB_all.AF_DS_R2*' 2>/dev/null | head -1)
if [ -n "$ORIG_FILE" ]; then
  echo "Found: $ORIG_FILE"
  echo "Size: $(ls -lh "$ORIG_FILE" | awk '{print $5}')"
  echo "Total lines: $(wc -l < "$ORIG_FILE")"
  echo "Header:"
  head -1 "$ORIG_FILE"
  echo ""
  echo "R2 distribution (sample of 1M lines):"
  awk -F'\t' 'NR>1 && NR<=1000001 {
    r2 = $8+0
    if (r2 >= 0.9) a++
    else if (r2 >= 0.8) b++
    else if (r2 >= 0.3) c++
    else d++
  } END {
    print "R2>=0.9: " a
    print "R2 0.8-0.9: " b
    print "R2 0.3-0.8: " c
    print "R2<0.3: " d
  }' "$ORIG_FILE"
else
  echo "UZB_all.AF_DS_R2.tsv not found"
fi

echo ""
echo "=== 11. Check PLINK bed/bim/fam in filtered_clean ==="
find /staging/ALSU-analysis -maxdepth 6 -name 'UZB_imputed_HQ_clean*' -type f 2>/dev/null | head -20

echo ""
echo "=== 12. ALSU.vcf header check ==="
if [ -f /staging/ALSU-analysis/ALSU.vcf ]; then
  echo "ALSU.vcf exists, size: $(ls -lh /staging/ALSU-analysis/ALSU.vcf | awk '{print $5}')"
  echo "Sample count: $(head -1000 /staging/ALSU-analysis/ALSU.vcf | grep '^#CHROM' | awk '{print NF-9}')"
  echo "Total variant lines (approx from first 100K lines):"
  head -100000 /staging/ALSU-analysis/ALSU.vcf | grep -cv '^#'
fi

echo ""
echo "=== DONE ==="
