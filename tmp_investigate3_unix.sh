#!/bin/bash
# Investigation 3: The 837K extra variants â€” are they TYPED (non-imputed) SNPs?
UNZ=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz

echo "=== 1. Check UZB_all.HQ.query.chrposrefalt.tsv ==="
if [ -f "$UNZ/UZB_all.HQ.query.chrposrefalt.tsv" ]; then
  echo "File exists"
  echo "Lines: $(wc -l < "$UNZ/UZB_all.HQ.query.chrposrefalt.tsv")"
  echo "First 10 lines:"
  head -10 "$UNZ/UZB_all.HQ.query.chrposrefalt.tsv"
fi

echo ""
echo "=== 2. Check UZB_all.HQ_imputed.R2ge0p8.MAFge0p001.tsv columns ==="
echo "Lines: $(wc -l < "$UNZ/UZB_all.HQ_imputed.R2ge0p8.MAFge0p001.tsv")"
echo "Header:"
head -1 "$UNZ/UZB_all.HQ_imputed.R2ge0p8.MAFge0p001.tsv"
echo "First 5 data:"
head -6 "$UNZ/UZB_all.HQ_imputed.R2ge0p8.MAFge0p001.tsv" | tail -5

echo ""
echo "=== 3. Check the FULL frequency table for TYPED vs IMPUTED column ==="
echo "Header:"
head -1 "$UNZ/UZB_all.AF_DS_R2.tsv"
echo ""
echo "Columns (numbered):"
head -1 "$UNZ/UZB_all.AF_DS_R2.tsv" | awk -F'\t' '{for(i=1;i<=NF;i++) print i": "$i}'

echo ""
echo "=== 4. Check for info files from Michigan that have Genotyped/Typed column ==="
find "$UNZ" -maxdepth 1 -name '*.info*' -o -name '*.info.gz' 2>/dev/null | head -10
find "$UNZ/.." -maxdepth 1 -name '*.info*' -o -name '*.info.gz' 2>/dev/null | head -10

echo ""
echo "=== 5. List files in imputation_results/unz ==="
ls -la "$UNZ/" | head -40

echo ""
echo "=== 6. List files in imputation_results/ ==="
ls -la "$UNZ/../" | head -20

echo ""
echo "=== 7. Check for chr1.info.gz or similar ==="
find /staging/ALSU-analysis/winter2025 -maxdepth 5 -name '*.info*' 2>/dev/null | head -10

echo ""
echo "=== 8. Check dose VCF for ID field (first 20 variants of chr22) ==="
if [ -f "$UNZ/chr22.dose.vcf.gz" ]; then
  zcat "$UNZ/chr22.dose.vcf.gz" 2>/dev/null | grep -v '^#' | head -20 | cut -f1-5
else
  echo "chr22.dose.vcf.gz not found in unz"
  # Check parent dir
  PARENT="$UNZ/.."
  find "$PARENT" -maxdepth 2 -name 'chr22*dose*' 2>/dev/null | head -5
fi

echo ""
echo "=== 9. Check filtered_clean directory ==="
ls -la "$UNZ/filtered_clean/" 2>/dev/null | head -30

echo ""
echo "=== 10. Check the .bim file for variant count and ID patterns ==="
BIM="/staging/ALSU-analysis/winter2025/PLINK_301125_0312/step7_ld_pruning/UZB_imputed_HQ_clean.bim"
if [ -f "$BIM" ]; then
  echo "Total variants: $(wc -l < "$BIM")"
  echo "Variants with rsID: $(grep -c 'rs' "$BIM")"
  echo "Variants with dot ID: $(awk '$2=="." {c++} END {print c+0}' "$BIM")"
  echo "Variants with chr:pos format: $(grep -c ':' "$BIM" | head -1)"
  echo "First 10 lines:"
  head -10 "$BIM"
fi

echo ""
echo "=== 11. Check the UZB_all.AF_DS_Rsq.tsv â€” another freq table? ==="
if [ -f "$UNZ/UZB_all.AF_DS_Rsq.tsv" ]; then
  echo "Lines: $(wc -l < "$UNZ/UZB_all.AF_DS_Rsq.tsv")"
  echo "Header:"
  head -1 "$UNZ/UZB_all.AF_DS_Rsq.tsv"
  echo "Columns:"
  head -1 "$UNZ/UZB_all.AF_DS_Rsq.tsv" | awk -F'\t' '{for(i=1;i<=NF;i++) print i": "$i}'
fi

echo ""
echo "=== 12. Count distinct IMPUTED column values in AF_DS_R2 ==="
# Check if there's a column for typed/imputed status
# In Michigan imputation info files, column is typically labeled TYPED/IMPUTED/Genotyped
# The filter was IMPUTED==1. Let's see what column 10 looks like
echo "Column 10 sample:"
awk -F'\t' 'NR<=6 {print NR": "$10}' "$UNZ/UZB_all.AF_DS_R2.tsv"

echo ""
echo "=== 13. Count typed vs imputed in full freq table (first 100 lines) ==="
awk -F'\t' 'NR>1 && NR<=101 {print $10}' "$UNZ/UZB_all.AF_DS_R2.tsv" | sort | uniq -c

echo ""
echo "=== DONE ==="
