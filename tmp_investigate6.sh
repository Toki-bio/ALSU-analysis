#!/bin/bash
UNZ=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz
FC=$UNZ/filtered_clean
S7=$UNZ/step7_ld_pruning

echo "=== 1. Are info.gz and dose.vcf.gz the same file? ==="
ls -li $UNZ/chr22.info.gz $UNZ/chr22.dose.vcf.gz 2>/dev/null

echo ""
echo "=== 2. Quick TYPED count in chr22 dose VCF (grep only) ==="
typed_chr22=$(zcat $UNZ/chr22.dose.vcf.gz 2>/dev/null | grep -c 'TYPED')
imputed_chr22=$(zcat $UNZ/chr22.dose.vcf.gz 2>/dev/null | grep -v '^#' | wc -l)
echo "chr22 TYPED (contains TYPED): $typed_chr22"
echo "chr22 total data lines: $imputed_chr22"

echo ""
echo "=== 3. Sample TYPED-only variant lines from chr22 ==="
zcat $UNZ/chr22.dose.vcf.gz 2>/dev/null | grep -v '^#' | awk -F'\t' '$8 ~ /TYPED/ && $8 !~ /IMPUTED/' | head -5

echo ""
echo "=== 4. Sample TYPED+IMPUTED variant lines ==="
zcat $UNZ/chr22.dose.vcf.gz 2>/dev/null | grep -v '^#' | awk -F'\t' '$8 ~ /TYPED/ && $8 ~ /IMPUTED/' | head -5

echo ""
echo "=== 5. filtered_clean dir listing ==="
ls -la $FC/ 2>/dev/null | head -50

echo ""
echo "=== 6. Position BED ==="
wc -l $FC/chr*_positions.bed 2>/dev/null

echo ""
echo "=== 7. Ref log ==="
cat $FC/chr10_ref.log 2>/dev/null
echo "---"
cat $FC/chr1_ref.log 2>/dev/null

echo ""
echo "=== 8. Position BED content ==="
head -3 $FC/chr10_positions.bed 2>/dev/null

echo ""
echo "=== 9. Scripts in the pipeline ==="
find /staging/ALSU-analysis/winter2025 -maxdepth 3 -name "*.sh" -o -name "*.py" 2>/dev/null

echo ""
echo "=== 10. Bash history ==="
tail -100 ~/.bash_history 2>/dev/null

echo ""
echo "=== 11. BIM: chr22 typed-like variants ==="
echo "BIM chr22 total:"
awk '$1==22' $S7/UZB_imputed_HQ_clean.bim | wc -l
echo "BIM chr22 dot IDs:"
awk '$1==22 && $2=="."' $S7/UZB_imputed_HQ_clean.bim | wc -l
echo "BIM chr22 sample dots:"
awk '$1==22 && $2=="."' $S7/UZB_imputed_HQ_clean.bim | head -5

echo ""
echo "=== 12. Admixture FAM counts ==="
for f in $(find /staging/ALSU-analysis/winter2025/admixture_analysis -name "*.fam" 2>/dev/null); do
  echo "$(basename $f): $(wc -l < $f)"
done

echo ""
echo "=== 13. Check for README or notes ==="
find /staging/ALSU-analysis/winter2025 -maxdepth 2 -name "README*" -o -name "NOTES*" -o -name "*.md" -o -name "*.txt" 2>/dev/null

echo ""
echo "=== 14. File listing in unz root ==="
ls -la $UNZ/*.sh $UNZ/*.py $UNZ/*.txt $UNZ/*.log 2>/dev/null

echo ""
echo "=== DONE ==="
