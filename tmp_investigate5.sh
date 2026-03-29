#!/bin/bash
UNZ=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz
FC=$UNZ/filtered_clean
S7=$UNZ/step7_ld_pruning

echo "=== 1. Dose VCF INFO flags: TYPED vs IMPUTED counts across all chromosomes ==="
for chr in $(seq 1 22); do
  f=$UNZ/chr${chr}.dose.vcf.gz
  if [ -f "$f" ]; then
    counts=$(zcat "$f" 2>/dev/null | grep -v '^#' | awk -F'\t' '{
      info=$8
      typed=0; imputed=0
      if(info ~ /TYPED/) typed=1
      if(info ~ /IMPUTED/) imputed=1
      if(typed && imputed) both++
      else if(typed) tonly++
      else if(imputed) ionly++
    } END {print tonly+0, ionly+0, both+0}')
    echo "chr$chr: TypedOnly=$(echo $counts|awk '{print $1}') ImputedOnly=$(echo $counts|awk '{print $2}') Both=$(echo $counts|awk '{print $3}')"
  fi
done

echo ""
echo "=== 2. Dose VCF: how many '.' ID variants per chr? ==="
for chr in 1 10 22; do
  f=$UNZ/chr${chr}.dose.vcf.gz
  total=$(zcat "$f" 2>/dev/null | grep -cv '^#')
  dots=$(zcat "$f" 2>/dev/null | grep -v '^#' | awk -F'\t' '$3=="."' | wc -l)
  echo "chr$chr: total=$total dots=$dots"
done

echo ""
echo "=== 3. Check filtered_clean contents ==="
ls -la $FC/ 2>/dev/null | head -30
echo "---"
echo "Total files:"
ls $FC/ 2>/dev/null | wc -l

echo ""
echo "=== 4. Position BED file (chr10 sample) ==="
wc -l $FC/chr10_positions.bed 2>/dev/null
head -5 $FC/chr10_positions.bed 2>/dev/null

echo ""
echo "=== 5. Ref log file ==="
cat $FC/chr10_ref.log 2>/dev/null

echo ""
echo "=== 6. Check ALL position BED file totals ==="
total_pos=0
for chr in $(seq 1 22); do
  if [ -f "$FC/chr${chr}_positions.bed" ]; then
    n=$(wc -l < "$FC/chr${chr}_positions.bed")
    total_pos=$((total_pos + n))
    echo "chr$chr: $n"
  fi
done
echo "TOTAL positions in BED files: $total_pos"

echo ""
echo "=== 7. Compare: BIM variants vs BED positions ==="
echo "BIM total: $(wc -l < $S7/UZB_imputed_HQ_clean.bim)"
echo "BIM dot-IDs: $(awk '$2=="."' $S7/UZB_imputed_HQ_clean.bim | wc -l)"

echo ""
echo "=== 8. Bash history (last 50 lines) ==="
tail -50 ~/.bash_history 2>/dev/null

echo ""
echo "=== 9. Check for any .sh or .py scripts in the pipeline dirs ==="
find $UNZ -maxdepth 2 -name "*.sh" -o -name "*.py" 2>/dev/null
find /staging/ALSU-analysis/winter2025 -maxdepth 2 -name "*.sh" -o -name "*.py" 2>/dev/null

echo ""
echo "=== 10. ADMIXTURE dir: FAM file sample counts ==="
for f in /staging/ALSU-analysis/winter2025/admixture_analysis/*.fam /staging/ALSU-analysis/winter2025/admixture_analysis/**/*.fam; do
  if [ -f "$f" ]; then
    echo "$(basename $f): $(wc -l < $f) samples"
  fi
done

echo ""
echo "=== DONE ==="
