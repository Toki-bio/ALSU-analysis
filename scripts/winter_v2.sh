#!/bin/bash
set -eo pipefail
DIR=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz

echo "=== UZB_all.INFO_metrics.tsv header ==="
head -3 $DIR/UZB_all.INFO_metrics.tsv

echo ""
echo "=== UZB_all.INFO_metrics.tsv line count ==="
wc -l $DIR/UZB_all.INFO_metrics.tsv

echo ""
echo "=== UZB_all.AF_DS_R2.tsv header ==="
head -3 $DIR/UZB_all.AF_DS_R2.tsv

echo ""
echo "=== Winter samples from chr1.dose.vcf.gz ==="
bcftools query -l $DIR/chr1.dose.vcf.gz 2>/dev/null | wc -l

echo ""
echo "=== Winter R2 from INFO field of .info.gz VCFs ==="
echo "Extracting R2 from chr1.info.gz INFO field..."
zcat $DIR/chr1.info.gz | grep -v '^#' | head -1 | cut -f8

echo ""
echo "=== Extract R2 distribution from all 22 chr .info.gz ==="
for chr in $(seq 1 22); do
  f=$DIR/chr${chr}.info.gz
  if [ -f "$f" ]; then
    zcat "$f" | grep -v '^#' | sed 's/.*R2=\([0-9.]*\).*/\1/'
  fi
done | awk '
BEGIN { for (i=0; i<10; i++) b[i] = 0 }
{
  v = $1 + 0
  n++; s += v
  bin = int(v * 10)
  if (bin >= 10) bin = 9
  b[bin]++
  if (v >= 0.3) a++
  if (v >= 0.8) c++
  if (v >= 0.9) d++
}
END {
  printf "Total: %d\n", n
  printf "Mean R2: %.6f\n", s/n
  printf "\nHistogram:\n"
  for (i=0; i<10; i++)
    printf "  %.1f-%.1f: %12d  (%5.1f%%)\n", i/10, (i+1)/10, b[i], b[i]/n*100
  printf "\nThresholds:\n"
  printf "  >=0.30: %12d  (%5.1f%%)\n", a, a/n*100
  printf "  >=0.80: %12d  (%5.1f%%)\n", c, c/n*100
  printf "  >=0.90: %12d  (%5.1f%%)\n", d, d/n*100
}
'

echo ""
echo "=== Per-chr variant counts (from .info.gz) ==="
for chr in $(seq 1 22); do
  f=$DIR/chr${chr}.info.gz
  if [ -f "$f" ]; then
    n=$(zcat "$f" | grep -v '^#' | wc -l)
    printf "chr%-2s: %d\n" "$chr" "$n"
  fi
done

echo ""
echo "=== Typed vs Imputed (from INFO field) ==="
for chr in $(seq 1 22); do
  f=$DIR/chr${chr}.info.gz
  if [ -f "$f" ]; then
    zcat "$f" | grep -v '^#' | grep -oP '(TYPED|IMPUTED)' 
  fi
done | sort | uniq -c | sort -rn

echo ""
echo "=== Spring download check ==="
stat -c '%s' /tmp/michigan_results/chr_1.zip 2>/dev/null || echo none
ls /tmp/michigan_results/chr_*.zip 2>/dev/null | wc -l
tail -3 /tmp/michigan_analyze.log 2>/dev/null

echo "=== DONE ==="
