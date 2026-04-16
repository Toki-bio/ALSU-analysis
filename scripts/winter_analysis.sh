#!/bin/bash
set -eo pipefail
INFODIR=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz

echo "=== Winter 2025 .info.gz analysis ==="
echo "Files found: $(ls $INFODIR/*.info.gz 2>/dev/null | wc -l)"

echo "---HEADER---"
zcat $INFODIR/chr1.info.gz | head -2

echo ""
echo "---WINTER SAMPLES---"
bcftools query -l $INFODIR/chr1.dose.vcf.gz 2>/dev/null | wc -l

echo ""
echo "---WINTER PER-CHR---"
total=0
for chr in $(seq 1 22); do
  f=$INFODIR/chr${chr}.info.gz
  if [ -f "$f" ]; then
    n=$(zcat "$f" | tail -n +2 | wc -l)
    total=$((total + n))
    printf "chr%-2s: %d\n" "$chr" "$n"
  else
    echo "chr${chr}: MISSING"
  fi
done
echo "TOTAL: $total"

echo ""
echo "---WINTER R2 DISTRIBUTION---"
for chr in $(seq 1 22); do
  f=$INFODIR/chr${chr}.info.gz
  [ -f "$f" ] && zcat "$f" | tail -n +2 | awk '{print $7}'
done | awk '
BEGIN { for (i=0; i<10; i++) b[i] = 0 }
{
  n++; s += $1
  bin = int($1 * 10)
  if (bin >= 10) bin = 9
  b[bin]++
  if ($1 >= 0.3) a++
  if ($1 >= 0.8) c++
  if ($1 >= 0.9) d++
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
echo "---WINTER TYPED vs IMPUTED---"
for chr in $(seq 1 22); do
  f=$INFODIR/chr${chr}.info.gz
  [ -f "$f" ] && zcat "$f" | tail -n +2 | awk '{print $8}'
done | sort | uniq -c | sort -rn

echo ""
echo "---SPRING DOWNLOAD STATUS---"
stat -c '%s' /tmp/michigan_results/chr_1.zip 2>/dev/null
ls /tmp/michigan_results/chr_*.zip 2>/dev/null | wc -l
ps aux | grep 'curl.*imputationserver' | grep -v grep | wc -l
du -sh /tmp/michigan_results/ 2>/dev/null

echo ""
echo "=== DONE ==="
