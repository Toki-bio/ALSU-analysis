#!/bin/bash
set -eo pipefail
TSV=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz/UZB_all.INFO_metrics.tsv
OUT=/tmp/winter_stats.txt

echo "Processing $(wc -l < "$TSV") lines..." > "$OUT"

tail -n +2 "$TSV" | awk -F'\t' '
BEGIN { for (i=0; i<10; i++) b[i] = 0 }
{
  v = $8 + 0
  n++; s += v
  bin = int(v * 10)
  if (bin >= 10) bin = 9
  b[bin]++
  if (v >= 0.3) a++
  if (v >= 0.5) e++
  if (v >= 0.8) c++
  if (v >= 0.9) d++
  # typed/imputed
  if ($10 == 1) typed++
  else imputed++
}
END {
  printf "WINTER_2025_R2_DISTRIBUTION\n"
  printf "Total_variants: %d\n", n
  printf "Mean_R2: %.6f\n", s/n
  printf "Typed: %d\n", typed
  printf "Imputed: %d\n", imputed
  printf "\nHISTOGRAM\n"
  for (i=0; i<10; i++)
    printf "BIN %.1f-%.1f %d %.3f\n", i/10, (i+1)/10, b[i], b[i]/n*100
  printf "\nTHRESHOLDS\n"
  printf "R2_ge_030 %d %.3f\n", a, a/n*100
  printf "R2_ge_050 %d %.3f\n", e, e/n*100
  printf "R2_ge_080 %d %.3f\n", c, c/n*100
  printf "R2_ge_090 %d %.3f\n", d, d/n*100
}
' >> "$OUT"

# Per-chr counts
echo "" >> "$OUT"
echo "PER_CHR" >> "$OUT"
tail -n +2 "$TSV" | awk -F'\t' '{print $1}' | sort -V | uniq -c | awk '{printf "%s %d\n", $2, $1}' >> "$OUT"

echo "DONE" >> "$OUT"
