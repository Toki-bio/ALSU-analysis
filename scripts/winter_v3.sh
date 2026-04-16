#!/bin/bash
set -eo pipefail
TSV=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz/UZB_all.INFO_metrics.tsv

echo "=== Winter 2025 R2 from UZB_all.INFO_metrics.tsv ==="
echo "Column 8 = R2 (confirmed from header)"

# Use the pre-extracted TSV — skip header, R2 is column 8
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
}
END {
  printf "Total variants: %d\n", n
  printf "Mean R2: %.6f\n", s/n
  printf "\nHistogram (10 bins):\n"
  for (i=0; i<10; i++)
    printf "  %.1f-%.1f: %12d  (%5.1f%%)\n", i/10, (i+1)/10, b[i], b[i]/n*100
  printf "\nThresholds:\n"
  printf "  R2 >= 0.30: %12d  (%5.1f%%)\n", a, a/n*100
  printf "  R2 >= 0.50: %12d  (%5.1f%%)\n", e, e/n*100
  printf "  R2 >= 0.80: %12d  (%5.1f%%)\n", c, c/n*100
  printf "  R2 >= 0.90: %12d  (%5.1f%%)\n", d, d/n*100
}
'

echo ""
echo "=== Typed vs Imputed (columns 10,11) ==="
tail -n +2 "$TSV" | awk -F'\t' '
{
  if ($10 == 1) typed++
  else if ($11 == 1) imputed++
  else other++
  n++
}
END {
  printf "Typed: %d\n", typed
  printf "Imputed: %d\n", imputed
  if (other > 0) printf "Other: %d\n", other
  printf "Total: %d\n", n
}
'

echo ""
echo "=== Per-chromosome counts ==="
tail -n +2 "$TSV" | awk -F'\t' '{print $1}' | sort -V | uniq -c | awk '{printf "%-6s %d\n", $2, $1}'

echo ""
echo "=== Spring download status ==="
stat -c '%s' /tmp/michigan_results/chr_1.zip 2>/dev/null || echo none
ls /tmp/michigan_results/chr_*.zip 2>/dev/null | wc -l
tail -5 /tmp/michigan_analyze.log 2>/dev/null

echo "=== DONE ==="
