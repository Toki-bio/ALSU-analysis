#!/bin/bash
set -e
WINTER=/staging/ALSU-analysis/winter2025/PLINK_301125_0312

echo "=== 1. BASH HISTORY ==="
grep -i 'snpqc\|mind20_dedup.*make-bed\|--geno.*--maf.*--hwe' ~/.bash_history 2>/dev/null || echo "(no matches in root history)"
grep -i 'snpqc\|--geno.*--maf' /home/bioinfo/.bash_history 2>/dev/null || echo "(no bioinfo history)"
grep -i 'snpqc\|--geno.*--maf' /home/copilot/.bash_history 2>/dev/null || echo "(no copilot history)"

echo ""
echo "=== 2. OVERWRITTEN LOG (first 40 lines) ==="
head -40 ${WINTER}/ConvSK_mind20_dedup_snpqc.log 2>/dev/null || echo "(log not found)"

echo ""
echo "=== 3. REVERSE-ENGINEER FROM DATA ==="
# Check max missingness among the 473,081 retained variants
plink --bfile ${WINTER}/ConvSK_mind20_dedup \
  --extract <(awk '{print $2}' ${WINTER}/ConvSK_mind20_dedup_snpqc.bim) \
  --missing --out /tmp/winter_retained_check --noweb 2>&1 | tail -5

echo ""
echo "--- Max per-SNP missingness among retained ---"
awk 'NR>1 {if ($5+0 > max) max=$5} END {printf "Max lmiss: %.6f\n", max}' \
  /tmp/winter_retained_check.lmiss

echo ""
echo "--- Top 10 highest-missingness retained SNPs ---"
sort -k5 -g -r /tmp/winter_retained_check.lmiss | head -11

echo ""
# Check min MAF among retained
plink --bfile ${WINTER}/ConvSK_mind20_dedup \
  --extract <(awk '{print $2}' ${WINTER}/ConvSK_mind20_dedup_snpqc.bim) \
  --freq --out /tmp/winter_retained_check --noweb 2>&1 | tail -5

echo ""
echo "--- Min MAF among retained ---"
awk 'NR>1 && $5!="NA" {if (min=="" || $5+0 < min) min=$5} END {printf "Min MAF: %.6f\n", min}' \
  /tmp/winter_retained_check.frq

echo ""
echo "--- Bottom 10 lowest-MAF retained SNPs ---"
sort -k5 -g /tmp/winter_retained_check.frq | head -11

echo ""
echo "--- HWE check: were midp used? Check min HWE p among retained ---"
plink --bfile ${WINTER}/ConvSK_mind20_dedup \
  --extract <(awk '{print $2}' ${WINTER}/ConvSK_mind20_dedup_snpqc.bim) \
  --hardy --out /tmp/winter_retained_check --noweb 2>&1 | tail -3
awk 'NR>1 && $9!="NA" {if (min=="" || $9+0 < min) min=$9} END {printf "Min HWE p: %g\n", min}' \
  /tmp/winter_retained_check.hwe

echo ""
echo "=== DONE ==="
