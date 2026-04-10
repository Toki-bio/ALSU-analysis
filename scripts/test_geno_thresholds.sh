#!/bin/bash
# Test multiple --geno thresholds on winter and spring data
# Run on Biotech2024

set -e

echo "=== WINTER 2025 DATA ==="
WDIR="/staging/ALSU-analysis/winter2025/PLINK_301125_0312"
WTEST="/tmp/geno_test_winter"
mkdir -p "$WTEST"
cd "$WDIR"

# Confirm starting point
echo "Winter input:"
wc -l ConvSK_mind20_dedup.bim | awk '{print "  Variants in bim:", $1}'
wc -l ConvSK_mind20_dedup.fam | awk '{print "  Samples in fam:", $1}'

# Also show genotyping rate from a quick check
plink --bfile ConvSK_mind20_dedup --chr 1-22 --make-bed --out "${WTEST}/w_chr" --threads 4 2>/dev/null
WAUTO=$(wc -l < "${WTEST}/w_chr.bim")
echo "  Autosomal variants: $WAUTO"
grep "Total genotyping rate" "${WTEST}/w_chr.log"

# Test each threshold - chr filter already applied via w_chr
for G in 0.005 0.01 0.02 0.03 0.05 0.10; do
  plink --bfile "${WTEST}/w_chr" \
    --geno $G --maf 0.01 --hwe 1e-6 \
    --make-bed --out "${WTEST}/w_geno_${G}" --threads 4 2>/dev/null
  N=$(wc -l < "${WTEST}/w_geno_${G}.bim")
  # Extract per-filter removal counts from log
  GREM=$(grep "variants removed due to missing genotype" "${WTEST}/w_geno_${G}.log" | awk '{print $1}')
  HREM=$(grep "variants removed due to Hardy-Weinberg" "${WTEST}/w_geno_${G}.log" | awk '{print $1}' | head -1)
  MREM=$(grep "variants removed due to minor allele" "${WTEST}/w_geno_${G}.log" | awk '{print $1}')
  echo "  --geno $G: $N variants kept (geno_rm=$GREM hwe_rm=$HREM maf_rm=$MREM)"
done

# Also test with NO --geno filter at all
plink --bfile "${WTEST}/w_chr" \
  --maf 0.01 --hwe 1e-6 \
  --make-bed --out "${WTEST}/w_geno_none" --threads 4 2>/dev/null
N=$(wc -l < "${WTEST}/w_geno_none.bim")
echo "  --geno NONE: $N variants kept"

# Check: does 473081 appear anywhere?
echo ""
echo "  Winter snpqc file (actual):"
wc -l "${WDIR}/ConvSK_mind20_dedup_snpqc.bim" | awk '{print "  Variants:", $1}'

echo ""
echo "=== SPRING 2026 DATA ==="
SDIR="/staging/ALSU-analysis/spring2026"
STEST="/tmp/geno_test_spring"
mkdir -p "$STEST"
cd "$SDIR"

echo "Spring input:"
wc -l ConvSK_mind20_dedup.bim | awk '{print "  Variants in bim:", $1}'
wc -l ConvSK_mind20_dedup.fam | awk '{print "  Samples in fam:", $1}'

plink --bfile ConvSK_mind20_dedup --chr 1-22 --make-bed --out "${STEST}/s_chr" --threads 4 2>/dev/null
SAUTO=$(wc -l < "${STEST}/s_chr.bim")
echo "  Autosomal variants: $SAUTO"
grep "Total genotyping rate" "${STEST}/s_chr.log"

for G in 0.005 0.01 0.02 0.03 0.05 0.10; do
  plink --bfile "${STEST}/s_chr" \
    --geno $G --maf 0.01 --hwe 1e-6 \
    --make-bed --out "${STEST}/s_geno_${G}" --threads 4 2>/dev/null
  N=$(wc -l < "${STEST}/s_geno_${G}.bim")
  GREM=$(grep "variants removed due to missing genotype" "${STEST}/s_geno_${G}.log" | awk '{print $1}')
  HREM=$(grep "variants removed due to Hardy-Weinberg" "${STEST}/s_geno_${G}.log" | awk '{print $1}' | head -1)
  MREM=$(grep "variants removed due to minor allele" "${STEST}/s_geno_${G}.log" | awk '{print $1}')
  echo "  --geno $G: $N variants kept (geno_rm=$GREM hwe_rm=$HREM maf_rm=$MREM)"
done

plink --bfile "${STEST}/s_chr" \
  --maf 0.01 --hwe 1e-6 \
  --make-bed --out "${STEST}/s_geno_none" --threads 4 2>/dev/null
N=$(wc -l < "${STEST}/s_geno_none.bim")
echo "  --geno NONE: $N variants kept"

echo ""
echo "=== DONE ==="
