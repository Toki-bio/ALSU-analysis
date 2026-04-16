#!/bin/bash
# Post-imputation pipeline: spring 2026
# Combines Step 5 (ID normalization) + Step 6 (Final QC)
# Run on DRAGEN: /staging/ALSU-analysis/spring2026/
set -euo pipefail

WORKDIR=/staging/ALSU-analysis/spring2026
IMPDIR=${WORKDIR}/imputation
POSTDIR=${WORKDIR}/post_imputation
mkdir -p "$POSTDIR"

cd "$WORKDIR"

echo "========================================="
echo "STEP 5: Sample ID Check + Normalization"
echo "========================================="

# Step 3e already fixed Cyrillic IDs BEFORE Michigan upload.
# Michigan v2 may or may not add numeric prefixes to sample IDs.
# Check and strip if present.

echo "=== Checking Michigan sample ID format ==="
# Save sample list to file first (avoids SIGPIPE with head + pipefail)
bcftools query -l ${IMPDIR}/chr1.dose.vcf.gz > ${POSTDIR}/sample_list.txt
NSAMPLES=$(wc -l < ${POSTDIR}/sample_list.txt)
echo "Total samples in imputed VCF: $NSAMPLES"
head -5 ${POSTDIR}/sample_list.txt

NPFX=$(grep -cP '^\d+_' ${POSTDIR}/sample_list.txt || true)
echo "Samples with numeric prefix: $NPFX"

NASCII=$(grep -cP '[^\x00-\x7F]' ${POSTDIR}/sample_list.txt || true)
echo "Non-ASCII sample IDs: $NASCII"

# Michigan v2 quirk: chr1-9 have numeric prefix (1_sampleID), chr10-22 don't.
# Must check each chromosome individually.
if [ "$NPFX" -gt 0 ]; then
  echo "=== Stripping Michigan numeric prefix (per-chromosome) ==="
  # Create rename mapping from chr1's prefixed IDs
  awk -F'_' 'NF>1 && $1~/^[0-9]+$/ {
    old=$0; sub(/^[0-9]+_/,""); print old"\t"$0
  }' ${POSTDIR}/sample_list.txt > ${POSTDIR}/strip_prefix.tsv
  echo "Prefix mappings: $(wc -l < ${POSTDIR}/strip_prefix.tsv)"

  for chr in $(seq 1 22); do
    first=$(bcftools query -l ${IMPDIR}/chr${chr}.dose.vcf.gz | head -1 || true)
    if echo "$first" | grep -qP '^\d+_'; then
      echo -n "chr${chr}: stripping prefix... "
      bcftools reheader -s ${POSTDIR}/strip_prefix.tsv \
        -o ${POSTDIR}/chr${chr}.dose.vcf.gz \
        ${IMPDIR}/chr${chr}.dose.vcf.gz
      echo "done"
    else
      echo "chr${chr}: no prefix — linking"
      ln -sf ${IMPDIR}/chr${chr}.dose.vcf.gz ${POSTDIR}/chr${chr}.dose.vcf.gz
    fi
  done
  VCFDIR=${POSTDIR}
else
  echo "No numeric prefix — using imputation VCFs directly"
  VCFDIR=${IMPDIR}
fi

echo ""
echo "========================================="
echo "STEP 6: Post-Imputation QC"
echo "========================================="

echo "=== 6a: Filter variants (R2 >= 0.80, MAF >= 0.001) ==="
# R2 and MAF are in the dose VCF INFO field — filter directly with expression.
# bcftools view -i streams through the VCF once (no index needed, no targets file).
# This is ~10x faster than the old -T targets approach.

mkdir -p ${POSTDIR}/filtered

for chr in $(seq 1 22); do
  echo -n "chr${chr}: "

  bcftools view -i 'INFO/R2>=0.80 && INFO/MAF>=0.001' \
    -Oz -o ${POSTDIR}/filtered/chr${chr}.clean.vcf.gz \
    ${VCFDIR}/chr${chr}.dose.vcf.gz

  bcftools index -f ${POSTDIR}/filtered/chr${chr}.clean.vcf.gz
  NOUT=$(bcftools index -n ${POSTDIR}/filtered/chr${chr}.clean.vcf.gz)
  echo "${NOUT} pass"
done

echo ""
echo "=== 6b: Count total filtered variants ==="
TOTAL=0
for chr in $(seq 1 22); do
  n=$(bcftools index -n ${POSTDIR}/filtered/chr${chr}.clean.vcf.gz)
  TOTAL=$((TOTAL + n))
done
echo "Total variants (R2>=0.80, MAF>=0.001): $TOTAL"

echo ""
echo "=== 6c: Convert to PLINK binary ==="
# Concatenate all chromosomes into a single VCF, then convert
bcftools concat ${POSTDIR}/filtered/chr{1..22}.clean.vcf.gz -Oz \
  -o ${POSTDIR}/UZB_imputed_HQ.vcf.gz
bcftools index ${POSTDIR}/UZB_imputed_HQ.vcf.gz

plink2 --vcf ${POSTDIR}/UZB_imputed_HQ.vcf.gz dosage=DS \
  --make-bed \
  --out ${POSTDIR}/UZB_imputed_HQ

echo ""
echo "=== 6d: Sample QC — heterozygosity check ==="
plink2 --bfile ${POSTDIR}/UZB_imputed_HQ \
  --het \
  --out ${POSTDIR}/het_check

# Flag samples with |F| > 0.2 (inbreeding coefficient)
awk 'NR>1 && ($6 > 0.2 || $6 < -0.2) {print $1, $2}' \
  ${POSTDIR}/het_check.het > ${POSTDIR}/het_outliers.txt
NHET=$(wc -l < ${POSTDIR}/het_outliers.txt)
echo "Heterozygosity outliers (|F|>0.2): $NHET"

if [ "$NHET" -gt 0 ]; then
  echo "=== Removing het outliers ==="
  plink2 --bfile ${POSTDIR}/UZB_imputed_HQ \
    --remove ${POSTDIR}/het_outliers.txt \
    --make-bed \
    --out ${POSTDIR}/UZB_imputed_HQ_clean
else
  echo "No het outliers — copying files"
  cp ${POSTDIR}/UZB_imputed_HQ.bed ${POSTDIR}/UZB_imputed_HQ_clean.bed
  cp ${POSTDIR}/UZB_imputed_HQ.bim ${POSTDIR}/UZB_imputed_HQ_clean.bim
  cp ${POSTDIR}/UZB_imputed_HQ.fam ${POSTDIR}/UZB_imputed_HQ_clean.fam
fi

echo ""
echo "=== FINAL SUMMARY ==="
echo "Samples: $(wc -l < ${POSTDIR}/UZB_imputed_HQ_clean.fam)"
echo "Variants: $(wc -l < ${POSTDIR}/UZB_imputed_HQ_clean.bim)"
echo "Genotyping rate: $(head -20 ${POSTDIR}/UZB_imputed_HQ_clean.log | grep 'genotyping rate' || echo 'N/A')"

echo ""
echo "=== DISK USAGE ==="
du -sh ${POSTDIR}/
df -h /staging/

echo "=== DONE ==="
