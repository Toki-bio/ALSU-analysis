#!/bin/bash
# Step 6 only — Post-Imputation QC (spring 2026)
# Prereq: Step 5 completed (prefix-stripped VCFs in post_imputation/)
# Optimization: uses bcftools view -i expression (streaming, no index needed)
set -euo pipefail

WORKDIR=/staging/ALSU-analysis/spring2026
POSTDIR=${WORKDIR}/post_imputation
VCFDIR=${POSTDIR}  # step 5 placed reheadered/linked VCFs here

cd "$WORKDIR"

echo "========================================="
echo "STEP 6: Post-Imputation QC"
echo "========================================="
date

# Clean up any partial output from prior run
rm -f ${POSTDIR}/filtered/chr*.clean.vcf.gz ${POSTDIR}/filtered/chr*.clean.vcf.gz.csi
rm -f ${POSTDIR}/filtered/pass_chr*.txt
mkdir -p ${POSTDIR}/filtered

echo "=== 6a: Filter variants (R2 >= 0.80, MAF >= 0.001) ==="
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
  echo "No het outliers — renaming files"
  for ext in bed bim fam log; do
    cp ${POSTDIR}/UZB_imputed_HQ.${ext} ${POSTDIR}/UZB_imputed_HQ_clean.${ext}
  done
fi

echo ""
echo "=== FINAL SUMMARY ==="
echo "Samples: $(wc -l < ${POSTDIR}/UZB_imputed_HQ_clean.fam)"
echo "Variants: $(wc -l < ${POSTDIR}/UZB_imputed_HQ_clean.bim)"
echo "Genotyping rate: $(head -20 ${POSTDIR}/UZB_imputed_HQ_clean.log | grep -i 'genotyping rate' || echo 'N/A')"
echo ""
echo "=== DISK USAGE ==="
du -sh ${POSTDIR}/
df -h /staging/
date
echo "=== DONE ==="
