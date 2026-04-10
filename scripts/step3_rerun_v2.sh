#!/bin/bash
# Step 3 full re-run with --geno 0.02
# All output to /tmp/step3/ to avoid permission issues
set -e

INDIR="/staging/ALSU-analysis/spring2026"
OUTDIR="/tmp/step3_geno02"
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=== STEP 3a: PLINK SNP QC with --geno 0.02 ==="
plink --bfile ${INDIR}/ConvSK_mind20_dedup \
  --maf 0.01 \
  --hwe 1e-6 \
  --geno 0.02 \
  --chr 1-22 \
  --make-bed \
  --out ${OUTDIR}/ConvSK_mind20_dedup_snpqc

echo "SNP QC result:"
wc -l ${OUTDIR}/ConvSK_mind20_dedup_snpqc.bim | awk '{print "  Variants:", $1}'
wc -l ${OUTDIR}/ConvSK_mind20_dedup_snpqc.fam | awk '{print "  Samples:", $1}'

echo ""
echo "=== STEP 3b: Histogram data extraction ==="
plink --bfile ${INDIR}/ConvSK_mind20_dedup --chr 1-22 \
  --freq --hardy --missing \
  --out ${OUTDIR}/pre_metrics

echo "--- MAF histogram (autosomal variants on 1,093 samples) ---"
awk 'NR>1 {
  v=$5+0;
  if      (v<0.01)  b[0]++
  else if (v<0.05)  b[1]++
  else if (v<0.10)  b[2]++
  else if (v<0.15)  b[3]++
  else if (v<0.20)  b[4]++
  else if (v<0.25)  b[5]++
  else if (v<0.30)  b[6]++
  else if (v<0.35)  b[7]++
  else if (v<0.40)  b[8]++
  else if (v<0.45)  b[9]++
  else              b[10]++
}
END {
  split("0.00-0.01,0.01-0.05,0.05-0.10,0.10-0.15,0.15-0.20,0.20-0.25,0.25-0.30,0.30-0.35,0.35-0.40,0.40-0.45,0.45-0.50", labels, ",")
  for (i=0; i<=10; i++) printf "MAF_BIN\t%s\t%d\n", labels[i+1], b[i]+0
}' ${OUTDIR}/pre_metrics.frq

echo "--- HWE -log10(p) histogram ---"
awk 'NR>1 && $9 != "NA" {
  p=$9+0;
  if (p<=0) lp=20; else lp=-log(p)/log(10);
  if      (lp<1)  b[0]++
  else if (lp<2)  b[1]++
  else if (lp<3)  b[2]++
  else if (lp<4)  b[3]++
  else if (lp<5)  b[4]++
  else if (lp<6)  b[5]++
  else if (lp<10) b[6]++
  else            b[7]++
}
END {
  split("0-1,1-2,2-3,3-4,4-5,5-6,6-10,10-20", labels, ",")
  for (i=0; i<=7; i++) printf "HWE_BIN\t%s\t%d\n", labels[i+1], b[i]+0
}' ${OUTDIR}/pre_metrics.hwe

echo "--- Per-SNP missingness histogram ---"
awk 'NR>1 {
  v=$5+0;
  if      (v<0.005) b[0]++
  else if (v<0.010) b[1]++
  else if (v<0.020) b[2]++
  else if (v<0.030) b[3]++
  else if (v<0.050) b[4]++
  else if (v<0.100) b[5]++
  else if (v<0.150) b[6]++
  else if (v<0.200) b[7]++
  else if (v<0.300) b[8]++
  else              b[9]++
}
END {
  split("0.000-0.005,0.005-0.010,0.010-0.020,0.020-0.030,0.030-0.050,0.050-0.100,0.100-0.150,0.150-0.200,0.200-0.300,0.300-1.000", labels, ",")
  for (i=0; i<=9; i++) printf "GENO_BIN\t%s\t%d\n", labels[i+1], b[i]+0
}' ${OUTDIR}/pre_metrics.lmiss

echo ""
echo "=== STEP 3c: I/D allele check ==="
ID_COUNT=$(awk '($5=="I"||$5=="D"||$6=="I"||$6=="D")' ${OUTDIR}/ConvSK_mind20_dedup_snpqc.bim | wc -l)
echo "I/D allele variants: $ID_COUNT"

echo ""
echo "=== STEP 3d: Build VCF-safe SNP list and export per-chromosome VCFs ==="
awk '($5!="I" && $5!="D" && $6!="I" && $6!="D") {print $2}' \
  ${OUTDIR}/ConvSK_mind20_dedup_snpqc.bim > ${OUTDIR}/vcf_ok_snps.txt
echo "VCF-safe SNPs: $(wc -l < ${OUTDIR}/vcf_ok_snps.txt)"

for chr in $(seq 1 22); do
  plink --bfile ${OUTDIR}/ConvSK_mind20_dedup_snpqc \
    --chr $chr --extract ${OUTDIR}/vcf_ok_snps.txt \
    --recode vcf bgz --out ${OUTDIR}/chr${chr} 2>/dev/null
  tabix -f -p vcf ${OUTDIR}/chr${chr}.vcf.gz
done
echo "Per-chromosome VCFs exported (chr1-chr22)"

# Per-chr variant counts
echo "--- Per-chromosome variant counts ---"
for chr in $(seq 1 22); do
  N=$(bcftools index -n ${OUTDIR}/chr${chr}.vcf.gz)
  echo "CHR\t$chr\t$N"
done

echo ""
echo "=== STEP 3e: Concatenate + clean VCF ==="
bcftools concat -Oz -o ${OUTDIR}/impute_in.autosomes.vcf.gz \
  ${OUTDIR}/chr{1..22}.vcf.gz
tabix -f -p vcf ${OUTDIR}/impute_in.autosomes.vcf.gz

CONCAT_N=$(bcftools index -n ${OUTDIR}/impute_in.autosomes.vcf.gz)
echo "After concat: $CONCAT_N variants"

bcftools view -m2 -M2 -v snps -Oz \
  -o ${OUTDIR}/impute_in.autosomes.clean.vcf.gz \
  ${OUTDIR}/impute_in.autosomes.vcf.gz
tabix -f -p vcf ${OUTDIR}/impute_in.autosomes.clean.vcf.gz

BIALLELIC_N=$(bcftools index -n ${OUTDIR}/impute_in.autosomes.clean.vcf.gz)
echo "After biallelic filter: $BIALLELIC_N variants"

bcftools norm -d all -Oz \
  -o ${OUTDIR}/impute_in.autosomes.clean.nodup.vcf.gz \
  ${OUTDIR}/impute_in.autosomes.clean.vcf.gz
tabix -f -p vcf ${OUTDIR}/impute_in.autosomes.clean.nodup.vcf.gz

NODUP_N=$(bcftools index -n ${OUTDIR}/impute_in.autosomes.clean.nodup.vcf.gz)
DUP_REMOVED=$((BIALLELIC_N - NODUP_N))
echo "After dedup: $NODUP_N variants ($DUP_REMOVED duplicates removed)"

echo ""
echo "=== STEP 3f: Reheader (fix Cyrillic IDs) ==="
NON_ASCII=$(bcftools query -l ${OUTDIR}/impute_in.autosomes.clean.nodup.vcf.gz | grep -cP '[^\x00-\x7F]' || true)
echo "Non-ASCII sample IDs found: $NON_ASCII"

if [ "$NON_ASCII" -gt 0 ]; then
  bcftools query -l ${OUTDIR}/impute_in.autosomes.clean.nodup.vcf.gz | grep -P '[^\x00-\x7F]'
  printf '808_03-25\xd0\xbc\t808_03-25m\n' > ${OUTDIR}/rename_samples.tsv
  printf '1038_08-176\xd0\xa5-00006\t1038_08-176X-00006\n' >> ${OUTDIR}/rename_samples.tsv
  bcftools reheader -s ${OUTDIR}/rename_samples.tsv \
    -o ${OUTDIR}/impute_in.final.vcf.gz \
    ${OUTDIR}/impute_in.autosomes.clean.nodup.vcf.gz
else
  cp ${OUTDIR}/impute_in.autosomes.clean.nodup.vcf.gz ${OUTDIR}/impute_in.final.vcf.gz
fi

bcftools index ${OUTDIR}/impute_in.final.vcf.gz
FINAL_N=$(bcftools index -n ${OUTDIR}/impute_in.final.vcf.gz)
FINAL_S=$(bcftools query -l ${OUTDIR}/impute_in.final.vcf.gz | wc -l)
echo "Final VCF: $FINAL_N variants, $FINAL_S samples"

echo ""
echo "=== STEP 3g: Verification ==="
plink --bfile ${OUTDIR}/ConvSK_mind20_dedup_snpqc \
  --freq --hardy --missing \
  --out ${OUTDIR}/verify 2>/dev/null

MAF_FAIL=$(awk 'NR>1 && $5+0 < 0.01' ${OUTDIR}/verify.frq | wc -l)
HWE_FAIL=$(awk 'NR>1 && $9+0 < 1e-6 && $9 != "NA"' ${OUTDIR}/verify.hwe | wc -l)
GENO_FAIL=$(awk 'NR>1 && $5+0 > 0.02' ${OUTDIR}/verify.lmiss | wc -l)
echo "Verification: MAF<0.01 violations=$MAF_FAIL, HWE<1e-6 violations=$HWE_FAIL, geno>0.02 violations=$GENO_FAIL"

FINAL_NONASCII=$(bcftools query -l ${OUTDIR}/impute_in.final.vcf.gz | grep -cP '[^\x00-\x7F]' || true)
echo "Non-ASCII IDs in final VCF: $FINAL_NONASCII"

echo ""
echo "=== SUMMARY ==="
echo "Input: ConvSK_mind20_dedup (654,027 variants x 1,093 samples)"
echo "Autosomal: 620,492 variants"
SNPQC_N=$(wc -l < ${OUTDIR}/ConvSK_mind20_dedup_snpqc.bim)
echo "After PLINK QC (--maf 0.01 --hwe 1e-6 --geno 0.02 --chr 1-22): $SNPQC_N variants"
echo "I/D alleles removed: $ID_COUNT"
echo "Duplicate positions removed: $DUP_REMOVED"
echo "Non-ASCII IDs fixed: $NON_ASCII"
echo "Final VCF for imputation: $FINAL_N variants x $FINAL_S samples"
echo "Output location: $OUTDIR"
echo "=== DONE ==="
