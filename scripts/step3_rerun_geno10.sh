#!/bin/bash
set -e
SRC=/staging/ALSU-analysis/spring2026
OUT=/tmp/step3_geno10
mkdir -p ${OUT}
cd ${OUT}

echo "=== Step 3a: PLINK QC (recovered winter params) ==="
plink --bfile ${SRC}/ConvSK_mind20_dedup \
  --maf 0.01 \
  --hwe 1e-6 \
  --geno 0.10 \
  --chr 1-22 \
  --make-bed \
  --out ${OUT}/ConvSK_mind20_dedup_snpqc 2>&1 | tail -15

echo ""
echo "=== Step 3b: VCF-safe SNP list ==="
awk '($5!="I" && $5!="D" && $6!="I" && $6!="D") {print $2}' \
  ${OUT}/ConvSK_mind20_dedup_snpqc.bim > ${OUT}/vcf_ok_snps.txt
TOTAL_BIM=$(wc -l < ${OUT}/ConvSK_mind20_dedup_snpqc.bim)
VCF_OK=$(wc -l < ${OUT}/vcf_ok_snps.txt)
ID_REMOVED=$((TOTAL_BIM - VCF_OK))
echo "Total after PLINK QC: ${TOTAL_BIM}"
echo "I/D alleles removed: ${ID_REMOVED}"
echo "VCF-safe SNPs: ${VCF_OK}"

echo ""
echo "=== Step 3c: Per-chromosome VCF export ==="
for chr in $(seq 1 22); do
  plink --bfile ${OUT}/ConvSK_mind20_dedup_snpqc \
    --chr $chr --extract ${OUT}/vcf_ok_snps.txt \
    --recode vcf bgz --out ${OUT}/chr${chr} 2>/dev/null
  tabix -f -p vcf ${OUT}/chr${chr}.vcf.gz
  N=$(bcftools view -H ${OUT}/chr${chr}.vcf.gz | wc -l)
  echo "chr${chr}: ${N} variants"
done

echo ""
echo "=== Step 3d: Concatenate and clean ==="
bcftools concat -Oz -o ${OUT}/impute_in.autosomes.vcf.gz ${OUT}/chr{1..22}.vcf.gz
tabix -f -p vcf ${OUT}/impute_in.autosomes.vcf.gz
echo "After concat: $(bcftools index -n ${OUT}/impute_in.autosomes.vcf.gz)"

bcftools view -m2 -M2 -v snps -Oz \
  -o ${OUT}/impute_in.autosomes.clean.vcf.gz ${OUT}/impute_in.autosomes.vcf.gz
tabix -f -p vcf ${OUT}/impute_in.autosomes.clean.vcf.gz
echo "After biallelic filter: $(bcftools index -n ${OUT}/impute_in.autosomes.clean.vcf.gz)"

bcftools norm -d all -Oz \
  -o ${OUT}/impute_in.autosomes.clean.nodup.vcf.gz \
  ${OUT}/impute_in.autosomes.clean.vcf.gz 2>&1
tabix -f -p vcf ${OUT}/impute_in.autosomes.clean.nodup.vcf.gz
echo "After dedup: $(bcftools index -n ${OUT}/impute_in.autosomes.clean.nodup.vcf.gz)"

echo ""
echo "=== Step 3e: Fix Cyrillic IDs ==="
printf '808_03-25\xd0\xbc\t808_03-25m\n' > ${OUT}/rename_samples.tsv
printf '1038_08-176\xd0\xa5-00006\t1038_08-176X-00006\n' >> ${OUT}/rename_samples.tsv

bcftools reheader -s ${OUT}/rename_samples.tsv \
  -o ${OUT}/impute_in.final.vcf.gz \
  ${OUT}/impute_in.autosomes.clean.nodup.vcf.gz
tabix -f -p vcf ${OUT}/impute_in.final.vcf.gz

echo "Final variants: $(bcftools index -n ${OUT}/impute_in.final.vcf.gz)"
echo "Final samples: $(bcftools query -l ${OUT}/impute_in.final.vcf.gz | wc -l)"
echo "Non-ASCII IDs: $(bcftools query -l ${OUT}/impute_in.final.vcf.gz | grep -cP '[^\x00-\x7F]' || true)"

echo ""
echo "=== Verification ==="
plink --bfile ${OUT}/ConvSK_mind20_dedup_snpqc \
  --freq --hardy --missing \
  --out ${OUT}/snpqc_verify 2>/dev/null
echo "MAF violations: $(awk 'NR>1 && $5+0 < 0.01' ${OUT}/snpqc_verify.frq | wc -l)"
echo "HWE violations: $(awk 'NR>1 && $9+0 < 1e-6 && $9 != "NA"' ${OUT}/snpqc_verify.hwe | wc -l)"
echo "Geno violations: $(awk 'NR>1 && $5+0 > 0.10' ${OUT}/snpqc_verify.lmiss | wc -l)"

echo ""
echo "=== Per-chr final VCF files for Michigan ==="
# Split final into per-chromosome for upload
for chr in $(seq 1 22); do
  bcftools view -r ${chr} -Oz -o ${OUT}/impute_chr${chr}.vcf.gz ${OUT}/impute_in.final.vcf.gz
  tabix -f -p vcf ${OUT}/impute_chr${chr}.vcf.gz
  echo "impute_chr${chr}.vcf.gz: $(bcftools index -n ${OUT}/impute_chr${chr}.vcf.gz) variants, $(stat -c%s ${OUT}/impute_chr${chr}.vcf.gz) bytes"
done

echo ""
echo "=== TOTAL SIZE ==="
du -sh ${OUT}/impute_chr*.vcf.gz | tail -1
ls -lh ${OUT}/impute_chr*.vcf.gz | awk '{sum+=$5} END {print "Total: " sum " bytes"}'
du -ch ${OUT}/impute_chr*.vcf.gz | tail -1

echo ""
echo "=== DONE ==="
