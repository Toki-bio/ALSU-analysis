#!/bin/bash
# Steps 7 + 8 + 11: PCA + 1000G merge + ADMIXTURE (spring 2026)
# Prereq: Step 6 completed — UZB_imputed_HQ_clean.{bed,bim,fam} exists
set -euo pipefail

WORKDIR=/staging/ALSU-analysis/spring2026/post_imputation
REFDIR=/staging/ALSU-analysis/Fst_analysis
KGDIR=${REFDIR}/1000G_data
cd "$WORKDIR"

echo "========================================="
echo "STEP 7: Local PCA (Uzbek-only)"
echo "========================================="
date

# 7a: Basic QC on imputed dataset
echo "=== 7a: QC filters (geno 0.05, mind 0.05, maf 0.01) ==="
plink2 --bfile UZB_imputed_HQ_clean \
  --geno 0.05 --mind 0.05 --maf 0.01 \
  --make-bed --out UZB_imputed_HQ_qc \
  --threads 8 --memory 200000

echo "QC result:"
wc -l < UZB_imputed_HQ_qc.fam
wc -l < UZB_imputed_HQ_qc.bim

# 7b: Assign unique variant IDs
echo "=== 7b: Unique variant IDs ==="
plink2 --bfile UZB_imputed_HQ_qc \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 50 \
  --make-bed --out UZB_imputed_HQ_unique \
  --threads 8

# 7c: LD pruning (aggressive for PCA)
echo "=== 7c: LD pruning (1000kb, step 1, r²<0.05) ==="
plink2 --bfile UZB_imputed_HQ_unique \
  --indep-pairwise 1000kb 1 0.05 \
  --out UZB_pruned --threads 8

echo "Pruned variants:"
wc -l < UZB_pruned.prune.in

# 7d: PCA
echo "=== 7d: PCA (10 PCs) ==="
plink2 --bfile UZB_imputed_HQ_unique \
  --extract UZB_pruned.prune.in \
  --pca 10 --out UZB_final_pca --threads 8

echo "PCA eigenvalues:"
cat UZB_final_pca.eigenval

echo ""
echo "========================================="
echo "STEP 8: Global PCA with 1000 Genomes"
echo "========================================="
date

# 8a: Extract matching SNPs from 1000G
echo "=== 8a: Extract pruned SNPs from 1000G ==="
mkdir -p ${WORKDIR}/1kg_merge

# Convert pruned variant list to BED format for bcftools
awk -F':' '{
  chr=$1; pos=$2
  # Handle both "chr1" and "1" formats
  if(chr !~ /^chr/) chr="chr"chr
  print chr"\t"(pos-1)"\t"pos
}' UZB_pruned.prune.in | sort -k1,1 -k2,2n > ${WORKDIR}/1kg_merge/pruned_targets.bed

for chr in $(seq 1 22); do
  KG_VCF=$(ls ${KGDIR}/ALL.chr${chr}.*.vcf.gz 2>/dev/null | head -1)
  if [ -z "$KG_VCF" ]; then
    echo "WARN: No 1000G VCF for chr${chr}"
    continue
  fi

  # Extract only our targets, biallelic SNPs
  bcftools view -T ${WORKDIR}/1kg_merge/pruned_targets.bed \
    -m2 -M2 -v snps \
    "$KG_VCF" \
    -Oz -o ${WORKDIR}/1kg_merge/chr${chr}_subset.vcf.gz 2>/dev/null || true

  if [ -s "${WORKDIR}/1kg_merge/chr${chr}_subset.vcf.gz" ]; then
    plink2 --vcf ${WORKDIR}/1kg_merge/chr${chr}_subset.vcf.gz \
      --set-all-var-ids '@:#:$r:$a' \
      --new-id-max-allele-len 50 \
      --make-bed --out ${WORKDIR}/1kg_merge/chr${chr}_ref \
      --threads 8 2>/dev/null
    echo "chr${chr}: $(wc -l < ${WORKDIR}/1kg_merge/chr${chr}_ref.bim) variants"
  fi
done

# 8b: Merge 1000G chromosomes
echo "=== 8b: Merge 1000G chromosomes ==="
ls ${WORKDIR}/1kg_merge/chr*_ref.bed | sed 's/.bed$//' | tail -n +2 > ${WORKDIR}/1kg_merge/merge_list.txt
FIRST=$(ls ${WORKDIR}/1kg_merge/chr*_ref.bed | sed 's/.bed$//' | head -1)

plink --bfile "$FIRST" \
  --merge-list ${WORKDIR}/1kg_merge/merge_list.txt \
  --make-bed --out ${WORKDIR}/1kg_merge/KG_reference_final \
  --allow-no-sex

echo "1000G reference:"
wc -l < ${WORKDIR}/1kg_merge/KG_reference_final.fam
wc -l < ${WORKDIR}/1kg_merge/KG_reference_final.bim

# 8c: Merge Uzbek + 1000G
echo "=== 8c: Merge Uzbek + 1000G ==="
plink --bfile UZB_imputed_HQ_unique \
  --bmerge ${WORKDIR}/1kg_merge/KG_reference_final \
  --make-bed --out ${WORKDIR}/UZB_1kG_merged \
  --allow-no-sex

echo "Merged dataset:"
wc -l < ${WORKDIR}/UZB_1kG_merged.fam
wc -l < ${WORKDIR}/UZB_1kG_merged.bim

# 8d: Global PCA
echo "=== 8d: Global PCA (10 PCs) ==="
plink2 --bfile ${WORKDIR}/UZB_1kG_merged \
  --pca 10 --out ${WORKDIR}/GLOBAL_PCA --threads 8

echo "Global PCA eigenvalues:"
cat ${WORKDIR}/GLOBAL_PCA.eigenval

echo ""
echo "========================================="
echo "STEP 11: ADMIXTURE Analysis"
echo "========================================="
date

# 11a: QC on merged dataset
echo "=== 11a: QC for ADMIXTURE ==="
plink --bfile ${WORKDIR}/UZB_1kG_merged \
  --maf 0.01 --geno 0.02 --hwe 1e-6 \
  --make-bed --out ${WORKDIR}/global_qc \
  --allow-no-sex

# 11b: LD pruning (ADMIXTURE-specific, less aggressive)
echo "=== 11b: LD pruning for ADMIXTURE (50, 10, 0.1) ==="
plink --bfile ${WORKDIR}/global_qc \
  --indep-pairwise 50 10 0.1 \
  --out ${WORKDIR}/global_prune \
  --allow-no-sex
plink --bfile ${WORKDIR}/global_qc \
  --extract ${WORKDIR}/global_prune.prune.in \
  --make-bed --out ${WORKDIR}/global_for_admixture \
  --allow-no-sex

echo "ADMIXTURE input:"
wc -l < ${WORKDIR}/global_for_admixture.fam
wc -l < ${WORKDIR}/global_for_admixture.bim

# 11c: Run ADMIXTURE K=2-8
echo "=== 11c: ADMIXTURE K=2 to K=8 ==="
cd ${WORKDIR}
for K in $(seq 2 8); do
  echo "--- K=$K ---"
  admixture --cv -j32 global_for_admixture.bed $K 2>&1 | tee ${WORKDIR}/admixture_log_K${K}.out
done

# 11d: Extract CV errors
echo ""
echo "=== CV Error Summary ==="
for K in $(seq 2 8); do
  grep -i "cv error" ${WORKDIR}/admixture_log_K${K}.out || echo "K=$K: no CV line"
done

echo ""
echo "=== ADMIXTURE DONE ==="
date
echo "========================================="
