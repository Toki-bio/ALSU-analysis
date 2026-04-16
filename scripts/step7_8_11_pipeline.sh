#!/bin/bash
# Spring 2026: Steps 7 (Local PCA) + 8 (Global PCA) + 11 (ADMIXTURE)
# Run AFTER step6 produces UZB_imputed_HQ_clean.{bed,bim,fam}
# Run on DRAGEN: /staging/ALSU-analysis/spring2026/
set -euo pipefail

WORKDIR=/staging/ALSU-analysis/spring2026
POSTDIR=${WORKDIR}/post_imputation
PCADIR=${WORKDIR}/pca
GLOBALDIR=${WORKDIR}/global_pca
ADMIXDIR=${WORKDIR}/admixture
REF1KG=/staging/ALSU-analysis/Fst_analysis/1000G_data
PANEL=${REF1KG}/integrated_call_samples_v3.20130502.ALL.panel
GENOME=/staging/Genomes/Human/hg38.fa

mkdir -p "$PCADIR" "$GLOBALDIR" "$ADMIXDIR"

echo "============================================="
echo "STEP 7: Local PCA on Uzbek Cohort"
echo "============================================="
echo "$(date)"

# 7a: Verify input
echo "=== 7a: Input validation ==="
NSAMPLES=$(wc -l < ${POSTDIR}/UZB_imputed_HQ_clean.fam)
NVARS=$(wc -l < ${POSTDIR}/UZB_imputed_HQ_clean.bim)
echo "Input: ${NSAMPLES} samples × ${NVARS} variants"

# 7b: Basic QC
echo "=== 7b: Basic QC (--geno 0.05 --mind 0.05 --maf 0.01) ==="
plink2 --bfile ${POSTDIR}/UZB_imputed_HQ_clean \
  --geno 0.05 --mind 0.05 --maf 0.01 \
  --make-bed \
  --out ${PCADIR}/UZB_qc

echo "After QC: $(wc -l < ${PCADIR}/UZB_qc.fam) samples × $(wc -l < ${PCADIR}/UZB_qc.bim) variants"

# 7c: Unique variant IDs (CHR:POS:REF:ALT)
echo "=== 7c: Assign unique variant IDs ==="
plink2 --bfile ${PCADIR}/UZB_qc \
  --set-all-var-ids '@:#:$r:$a' \
  --make-bed \
  --out ${PCADIR}/UZB_unique

# 7d: LD pruning
echo "=== 7d: LD pruning (window=1000kb step=1 r²<0.05) ==="
plink2 --bfile ${PCADIR}/UZB_unique \
  --indep-pairwise 1000kb 1 0.05 \
  --out ${PCADIR}/UZB_pruned

NPRUNED=$(wc -l < ${PCADIR}/UZB_pruned.prune.in)
echo "Independent SNPs: ${NPRUNED}"

# 7e: PCA
echo "=== 7e: PCA (10 components) ==="
plink2 --bfile ${PCADIR}/UZB_unique \
  --extract ${PCADIR}/UZB_pruned.prune.in \
  --pca 10 \
  --out ${PCADIR}/UZB_final_pca

echo "PCA eigenvalues:"
cat ${PCADIR}/UZB_final_pca.eigenval 2>/dev/null || echo "(see .eigenval file)"

echo ""
echo "============================================="
echo "STEP 8: Global PCA with 1000 Genomes"
echo "============================================="
echo "$(date)"

# 8a: Get spring UZB variant positions for 1000G extraction
echo "=== 8a: Extract UZB SNP positions ==="
awk '{print $1, $4}' ${PCADIR}/UZB_unique.bim > ${GLOBALDIR}/uzb_positions.txt
echo "UZB positions: $(wc -l < ${GLOBALDIR}/uzb_positions.txt)"

# 8b: Extract matching variants from 1000G per-chr PLINK files
echo "=== 8b: Extract from 1000G ==="
> ${GLOBALDIR}/merge_list.txt
for chr in $(seq 1 22); do
  BED=${REF1KG}/1000G_EUR_chr${chr}
  if [ -f "${BED}.bed" ]; then
    # Assign same ID format to 1000G
    plink2 --bfile "$BED" \
      --set-all-var-ids '@:#:$r:$a' \
      --make-bed \
      --out ${GLOBALDIR}/KG_chr${chr}_id

    # Extract overlapping positions
    awk -v c="chr${chr}" '$1==c {print $1":"$2}' ${GLOBALDIR}/uzb_positions.txt \
      > ${GLOBALDIR}/extract_chr${chr}.txt

    if [ -s ${GLOBALDIR}/extract_chr${chr}.txt ]; then
      plink2 --bfile ${GLOBALDIR}/KG_chr${chr}_id \
        --extract range ${GLOBALDIR}/extract_chr${chr}.txt \
        --make-bed \
        --out ${GLOBALDIR}/KG_chr${chr}_filt
      echo "${GLOBALDIR}/KG_chr${chr}_filt" >> ${GLOBALDIR}/merge_list.txt
    fi
  fi
done

# 8c: Merge 1000G chromosomes
echo "=== 8c: Merge 1000G chromosomes ==="
FIRST=$(head -1 ${GLOBALDIR}/merge_list.txt)
tail -n +2 ${GLOBALDIR}/merge_list.txt | awk '{print $0".bed",$0".bim",$0".fam"}' \
  > ${GLOBALDIR}/kg_merge_list.txt
plink --bfile "$FIRST" \
  --merge-list ${GLOBALDIR}/kg_merge_list.txt \
  --make-bed \
  --out ${GLOBALDIR}/KG_reference_final

echo "1000G merged: $(wc -l < ${GLOBALDIR}/KG_reference_final.fam) samples × $(wc -l < ${GLOBALDIR}/KG_reference_final.bim) variants"

# 8d: Find common variants and merge with UZB
echo "=== 8d: Find common SNPs ==="
# Get variant IDs in both datasets
awk '{print $2}' ${PCADIR}/UZB_unique.bim | sort > ${GLOBALDIR}/uzb_snps.txt
awk '{print $2}' ${GLOBALDIR}/KG_reference_final.bim | sort > ${GLOBALDIR}/kg_snps.txt
comm -12 ${GLOBALDIR}/uzb_snps.txt ${GLOBALDIR}/kg_snps.txt > ${GLOBALDIR}/common_snps.txt
echo "Common SNPs: $(wc -l < ${GLOBALDIR}/common_snps.txt)"

# Extract common SNPs from both
plink2 --bfile ${PCADIR}/UZB_unique \
  --extract ${GLOBALDIR}/common_snps.txt \
  --make-bed \
  --out ${GLOBALDIR}/UZB_common

plink2 --bfile ${GLOBALDIR}/KG_reference_final \
  --extract ${GLOBALDIR}/common_snps.txt \
  --make-bed \
  --out ${GLOBALDIR}/KG_common

# 8e: Merge UZB + 1000G
echo "=== 8e: Merge UZB + 1000G ==="
plink --bfile ${GLOBALDIR}/UZB_common \
  --bmerge ${GLOBALDIR}/KG_common.bed ${GLOBALDIR}/KG_common.bim ${GLOBALDIR}/KG_common.fam \
  --make-bed \
  --out ${GLOBALDIR}/UZB_1kG_merged

echo "Merged: $(wc -l < ${GLOBALDIR}/UZB_1kG_merged.fam) samples × $(wc -l < ${GLOBALDIR}/UZB_1kG_merged.bim) variants"

# 8f: Global PCA
echo "=== 8f: Global PCA (10 components) ==="
plink2 --bfile ${GLOBALDIR}/UZB_1kG_merged \
  --pca 10 \
  --out ${GLOBALDIR}/GLOBAL_PCA

echo "Global PCA eigenvalues:"
cat ${GLOBALDIR}/GLOBAL_PCA.eigenval 2>/dev/null || echo "(see .eigenval file)"

# 8g: Create population mapping
echo "=== 8g: Population mapping ==="
# UZB samples: label as UZBEK/UZB
awk '{print $1, $2, "UZBEK", "UZB"}' ${GLOBALDIR}/UZB_common.fam > ${GLOBALDIR}/pop_labels.txt
# 1000G samples: extract from panel file
if [ -f "$PANEL" ]; then
  awk -F'\t' 'NR>1 {print $1, $1, $2, $3}' "$PANEL" | \
    grep -f <(awk '{print $2}' ${GLOBALDIR}/KG_common.fam) \
    >> ${GLOBALDIR}/pop_labels.txt || true
fi
echo "Pop labels: $(wc -l < ${GLOBALDIR}/pop_labels.txt)"

echo ""
echo "============================================="
echo "STEP 11: ADMIXTURE with Cross-Validation"
echo "============================================="
echo "$(date)"

# 11a: Prepare data — LD-prune the merged dataset
echo "=== 11a: LD pruning for ADMIXTURE ==="
plink2 --bfile ${GLOBALDIR}/UZB_1kG_merged \
  --maf 0.01 --geno 0.02 --hwe 1e-6 \
  --make-bed \
  --out ${ADMIXDIR}/global_qc

plink2 --bfile ${ADMIXDIR}/global_qc \
  --indep-pairwise 50 10 0.1 \
  --out ${ADMIXDIR}/global_ld

plink2 --bfile ${ADMIXDIR}/global_qc \
  --extract ${ADMIXDIR}/global_ld.prune.in \
  --make-bed \
  --out ${ADMIXDIR}/global_for_admixture

echo "ADMIXTURE input: $(wc -l < ${ADMIXDIR}/global_for_admixture.fam) samples × $(wc -l < ${ADMIXDIR}/global_for_admixture.bim) variants"

# 11b: Run ADMIXTURE K=2 to K=8 with cross-validation
echo "=== 11b: ADMIXTURE K=2..8 ==="
cd ${ADMIXDIR}
for K in $(seq 2 8); do
  echo -n "K=${K}: "
  admixture --cv -j32 ${ADMIXDIR}/global_for_admixture.bed $K 2>&1 | tee ${ADMIXDIR}/admixture_K${K}.log | grep -i "cv error"
done

# 11c: Summary
echo ""
echo "=== CV Error Summary ==="
for K in $(seq 2 8); do
  grep -i "cv error" ${ADMIXDIR}/admixture_K${K}.log 2>/dev/null || echo "K=${K}: (not found)"
done

echo ""
echo "=== FINAL SUMMARY ==="
echo "Local PCA: $(wc -l < ${PCADIR}/UZB_final_pca.eigenvec) samples, ${NPRUNED} LD-pruned SNPs"
echo "Global PCA: $(wc -l < ${GLOBALDIR}/GLOBAL_PCA.eigenvec) samples"
echo "ADMIXTURE: $(wc -l < ${ADMIXDIR}/global_for_admixture.fam) samples × $(wc -l < ${ADMIXDIR}/global_for_admixture.bim) variants"
echo ""
echo "Output locations:"
echo "  Local PCA: ${PCADIR}/"
echo "  Global PCA: ${GLOBALDIR}/"
echo "  ADMIXTURE: ${ADMIXDIR}/"

echo ""
echo "=== DISK ==="
du -sh ${PCADIR} ${GLOBALDIR} ${ADMIXDIR}
df -h /staging/

echo "=== DONE ==="
echo "$(date)"
