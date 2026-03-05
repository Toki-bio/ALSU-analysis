#!/bin/bash
set -e

###############################################################################
# Global ADMIXTURE Pipeline: Uzbek Cohort + 1000 Genomes Reference Populations
# 
# Populations included:
#   UZBEK: 1,199 samples (our cohort, lifted over to hg19)
#   EUR:   503 samples (CEU=99, GBR=91, FIN=99, IBS=107, TSI=107)
#   SAS:   199 samples (GIH=103, PJL=96)
#   EAS:   207 samples (CHB=103, JPT=104)
#   AFR:   108 samples (YRI=108)
#
# Total: ~2,216 samples
# Uses 1000G Phase 3 VCFs (GRCh37/hg19) + Uzbek data lifted to hg19
###############################################################################

WORKDIR=/staging/ALSU-analysis/admixture_analysis/global_admixture
FSTDIR=/staging/ALSU-analysis/Fst_analysis
VCF_DIR=${FSTDIR}/1000G_data
PANEL=${VCF_DIR}/integrated_call_samples_v3.20130502.ALL.panel
UZB_DATA=${FSTDIR}/uzbek_data_hg19
THREADS=32

mkdir -p ${WORKDIR}/temp ${WORKDIR}/admix_results
cd ${WORKDIR}

LOG=${WORKDIR}/pipeline.log
exec > >(tee -a ${LOG}) 2>&1

echo "============================================================"
echo " Global ADMIXTURE Pipeline"
echo " Start: $(date)"
echo "============================================================"

#==============================================================================
# STEP 1: Create population sample lists from 1000G panel
#==============================================================================
echo ""
echo "[Step 1] Creating population sample lists..."

# EUR: CEU, GBR, FIN, IBS, TSI (all 5 European populations)
awk 'NR>1 && ($3=="EUR"){print $1,$1}' ${PANEL} > temp/EUR_samples.txt

# SAS: GIH (Gujarati) + PJL (Punjabi) — closest South Asians to Central Asia
awk 'NR>1 && ($2=="GIH" || $2=="PJL"){print $1,$1}' ${PANEL} > temp/SAS_samples.txt

# EAS: CHB (Han Chinese Beijing) + JPT (Japanese Tokyo)
awk 'NR>1 && ($2=="CHB" || $2=="JPT"){print $1,$1}' ${PANEL} > temp/EAS_samples.txt

# AFR: YRI (Yoruba Ibadan) — outgroup baseline
awk 'NR>1 && ($2=="YRI"){print $1,$1}' ${PANEL} > temp/AFR_samples.txt

# Combine all reference samples
cat temp/EUR_samples.txt temp/SAS_samples.txt temp/EAS_samples.txt temp/AFR_samples.txt > temp/all_ref_samples.txt

echo "  EUR: $(wc -l < temp/EUR_samples.txt) samples"
echo "  SAS: $(wc -l < temp/SAS_samples.txt) samples"
echo "  EAS: $(wc -l < temp/EAS_samples.txt) samples"
echo "  AFR: $(wc -l < temp/AFR_samples.txt) samples"
echo "  Total reference: $(wc -l < temp/all_ref_samples.txt) samples"
echo "  Uzbek: $(wc -l < ${UZB_DATA}.fam) samples"

#==============================================================================
# STEP 2: Process each chromosome — extract, match, merge
#==============================================================================
echo ""
echo "[Step 2] Processing chromosomes..."

for chr in $(seq 1 22); do
    echo ""
    echo "  --- Chromosome ${chr} ---"
    
    VCF=${VCF_DIR}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
    
    if [ ! -f "${VCF}" ]; then
        echo "  WARNING: VCF for chr${chr} not found, skipping!"
        continue
    fi
    
    # 2a. Extract reference samples from VCF → plink
    echo "  [a] Extracting reference samples from VCF..."
    plink --vcf ${VCF} \
          --keep temp/all_ref_samples.txt \
          --snps-only just-acgt \
          --biallelic-only strict \
          --make-bed \
          --out temp/ref_chr${chr} \
          --silent
    
    # 2b. Rename SNPs to chr:pos format for matching
    echo "  [b] Renaming SNPs to chr:pos format..."
    awk '{$2=$1":"$4; print}' OFS='\t' temp/ref_chr${chr}.bim > temp/ref_chr${chr}_renamed.bim
    cp temp/ref_chr${chr}.bed temp/ref_chr${chr}_renamed.bed
    cp temp/ref_chr${chr}.fam temp/ref_chr${chr}_renamed.fam
    
    # Remove duplicate positions (keep first)
    awk '!seen[$2]++' temp/ref_chr${chr}_renamed.bim > temp/ref_chr${chr}_dedup_snps.txt
    cut -f2 temp/ref_chr${chr}_dedup_snps.txt > temp/ref_chr${chr}_keep.txt
    plink --bfile temp/ref_chr${chr}_renamed \
          --extract temp/ref_chr${chr}_keep.txt \
          --make-bed \
          --out temp/ref_chr${chr}_dedup \
          --silent
    
    # 2c. Extract Uzbek data for this chromosome
    echo "  [c] Extracting Uzbek chr${chr}..."
    plink --bfile ${UZB_DATA} \
          --chr ${chr} \
          --make-bed \
          --out temp/uzb_chr${chr} \
          --silent
    
    # Rename Uzbek SNPs to chr:pos
    awk '{$2=$1":"$4; print}' OFS='\t' temp/uzb_chr${chr}.bim > temp/uzb_chr${chr}_renamed.bim
    cp temp/uzb_chr${chr}.bed temp/uzb_chr${chr}_renamed.bed
    cp temp/uzb_chr${chr}.fam temp/uzb_chr${chr}_renamed.fam
    
    # Remove Uzbek duplicate positions
    awk '!seen[$2]++' temp/uzb_chr${chr}_renamed.bim > temp/uzb_chr${chr}_dedup_snps.txt
    cut -f2 temp/uzb_chr${chr}_dedup_snps.txt > temp/uzb_chr${chr}_keep.txt
    plink --bfile temp/uzb_chr${chr}_renamed \
          --extract temp/uzb_chr${chr}_keep.txt \
          --make-bed \
          --out temp/uzb_chr${chr}_dedup \
          --silent
    
    # 2d. Find overlapping positions
    echo "  [d] Finding overlapping SNPs..."
    awk '{print $2}' temp/ref_chr${chr}_dedup.bim | sort > temp/ref_chr${chr}_snplist.txt
    awk '{print $2}' temp/uzb_chr${chr}_dedup.bim | sort > temp/uzb_chr${chr}_snplist.txt
    comm -12 temp/ref_chr${chr}_snplist.txt temp/uzb_chr${chr}_snplist.txt > temp/overlap_chr${chr}.txt
    
    OVERLAP=$(wc -l < temp/overlap_chr${chr}.txt)
    echo "    Ref SNPs: $(wc -l < temp/ref_chr${chr}_dedup.bim)"
    echo "    Uzb SNPs: $(wc -l < temp/uzb_chr${chr}_dedup.bim)"
    echo "    Overlap:  ${OVERLAP}"
    
    if [ "${OVERLAP}" -eq 0 ]; then
        echo "    WARNING: No overlap for chr${chr}, skipping!"
        continue
    fi
    
    # 2e. Extract overlapping SNPs from both datasets
    plink --bfile temp/uzb_chr${chr}_dedup \
          --extract temp/overlap_chr${chr}.txt \
          --make-bed \
          --out temp/uzb_chr${chr}_matched \
          --silent
    
    plink --bfile temp/ref_chr${chr}_dedup \
          --extract temp/overlap_chr${chr}.txt \
          --make-bed \
          --out temp/ref_chr${chr}_matched \
          --silent
    
    # 2f. Merge (with conflict resolution)
    echo "  [e] Merging..."
    plink --bfile temp/uzb_chr${chr}_matched \
          --bmerge temp/ref_chr${chr}_matched \
          --make-bed \
          --out temp/merged_chr${chr} \
          --silent 2>&1 || true
    
    # Check for merge conflicts
    if [ -f temp/merged_chr${chr}-merge.missnp ]; then
        CONFLICTS=$(wc -l < temp/merged_chr${chr}-merge.missnp)
        echo "    Resolving ${CONFLICTS} allele conflicts..."
        
        plink --bfile temp/uzb_chr${chr}_matched \
              --exclude temp/merged_chr${chr}-merge.missnp \
              --make-bed \
              --out temp/uzb_chr${chr}_clean \
              --silent
        
        plink --bfile temp/ref_chr${chr}_matched \
              --exclude temp/merged_chr${chr}-merge.missnp \
              --make-bed \
              --out temp/ref_chr${chr}_clean \
              --silent
        
        plink --bfile temp/uzb_chr${chr}_clean \
              --bmerge temp/ref_chr${chr}_clean \
              --make-bed \
              --out temp/merged_chr${chr} \
              --silent
    fi
    
    MERGED_SNPS=$(wc -l < temp/merged_chr${chr}.bim)
    MERGED_SAMPLES=$(wc -l < temp/merged_chr${chr}.fam)
    echo "    ✓ Chr${chr}: ${MERGED_SNPS} SNPs × ${MERGED_SAMPLES} samples"
done

#==============================================================================
# STEP 3: Merge all chromosomes into one dataset
#==============================================================================
echo ""
echo "[Step 3] Merging all chromosomes..."

# Create merge list (chr2-22)
> temp/merge_list.txt
for chr in $(seq 2 22); do
    if [ -f temp/merged_chr${chr}.bed ]; then
        echo "temp/merged_chr${chr}" >> temp/merge_list.txt
    fi
done

plink --bfile temp/merged_chr1 \
      --merge-list temp/merge_list.txt \
      --make-bed \
      --out global_merged \
      --silent

TOTAL_SAMPLES=$(wc -l < global_merged.fam)
TOTAL_SNPS=$(wc -l < global_merged.bim)
echo "  Merged dataset: ${TOTAL_SAMPLES} samples × ${TOTAL_SNPS} SNPs"

#==============================================================================
# STEP 4: QC — remove high-missingness SNPs and samples
#==============================================================================
echo ""
echo "[Step 4] QC filtering..."

plink --bfile global_merged \
      --geno 0.05 \
      --mind 0.1 \
      --maf 0.01 \
      --make-bed \
      --out global_qc \
      --silent

QC_SAMPLES=$(wc -l < global_qc.fam)
QC_SNPS=$(wc -l < global_qc.bim)
echo "  After QC: ${QC_SAMPLES} samples × ${QC_SNPS} SNPs"

#==============================================================================
# STEP 5: LD pruning
#==============================================================================
echo ""
echo "[Step 5] LD pruning..."

plink --bfile global_qc \
      --indep-pairwise 50 5 0.2 \
      --out global_ld \
      --silent

plink --bfile global_qc \
      --extract global_ld.prune.in \
      --make-bed \
      --out global_for_admixture \
      --silent

LD_SNPS=$(wc -l < global_for_admixture.bim)
echo "  After LD pruning: ${LD_SNPS} SNPs retained"

#==============================================================================
# STEP 6: Create population labels file
#==============================================================================
echo ""
echo "[Step 6] Creating population labels..."

# Build lookup: sample → population (from 1000G panel)
awk 'NR>1{print $1,$2,$3}' ${PANEL} > temp/panel_lookup.txt

# For each sample in .fam, assign population
awk '
BEGIN {
    # Load panel lookup
    while ((getline line < "temp/panel_lookup.txt") > 0) {
        split(line, a, " ")
        pop[a[1]] = a[2]
        superpop[a[1]] = a[3]
    }
}
{
    iid = $2
    if (iid in pop) {
        printf "%s %s %s %s\n", $1, $2, pop[iid], superpop[iid]
    } else {
        printf "%s %s UZBEK UZB\n", $1, $2
    }
}
' global_for_admixture.fam > global_pop_labels.txt

echo "  Population counts:"
awk '{print $3}' global_pop_labels.txt | sort | uniq -c | sort -rn
echo ""
echo "  Super-population counts:"
awk '{print $4}' global_pop_labels.txt | sort | uniq -c | sort -rn

# Create .pop file for supervised ADMIXTURE (optional)
# "-" means unknown (let ADMIXTURE infer), population name means fixed
# For unsupervised: we don't need this, but create for reference
awk '{print $4}' global_pop_labels.txt > global_superpop.pop

#==============================================================================
# STEP 7: Run ADMIXTURE K=2 through K=8 (unsupervised)
#==============================================================================
echo ""
echo "[Step 7] Running ADMIXTURE K=2-8 (unsupervised)..."

cd admix_results

for K in $(seq 2 8); do
    echo ""
    echo "  >>> K=${K} starting at $(date)..."
    
    admixture --cv ${WORKDIR}/global_for_admixture.bed ${K} -j${THREADS} \
        > admixture_K${K}.log 2>&1
    
    CV=$(grep 'CV error' admixture_K${K}.log | awk '{print $NF}')
    LL=$(grep 'Loglikelihood' admixture_K${K}.log | tail -1 | awk '{print $2}')
    echo "  <<< K=${K} done. CV=${CV}, LogL=${LL}"
    
    # Copy Q/P files with clear names
    cp global_for_admixture.${K}.Q K${K}.Q
    cp global_for_admixture.${K}.P K${K}.P
done

cd ${WORKDIR}

#==============================================================================
# STEP 8: Summary
#==============================================================================
echo ""
echo "============================================================"
echo " Pipeline Complete!"
echo " End: $(date)"
echo "============================================================"
echo ""
echo "Results in: ${WORKDIR}/admix_results/"
echo ""
echo "CV errors:"
for K in $(seq 2 8); do
    grep 'CV error' admix_results/admixture_K${K}.log 2>/dev/null || echo "  K=${K}: not yet complete"
done
echo ""
echo "Files:"
echo "  global_for_admixture.bed/bim/fam — LD-pruned merged dataset"
echo "  global_pop_labels.txt — FID IID POP SUPERPOP"
echo "  admix_results/K{2-8}.Q — Q matrices"
echo "  admix_results/K{2-8}.P — P matrices (allele frequencies)"
echo ""
echo "Population composition:"
awk '{print $3}' global_pop_labels.txt | sort | uniq -c | sort -rn
