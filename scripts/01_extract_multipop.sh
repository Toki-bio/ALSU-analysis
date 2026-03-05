#!/bin/bash
# ============================================================================
# Step 1: Extract EAS, SAS, AFR sample lists and generate per-population
#         PLINK files from 1000 Genomes Phase 3 VCFs
#
# Prerequisites:
#   - 1000G VCFs in 1000G_data/
#   - Panel file: integrated_call_samples_v3.20130502.ALL.panel
#   - Uzbek hg19 data: uzbek_data_hg19.{bed,bim,fam}
#
# Output:
#   - Per-population merged PLINK files (all chromosomes)
#   - Allele frequency files for each population
# ============================================================================
set -e

BASEDIR="/staging/ALSU-analysis/Fst_analysis"
OUTDIR="/staging/ALSU-analysis/admixture_analysis/pop_specific"
PANEL="${BASEDIR}/1000G_data/integrated_call_samples_v3.20130502.ALL.panel"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "===== Multi-Population Extraction Pipeline ====="
echo "Start time: $(date)"
echo ""

# ---- Step 1: Create sample lists for each superpopulation ----
echo "[Step 1] Creating population sample lists..."

# EAS samples (CHB, JPT, CHS, CDX, KHV)
awk -F'\t' '$3=="EAS" {print $1, $1}' "${PANEL}" > eas_samples.txt
EAS_N=$(wc -l < eas_samples.txt)
echo "  EAS samples: ${EAS_N}"

# SAS samples (GIH, PJL, BEB, STU, ITU)
awk -F'\t' '$3=="SAS" {print $1, $1}' "${PANEL}" > sas_samples.txt
SAS_N=$(wc -l < sas_samples.txt)
echo "  SAS samples: ${SAS_N}"

# AFR samples (YRI, LWK, GWD, MSL, ESN, ACB, ASW)
awk -F'\t' '$3=="AFR" {print $1, $1}' "${PANEL}" > afr_samples.txt
AFR_N=$(wc -l < afr_samples.txt)
echo "  AFR samples: ${AFR_N}"

echo ""

# ---- Step 2: Get the list of positions that exist in UZB-EUR overlap ----
# Use the merged_all_chrs.bim which has 376K overlapping SNPs
echo "[Step 2] Extracting target SNP positions from UZB-EUR merged data..."
awk '{print $2}' "${BASEDIR}/merged_all_chrs.bim" > target_snps.txt
TARGET_N=$(wc -l < target_snps.txt)
echo "  Target SNPs: ${TARGET_N}"

# Also create per-chromosome position lists for bcftools extraction
for chr in {1..22}; do
    awk -v c=${chr} '$1==c {print $1":"$4}' "${BASEDIR}/merged_all_chrs.bim" > target_pos_chr${chr}.txt
done
echo ""

# ---- Step 3: Process each population across all chromosomes ----
for POP in EAS SAS AFR; do
    pop_lower=$(echo ${POP} | tr '[:upper:]' '[:lower:]')
    SAMPLE_FILE="${pop_lower}_samples.txt"
    
    echo "===== Processing ${POP} ====="
    
    for chr in {1..22}; do
        echo "  Chr${chr}: Converting 1000G VCF..."
        
        # Convert VCF to PLINK, keeping only this population's samples
        plink --vcf "${BASEDIR}/1000G_data/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
              --keep "${SAMPLE_FILE}" \
              --make-bed \
              --out "1000G_${POP}_chr${chr}_raw" \
              --silent
        
        # Rename SNP IDs to chr:pos format (matching UZB-EUR convention)
        awk '{$2=$1":"$4; print}' OFS='\t' "1000G_${POP}_chr${chr}_raw.bim" > "1000G_${POP}_chr${chr}.bim"
        cp "1000G_${POP}_chr${chr}_raw.bed" "1000G_${POP}_chr${chr}.bed"
        cp "1000G_${POP}_chr${chr}_raw.fam" "1000G_${POP}_chr${chr}.fam"
        
        # Extract only the target SNPs (matching UZB-EUR overlap set)
        plink --bfile "1000G_${POP}_chr${chr}" \
              --extract target_snps.txt \
              --make-bed \
              --out "1000G_${POP}_chr${chr}_matched" \
              --silent
        
        snps=$(wc -l < "1000G_${POP}_chr${chr}_matched.bim")
        echo "    ${POP} chr${chr}: ${snps} overlapping SNPs"
        
        # Cleanup raw files
        rm -f "1000G_${POP}_chr${chr}_raw".*
        rm -f "1000G_${POP}_chr${chr}.bed" "1000G_${POP}_chr${chr}.bim" "1000G_${POP}_chr${chr}.fam"
    done
    
    # ---- Merge all chromosomes for this population ----
    echo "  Merging all chromosomes for ${POP}..."
    
    # Create merge list (chr2-22)
    > "merge_list_${POP}.txt"
    for chr in {2..22}; do
        echo "1000G_${POP}_chr${chr}_matched" >> "merge_list_${POP}.txt"
    done
    
    plink --bfile "1000G_${POP}_chr1_matched" \
          --merge-list "merge_list_${POP}.txt" \
          --make-bed \
          --out "1000G_${POP}_all" \
          --silent
    
    POP_SNPS=$(wc -l < "1000G_${POP}_all.bim")
    POP_SAMPLES=$(wc -l < "1000G_${POP}_all.fam")
    echo "  ✓ ${POP} complete: ${POP_SAMPLES} samples × ${POP_SNPS} SNPs"
    
    # ---- Compute allele frequencies ----
    echo "  Computing ${POP} allele frequencies..."
    plink --bfile "1000G_${POP}_all" \
          --freq \
          --out "freq_${POP}" \
          --silent
    
    echo "  ✓ Frequencies written to freq_${POP}.frq"
    
    # Cleanup per-chromosome files
    for chr in {1..22}; do
        rm -f "1000G_${POP}_chr${chr}_matched".*
        rm -f "target_pos_chr${chr}.txt"
    done
    
    echo ""
done

# ---- Step 4: Compute UZB and EUR frequencies from existing data ----
echo "===== Computing UZB and EUR frequencies ====="

# UZB frequencies: extract UZB samples from merged_all_chrs
echo "  Computing UZB allele frequencies..."
plink --bfile "${BASEDIR}/merged_all_chrs" \
      --keep <(awk '$3=="UZBEK" {print $1, $2}' "${BASEDIR}/populations.txt") \
      --freq \
      --out "freq_UZB" \
      --silent
echo "  ✓ UZB frequencies written"

echo "  Computing EUR allele frequencies..."
plink --bfile "${BASEDIR}/merged_all_chrs" \
      --keep <(awk '$3=="EUR" {print $1, $2}' "${BASEDIR}/populations.txt") \
      --freq \
      --out "freq_EUR" \
      --silent
echo "  ✓ EUR frequencies written"

echo ""

# ---- Step 5: Compute pairwise Fst ----
echo "===== Computing Pairwise Fst ====="

# We need UZB+EAS and UZB+SAS merged files for PLINK --fst
# Approach: merge UZB (from merged_all_chrs) with each new population

# Extract UZB-only from merged_all_chrs
echo "  Extracting UZB-only PLINK files..."
plink --bfile "${BASEDIR}/merged_all_chrs" \
      --keep <(awk '$3=="UZBEK" {print $1, $2}' "${BASEDIR}/populations.txt") \
      --make-bed \
      --out "UZB_only" \
      --silent

# Extract EUR-only from merged_all_chrs
echo "  Extracting EUR-only PLINK files..."
plink --bfile "${BASEDIR}/merged_all_chrs" \
      --keep <(awk '$3=="EUR" {print $1, $2}' "${BASEDIR}/populations.txt") \
      --make-bed \
      --out "EUR_only" \
      --silent

for POP in EAS SAS AFR; do
    echo ""
    echo "  --- UZB vs ${POP} ---"
    
    # Find common SNPs
    comm -12 <(awk '{print $2}' UZB_only.bim | sort) \
             <(awk '{print $2}' "1000G_${POP}_all.bim" | sort) > "common_snps_UZB_${POP}.txt"
    COMMON_N=$(wc -l < "common_snps_UZB_${POP}.txt")
    echo "  Common SNPs: ${COMMON_N}"
    
    # Extract common SNPs from both
    plink --bfile UZB_only \
          --extract "common_snps_UZB_${POP}.txt" \
          --make-bed --out "UZB_for_${POP}" --silent
    
    plink --bfile "1000G_${POP}_all" \
          --extract "common_snps_UZB_${POP}.txt" \
          --make-bed --out "${POP}_for_merge" --silent
    
    # Merge
    plink --bfile "UZB_for_${POP}" \
          --bmerge "${POP}_for_merge" \
          --make-bed \
          --out "merged_UZB_${POP}" \
          --silent 2>&1 || {
        echo "  Merge conflicts detected, excluding problem SNPs..."
        plink --bfile "UZB_for_${POP}" \
              --exclude "merged_UZB_${POP}-merge.missnp" \
              --make-bed --out "UZB_for_${POP}_clean" --silent
        plink --bfile "${POP}_for_merge" \
              --exclude "merged_UZB_${POP}-merge.missnp" \
              --make-bed --out "${POP}_for_merge_clean" --silent
        plink --bfile "UZB_for_${POP}_clean" \
              --bmerge "${POP}_for_merge_clean" \
              --make-bed --out "merged_UZB_${POP}" --silent
        rm -f "UZB_for_${POP}_clean".* "${POP}_for_merge_clean".*
    }
    
    # Create population file
    awk '{print $1, $2, "UZBEK"}' "UZB_for_${POP}.fam" > "pop_UZB_${POP}.txt"
    awk -v p="${POP}" '{print $1, $2, p}' "${POP}_for_merge.fam" >> "pop_UZB_${POP}.txt"
    
    # Compute Fst
    plink --bfile "merged_UZB_${POP}" \
          --fst \
          --within "pop_UZB_${POP}.txt" \
          --out "fst_UZB_vs_${POP}" \
          --silent
    
    echo "  ✓ UZB vs ${POP} Fst computed"
    tail -3 "fst_UZB_vs_${POP}.log" | grep -i "fst"
    
    # Cleanup
    rm -f "UZB_for_${POP}".* "${POP}_for_merge".* "common_snps_UZB_${POP}.txt"
done

# EUR vs EAS Fst (needed for PBS triangle)
echo ""
echo "  --- EUR vs EAS (for PBS) ---"
comm -12 <(awk '{print $2}' EUR_only.bim | sort) \
         <(awk '{print $2}' 1000G_EAS_all.bim | sort) > common_snps_EUR_EAS.txt

plink --bfile EUR_only --extract common_snps_EUR_EAS.txt --make-bed --out EUR_for_EAS --silent
plink --bfile 1000G_EAS_all --extract common_snps_EUR_EAS.txt --make-bed --out EAS_for_EUR --silent

plink --bfile EUR_for_EAS --bmerge EAS_for_EUR --make-bed --out merged_EUR_EAS --silent 2>&1 || {
    plink --bfile EUR_for_EAS --exclude merged_EUR_EAS-merge.missnp --make-bed --out EUR_for_EAS_clean --silent
    plink --bfile EAS_for_EUR --exclude merged_EUR_EAS-merge.missnp --make-bed --out EAS_for_EUR_clean --silent
    plink --bfile EUR_for_EAS_clean --bmerge EAS_for_EUR_clean --make-bed --out merged_EUR_EAS --silent
    rm -f EUR_for_EAS_clean.* EAS_for_EUR_clean.*
}

awk '{print $1, $2, "EUR"}' EUR_for_EAS.fam > pop_EUR_EAS.txt
awk '{print $1, $2, "EAS"}' EAS_for_EUR.fam >> pop_EUR_EAS.txt

plink --bfile merged_EUR_EAS --fst --within pop_EUR_EAS.txt --out fst_EUR_vs_EAS --silent
echo "  ✓ EUR vs EAS Fst computed"
tail -3 fst_EUR_vs_EAS.log | grep -i "fst"

rm -f EUR_for_EAS.* EAS_for_EUR.* common_snps_EUR_EAS.txt

echo ""
echo "===== All pairwise Fst computed! ====="
echo "End time: $(date)"
echo ""
echo "Output files:"
echo "  freq_UZB.frq, freq_EUR.frq, freq_EAS.frq, freq_SAS.frq, freq_AFR.frq"
echo "  fst_UZB_vs_EAS.fst, fst_UZB_vs_SAS.fst, fst_UZB_vs_AFR.fst"
echo "  fst_EUR_vs_EAS.fst (for PBS)"
echo "  genomewide_uzbek_vs_eur_fst.fst (existing, in Fst_analysis/)"
echo ""
echo "Next: run calculate_pbs.py to compute PBS and identify Uzbek-specific SNPs"
