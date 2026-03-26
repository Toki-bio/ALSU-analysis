#!/bin/bash
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH

# ================================================================
# V2 FULL RE-ANALYSIS — All downstream steps
# ================================================================
# Phase A: Quick analyses (PCA local, PCA global, ROH, IBD)  ~1-2 hours
# Phase B: UZB ADMIXTURE K=2..8 with 5-fold CV               ~12-24 hours
# Phase C: Global ADMIXTURE K=2..8 with 5-fold CV             ~24-48 hours
# ================================================================

V2=$HOME/v2/plink
OUT=$HOME/v2
S9=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/1000g_grch38
T=8; M=16000
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "========== V2 FULL RE-ANALYSIS START =========="
log "V2 QC: $(wc -l < $V2/UZB_v2_qc.fam) samples, $(wc -l < $V2/UZB_v2_qc.bim) variants"
log "V2 LD-pruned: $(wc -l < $V2/UZB_v2_ldpruned.fam) samples, $(wc -l < $V2/UZB_v2_ldpruned.bim) variants"
log "Disk free: $(df -h /home | awk 'NR==2{print $4}')"

# ================================================================
# PHASE A1: Uzbek-Only PCA
# ================================================================
log "===== PHASE A1: Uzbek-Only PCA ====="
mkdir -p $OUT/pca
t0=$(date +%s)
plink2 --bfile $V2/UZB_v2_ldpruned \
  --pca 10 \
  --out $OUT/pca/UZB_v2_pca \
  --threads $T --memory $M
log "UZB PCA done in $(($(date +%s)-t0))s"
log "Eigenvec: $(wc -l < $OUT/pca/UZB_v2_pca.eigenvec) lines"

# ================================================================
# PHASE A2: Global PCA (UZB V2 + 1000G)
# ================================================================
log "===== PHASE A2: Global PCA ====="
mkdir -p $OUT/global_pca

# Extract common variants between V2 and KG_reference_final (by variant ID)
t0=$(date +%s)
awk '{print $2}' $V2/UZB_v2_ldpruned.bim | sort > $OUT/global_pca/v2_snps.txt
awk '{print $2}' $S9/KG_reference_final.bim | sort > $OUT/global_pca/kg_snps.txt
comm -12 $OUT/global_pca/v2_snps.txt $OUT/global_pca/kg_snps.txt > $OUT/global_pca/common_snps.txt
n_common=$(wc -l < $OUT/global_pca/common_snps.txt)
log "Common variants: $n_common (V2: $(wc -l < $OUT/global_pca/v2_snps.txt), KG: $(wc -l < $OUT/global_pca/kg_snps.txt))"

# Extract common variants from both datasets
plink2 --bfile $V2/UZB_v2_ldpruned \
  --extract $OUT/global_pca/common_snps.txt \
  --make-bed --out $OUT/global_pca/v2_common \
  --threads $T --memory $M

plink2 --bfile $S9/KG_reference_final \
  --extract $OUT/global_pca/common_snps.txt \
  --make-bed --out $OUT/global_pca/kg_common \
  --threads $T --memory $M

# Merge using plink2 --pmerge
plink2 --bfile $OUT/global_pca/v2_common \
  --pmerge $OUT/global_pca/kg_common \
  --make-bed --out $OUT/global_pca/UZB_1kG_v2_merged \
  --threads $T --memory $M || {
  log "plink2 merge failed, trying plink1.9 --bmerge with strand handling"
  # plink1.9 handles strand flips better
  plink --bfile $OUT/global_pca/v2_common \
    --bmerge $OUT/global_pca/kg_common.bed $OUT/global_pca/kg_common.bim $OUT/global_pca/kg_common.fam \
    --make-bed --out $OUT/global_pca/UZB_1kG_v2_merged || {
    # Exclude problematic SNPs (strand/triallelic)
    log "Excluding $(wc -l < $OUT/global_pca/UZB_1kG_v2_merged-merge.missnp) mismatched SNPs"
    plink --bfile $OUT/global_pca/v2_common \
      --exclude $OUT/global_pca/UZB_1kG_v2_merged-merge.missnp \
      --make-bed --out $OUT/global_pca/v2_common_clean
    plink --bfile $OUT/global_pca/kg_common \
      --exclude $OUT/global_pca/UZB_1kG_v2_merged-merge.missnp \
      --make-bed --out $OUT/global_pca/kg_common_clean
    plink --bfile $OUT/global_pca/v2_common_clean \
      --bmerge $OUT/global_pca/kg_common_clean.bed $OUT/global_pca/kg_common_clean.bim $OUT/global_pca/kg_common_clean.fam \
      --make-bed --out $OUT/global_pca/UZB_1kG_v2_merged
  }
}
log "Merged: $(wc -l < $OUT/global_pca/UZB_1kG_v2_merged.fam) samples, $(wc -l < $OUT/global_pca/UZB_1kG_v2_merged.bim) variants"

# PCA on merged set
plink2 --bfile $OUT/global_pca/UZB_1kG_v2_merged \
  --pca 10 \
  --out $OUT/pca/GLOBAL_PCA_v2 \
  --threads $T --memory $M
log "Global PCA done in $(($(date +%s)-t0))s"
log "Global eigenvec: $(wc -l < $OUT/pca/GLOBAL_PCA_v2.eigenvec) lines"

# Copy pop_mapping for reference
cp $S9/pop_mapping.txt $OUT/pca/pop_mapping.txt 2>/dev/null || true

# ================================================================
# PHASE A3: ROH Analysis
# ================================================================
log "===== PHASE A3: ROH Analysis ====="
mkdir -p $OUT/roh
t0=$(date +%s)
plink2 --bfile $V2/UZB_v2_qc \
  --make-bed --out $OUT/roh/UZB_v2_for_roh \
  --threads $T --memory $M
# Use plink1.9 for --homozyg (plink2 doesn't support it)
plink --bfile $OUT/roh/UZB_v2_for_roh \
  --homozyg \
  --homozyg-snp 50 \
  --homozyg-kb 1000 \
  --homozyg-density 50 \
  --homozyg-gap 1000 \
  --homozyg-window-snp 50 \
  --homozyg-window-het 1 \
  --homozyg-window-missing 5 \
  --homozyg-window-threshold 0.05 \
  --out $OUT/roh/UZB_v2_ROH
log "ROH done in $(($(date +%s)-t0))s"
log "ROH segments: $(wc -l < $OUT/roh/UZB_v2_ROH.hom) (excl header)"
log "ROH per-individual: $(wc -l < $OUT/roh/UZB_v2_ROH.hom.indiv) lines"
# Clean intermediate
rm -f $OUT/roh/UZB_v2_for_roh.*

# ================================================================
# PHASE A4: IBD Analysis
# ================================================================
log "===== PHASE A4: IBD Analysis ====="
mkdir -p $OUT/ibd
t0=$(date +%s)
plink --bfile $V2/UZB_v2_ldpruned \
  --genome \
  --min 0.05 \
  --out $OUT/ibd/UZB_v2_IBD
log "IBD done in $(($(date +%s)-t0))s"
log "IBD pairs (PI_HAT>=0.05): $(tail -n +2 $OUT/ibd/UZB_v2_IBD.genome | wc -l)"

log "===== PHASE A COMPLETE ====="
log "Disk free: $(df -h /home | awk 'NR==2{print $4}')"

# ================================================================
# PHASE B: UZB-Only ADMIXTURE (K=2..8)
# ================================================================
log "===== PHASE B: UZB ADMIXTURE K=2..8 ====="
mkdir -p $OUT/admixture
cd $OUT/admixture

# Copy input files to working dir (ADMIXTURE writes output in cwd)
cp $V2/UZB_v2_ldpruned.bed UZB_v2_admix.bed
cp $V2/UZB_v2_ldpruned.bim UZB_v2_admix.bim
cp $V2/UZB_v2_ldpruned.fam UZB_v2_admix.fam

for K in $(seq 2 8); do
  t0=$(date +%s)
  log "Starting ADMIXTURE K=$K (5-fold CV)..."
  admixture --cv=5 -j$T UZB_v2_admix.bed $K 2>&1 | tee UZB_v2_admix_K${K}.log
  elapsed=$(($(date +%s)-t0))
  log "ADMIXTURE K=$K done in ${elapsed}s"
  # Rename output files  
  mv UZB_v2_admix.${K}.Q UZB_v2_admix.${K}.Q 2>/dev/null || true
  mv UZB_v2_admix.${K}.P UZB_v2_admix.${K}.P 2>/dev/null || true
done
log "===== PHASE B COMPLETE: UZB ADMIXTURE ====="

# Extract CV errors
log "CV Errors:"
for K in $(seq 2 8); do
  cv=$(grep -i 'CV error' UZB_v2_admix_K${K}.log 2>/dev/null | tail -1)
  log "  K=$K: $cv"
done

# ================================================================
# PHASE C: Global ADMIXTURE (UZB V2 + 1000G)
# ================================================================
log "===== PHASE C: Global ADMIXTURE ====="
mkdir -p $OUT/global_admixture

# Use the already-merged global PCA dataset for ADMIXTURE
# It's already LD-pruned (intersection of V2_ldpruned and KG_reference_final)
cp $OUT/global_pca/UZB_1kG_v2_merged.bed $OUT/global_admixture/global_v2_admix.bed
cp $OUT/global_pca/UZB_1kG_v2_merged.bim $OUT/global_admixture/global_v2_admix.bim
cp $OUT/global_pca/UZB_1kG_v2_merged.fam $OUT/global_admixture/global_v2_admix.fam

cd $OUT/global_admixture
for K in $(seq 2 8); do
  t0=$(date +%s)
  log "Starting Global ADMIXTURE K=$K (5-fold CV)..."
  admixture --cv=5 -j$T global_v2_admix.bed $K 2>&1 | tee global_v2_admix_K${K}.log
  elapsed=$(($(date +%s)-t0))
  log "Global ADMIXTURE K=$K done in ${elapsed}s"
done
log "===== PHASE C COMPLETE: Global ADMIXTURE ====="

# Global CV errors
log "Global CV Errors:"
for K in $(seq 2 8); do
  cv=$(grep -i 'CV error' global_v2_admix_K${K}.log 2>/dev/null | tail -1)
  log "  K=$K: $cv"
done

# ================================================================
# SUMMARY
# ================================================================
log "========== FULL RE-ANALYSIS COMPLETE =========="
log "Results:"
log "  UZB PCA:       $OUT/pca/UZB_v2_pca.{eigenvec,eigenval}"
log "  Global PCA:    $OUT/pca/GLOBAL_PCA_v2.{eigenvec,eigenval}"
log "  ROH:           $OUT/roh/UZB_v2_ROH.{hom,hom.indiv,hom.summary}"
log "  IBD:           $OUT/ibd/UZB_v2_IBD.genome"
log "  UZB ADMIX:     $OUT/admixture/UZB_v2_admix.{2-8}.Q"
log "  Global ADMIX:  $OUT/global_admixture/global_v2_admix.{2-8}.Q"
log "Disk free: $(df -h /home | awk 'NR==2{print $4}')"
