#!/bin/bash
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH
V2_BASE=$HOME/v2
V2_PERCHR=$V2_BASE/perchr
V2_PLINK=$V2_BASE/plink
V1_DIR=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz/filtered_clean
REF_GENOME=/staging/Genomes/Human/chr/GRCh38.fa
THREADS=8
MEM_MB=16000
MERGE_LIST=$V2_PERCHR/merge_list.txt
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# Clean up failed merge attempt
rm -f $V2_PLINK/UZB_v2_clean.{bed,bim,fam,pgen,pvar,psam,log}

log "===== PHASE 3 (retry): Merge + assign IDs ====="
first=$(head -1 "$MERGE_LIST")
tail -n +2 "$MERGE_LIST" > "${MERGE_LIST}.rest"
t0=$(date +%s)
plink2 --pfile "$first" \
  --pmerge-list "${MERGE_LIST}.rest" \
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 50 \
  --fa "$REF_GENOME" --ref-from-fa force \
  --make-bed \
  --out "$V2_PLINK/UZB_v2_clean" \
  --threads $THREADS --memory $MEM_MB
t1=$(date +%s)
log "Merge done in $((t1-t0))s"
v2_variants=$(wc -l < "$V2_PLINK/UZB_v2_clean.bim")
v2_samples=$(wc -l < "$V2_PLINK/UZB_v2_clean.fam")
log "PLINK binary: $v2_samples samples, $v2_variants variants"
v2_dots=$(awk '$2=="."' "$V2_PLINK/UZB_v2_clean.bim" | wc -l)
log "Dot-ID variants: $v2_dots"

# Clean per-chr files to save space
log "Cleaning per-chromosome pgen files..."
rm -rf "$V2_PERCHR"
log "Freed per-chr space"

log "===== PHASE 4: Standard QC ====="
t0=$(date +%s)
plink2 --bfile "$V2_PLINK/UZB_v2_clean" \
  --geno 0.05 --mind 0.05 --maf 0.01 \
  --make-bed \
  --out "$V2_PLINK/UZB_v2_qc" \
  --threads $THREADS --memory $MEM_MB
t1=$(date +%s)
log "QC done in $((t1-t0))s"
v2_qc_s=$(wc -l < "$V2_PLINK/UZB_v2_qc.fam")
v2_qc_v=$(wc -l < "$V2_PLINK/UZB_v2_qc.bim")
log "After QC: $v2_qc_s samples, $v2_qc_v variants"
v1_s=$(wc -l < "$V1_DIR/UZB_imputed_HQ_qc.fam" 2>/dev/null || echo "?")
v1_v=$(wc -l < "$V1_DIR/UZB_imputed_HQ_qc.bim" 2>/dev/null || echo "?")
log "V1 (buggy):     $v1_s samples, $v1_v variants"
log "V2 (corrected): $v2_qc_s samples, $v2_qc_v variants"

log "===== PHASE 5: LD pruning ====="
t0=$(date +%s)
plink2 --bfile "$V2_PLINK/UZB_v2_qc" \
  --indep-pairwise 1000kb 1 0.05 \
  --out "$V2_PLINK/UZB_v2_pruned" \
  --threads $THREADS
t1=$(date +%s)
log "LD pruning done in $((t1-t0))s"
prune_in=$(wc -l < "$V2_PLINK/UZB_v2_pruned.prune.in")
log "LD-independent variants: $prune_in"

plink2 --bfile "$V2_PLINK/UZB_v2_qc" \
  --extract "$V2_PLINK/UZB_v2_pruned.prune.in" \
  --make-bed \
  --out "$V2_PLINK/UZB_v2_ldpruned" \
  --threads $THREADS
ldp_s=$(wc -l < "$V2_PLINK/UZB_v2_ldpruned.fam")
ldp_v=$(wc -l < "$V2_PLINK/UZB_v2_ldpruned.bim")
log "LD-pruned set: $ldp_s samples, $ldp_v variants"

log "===== SUMMARY ====="
log "Pre-QC:    $v2_variants variants, $v2_samples samples"
log "Post-QC:   $v2_qc_v variants, $v2_qc_s samples (V1: $v1_v / $v1_s)"
log "LD-pruned: $ldp_v variants, $ldp_s samples"
log "Dot-IDs:   $v2_dots"
log "===== ALL DONE: $(date) ====="