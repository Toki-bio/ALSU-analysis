#!/bin/bash
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH
V2=$HOME/v2; P=$V2/plink; PC=$V2/perchr; ML=$PC/merge_list.txt
V1=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz/filtered_clean
REF=/staging/Genomes/Human/chr/GRCh38.fa; T=8; M=16000
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

rm -f $P/UZB_v2_clean.{bed,bim,fam,pgen,pvar,psam,log}
log "===== PHASE 3 (retry): Merge + assign IDs ====="
first=$(head -1 "$ML"); tail -n +2 "$ML" > "${ML}.rest"
t0=$(date +%s)
plink2 --pfile "$first" --pmerge-list "${ML}.rest" \
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 \
  --fa "$REF" --ref-from-fa force --make-bed \
  --out "$P/UZB_v2_clean" --threads $T --memory $M
log "Merge done in $(($(date +%s)-t0))s"
log "PLINK: $(wc -l < $P/UZB_v2_clean.fam) samples, $(wc -l < $P/UZB_v2_clean.bim) variants"
log "Dot-IDs: $(awk '$2=="."' $P/UZB_v2_clean.bim | wc -l)"

log "Cleaning per-chr files..."; rm -rf "$PC"; log "Freed space"

log "===== PHASE 4: QC ====="
t0=$(date +%s)
plink2 --bfile "$P/UZB_v2_clean" --geno 0.05 --mind 0.05 --maf 0.01 \
  --make-bed --out "$P/UZB_v2_qc" --threads $T --memory $M
log "QC done in $(($(date +%s)-t0))s"
log "V2 QC: $(wc -l < $P/UZB_v2_qc.fam) samples, $(wc -l < $P/UZB_v2_qc.bim) variants"
log "V1 QC: $(wc -l < $V1/UZB_imputed_HQ_qc.fam 2>/dev/null || echo ?) samples, $(wc -l < $V1/UZB_imputed_HQ_qc.bim 2>/dev/null || echo ?) variants"

log "===== PHASE 5: LD pruning ====="
t0=$(date +%s)
plink2 --bfile "$P/UZB_v2_qc" --indep-pairwise 1000kb 1 0.05 --out "$P/UZB_v2_pruned" --threads $T
log "LD pruning done in $(($(date +%s)-t0))s: $(wc -l < $P/UZB_v2_pruned.prune.in) independent variants"
plink2 --bfile "$P/UZB_v2_qc" --extract "$P/UZB_v2_pruned.prune.in" --make-bed --out "$P/UZB_v2_ldpruned" --threads $T
log "LD-pruned: $(wc -l < $P/UZB_v2_ldpruned.fam) samples, $(wc -l < $P/UZB_v2_ldpruned.bim) variants"

log "===== ALL DONE: $(date) ====="