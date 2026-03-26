#!/bin/bash
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH
V2=$HOME/v2; P=$V2/plink; PC=$V2/perchr
V1=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz/filtered_clean
REF=/staging/Genomes/Human/chr/GRCh38.fa; T=8; M=16000
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "===== SPACE-EFFICIENT PHASE 3-5 (pgen→bed→merge) ====="
log "Disk before cleanup: $(df -h /home | awk 'NR==2{print $4}') free"

# Step 1: Delete failed merge output
rm -f $P/UZB_v2_clean.{pgen,pvar,psam,bed,bim,fam,log}
log "Disk after cleanup: $(df -h /home | awk 'NR==2{print $4}') free"

# Step 2: Convert per-chr pgen → bed (one by one, delete pgen after each)
# This shrinks ~15GB pgen (with dosage) to ~3GB bed (hard calls only)
log "--- Converting per-chr pgen → bed ---"
for chr in $(seq 1 22); do
    pf="$PC/chr${chr}"
    bf="$PC/chr${chr}_bed"
    if [ -f "${pf}.pgen" ]; then
        plink2 --pfile "$pf" \
            --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 \
            --make-bed --out "$bf" --threads $T --memory $M 2>&1 | tail -3
        nv=$(wc -l < "${bf}.bim")
        rm -f "${pf}.pgen" "${pf}.pvar" "${pf}.psam" "${pf}.log"
        log "chr${chr}: ${nv} variants → bed, pgen deleted"
    fi
done
log "Disk after conversion: $(df -h /home | awk 'NR==2{print $4}') free"

# Step 3: Merge all bed files
log "===== PHASE 3: Merge bed files ====="
ML="$PC/merge_bed_list.txt"
> "$ML"
first=""
for chr in $(seq 1 22); do
    bf="$PC/chr${chr}_bed"
    if [ -f "${bf}.bed" ]; then
        if [ -z "$first" ]; then
            first="$bf"
        else
            echo "$bf" >> "$ML"
        fi
    fi
done
nmerge=$(wc -l < "$ML")
log "Merging $first + $nmerge files"
t0=$(date +%s)
plink2 --bfile "$first" --pmerge-list "$ML" bfile \
    --fa "$REF" --ref-from-fa force \
    --make-bed --out "$P/UZB_v2_clean" --threads $T --memory $M
log "Merge done in $(($(date +%s)-t0))s"
log "PLINK: $(wc -l < $P/UZB_v2_clean.fam) samples, $(wc -l < $P/UZB_v2_clean.bim) variants"
log "Dot-IDs: $(awk '$2=="."' $P/UZB_v2_clean.bim | wc -l)"

# Clean per-chr bed files
rm -rf "$PC"
log "Disk after merge cleanup: $(df -h /home | awk 'NR==2{print $4}') free"

# Phase 4: QC
log "===== PHASE 4: QC ====="
t0=$(date +%s)
plink2 --bfile "$P/UZB_v2_clean" --geno 0.05 --mind 0.05 --maf 0.01 \
    --make-bed --out "$P/UZB_v2_qc" --threads $T --memory $M
log "QC done in $(($(date +%s)-t0))s"
log "V2 QC: $(wc -l < $P/UZB_v2_qc.fam) samples, $(wc -l < $P/UZB_v2_qc.bim) variants"
log "V1 QC: $(wc -l < $V1/UZB_imputed_HQ_qc.fam 2>/dev/null || echo ?) samples, $(wc -l < $V1/UZB_imputed_HQ_qc.bim 2>/dev/null || echo ?) variants"

# Phase 5: LD pruning
log "===== PHASE 5: LD pruning ====="
t0=$(date +%s)
plink2 --bfile "$P/UZB_v2_qc" --indep-pairwise 1000kb 1 0.05 --out "$P/UZB_v2_pruned" --threads $T
log "LD pruning done in $(($(date +%s)-t0))s: $(wc -l < $P/UZB_v2_pruned.prune.in) independent variants"
plink2 --bfile "$P/UZB_v2_qc" --extract "$P/UZB_v2_pruned.prune.in" --make-bed --out "$P/UZB_v2_ldpruned" --threads $T
log "LD-pruned: $(wc -l < $P/UZB_v2_ldpruned.fam) samples, $(wc -l < $P/UZB_v2_ldpruned.bim) variants"

log "===== ALL DONE: $(date) ====="
