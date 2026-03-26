#!/bin/bash
###############################################################################
# ALSU v2 Validation Script
#
# Run AFTER recalc_v2.sh completes.
# Performs detailed comparison between v1 (buggy) and v2 (corrected) datasets.
#
# Usage: bash validate_v2.sh > validate_v2.log 2>&1
###############################################################################
set -euo pipefail

export PATH=/staging/conda/envs/bioinfo/bin:$PATH

# Paths
V1=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz/filtered_clean
V2=/staging/ALSU-analysis/winter2025/v2_position_filtered/plink

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "===== V1 vs V2 Validation ====="
echo ""

###############################################################################
# 1. Variant count comparison
###############################################################################
log "--- 1. Variant Counts ---"
echo ""
printf "%-35s %12s %12s\n" "Dataset" "V1 (buggy)" "V2 (fixed)"
printf "%-35s %12s %12s\n" "---" "---" "---"

v1_clean=$(wc -l < "$V1/UZB_imputed_HQ_clean.bim")
v2_clean=$(wc -l < "$V2/UZB_v2_clean.bim")
printf "%-35s %12d %12d\n" "Pre-QC (after VCF filter)" "$v1_clean" "$v2_clean"

v1_qc=$(wc -l < "$V1/UZB_imputed_HQ_qc.bim")
v2_qc=$(wc -l < "$V2/UZB_v2_qc.bim")
printf "%-35s %12d %12d\n" "Post-QC" "$v1_qc" "$v2_qc"

v1_samples=$(wc -l < "$V1/UZB_imputed_HQ_qc.fam")
v2_samples=$(wc -l < "$V2/UZB_v2_qc.fam")
printf "%-35s %12d %12d\n" "Post-QC samples" "$v1_samples" "$v2_samples"

echo ""
extra_v1=$((v1_clean - v2_clean))
log "Extra variants in V1 (the bug): $extra_v1 (expected ~837,039)"

###############################################################################
# 2. Dot-ID check
###############################################################################
log ""
log "--- 2. Dot-ID Contamination Check ---"

v1_dots=$(awk '$2=="."' "$V1/UZB_imputed_HQ_clean.bim" | wc -l)
v2_dots=$(awk '$2=="."' "$V2/UZB_v2_clean.bim" | wc -l)
log "Dot-IDs in V1 pre-QC BIM: $v1_dots (expected ~813,191)"
log "Dot-IDs in V2 pre-QC BIM: $v2_dots (expected 0 or near-0)"

if [ "$v2_dots" -eq 0 ]; then
  log "PASS: No dot-ID contamination in V2"
else
  log "WARNING: $v2_dots dot-IDs still present — investigate"
fi

###############################################################################
# 3. Sample rescue check
###############################################################################
log ""
log "--- 3. Sample Rescue Analysis ---"

if [ "$v2_samples" -gt "$v1_samples" ]; then
  rescued=$((v2_samples - v1_samples))
  log "RESCUED $rescued samples that were removed in V1 by --mind 0.05"
  log "Explanation: V1's 837K spurious variants had high missingness in some"
  log "  samples. With those variants removed, those samples pass --mind 0.05."

  # Find which samples were rescued
  v2_ids=$(awk '{print $2}' "$V2/UZB_v2_qc.fam" | sort)
  v1_ids=$(awk '{print $2}' "$V1/UZB_imputed_HQ_qc.fam" | sort)
  log "Rescued samples:"
  comm -23 <(echo "$v2_ids") <(echo "$v1_ids") | while read -r s; do
    echo "  + $s"
  done
elif [ "$v2_samples" -eq "$v1_samples" ]; then
  log "Same sample count in V1 and V2 (no samples rescued)"
else
  diff=$((v1_samples - v2_samples))
  log "WARNING: V2 has $diff FEWER samples than V1 — investigate"
fi

###############################################################################
# 4. Quality metrics of spurious variants
###############################################################################
log ""
log "--- 4. Quality Profile of Spurious Variants ---"

# Extract V1-only variants (spurious)
log "Extracting V1-only variants..."
awk '{print $1":"$4}' "$V1/UZB_imputed_HQ_clean.bim" | sort > /tmp/v1_chrpos.txt
awk '{print $1":"$4}' "$V2/UZB_v2_clean.bim" | sort > /tmp/v2_chrpos.txt
comm -23 /tmp/v1_chrpos.txt /tmp/v2_chrpos.txt > /tmp/v1_only_chrpos.txt
spurious=$(wc -l < /tmp/v1_only_chrpos.txt)
log "V1-only positions: $spurious"

# How many of the spurious survived QC?
if [ -f "$V1/UZB_imputed_HQ_qc.bim" ]; then
  awk '{print $1":"$4}' "$V1/UZB_imputed_HQ_qc.bim" | sort > /tmp/v1_qc_chrpos.txt
  survived=$(comm -12 /tmp/v1_only_chrpos.txt /tmp/v1_qc_chrpos.txt | wc -l)
  log "Spurious variants that survived V1 QC: $survived"
  if [ "$spurious" -gt 0 ]; then
    pct=$(echo "scale=1; $survived * 100 / $spurious" | bc)
    log "  ($pct% of spurious variants passed --geno/--maf/--mind filters)"
  fi
fi

###############################################################################
# 5. LD pruning comparison
###############################################################################
log ""
log "--- 5. LD Pruning Summary ---"

if [ -f "$V2/UZB_v2_ldpruned.bim" ]; then
  v2_ldp=$(wc -l < "$V2/UZB_v2_ldpruned.bim")
  log "V2 LD-pruned variants: $v2_ldp"
fi

# V1 LD-pruned (if accessible)
if [ -f "$V1/UZB_imputed_HQ_unique.bim" ]; then
  v1_unique=$(wc -l < "$V1/UZB_imputed_HQ_unique.bim")
  log "V1 unique-ID variants: $v1_unique"
fi

###############################################################################
# 6. Quick PCA sanity check
###############################################################################
log ""
log "--- 6. PCA (quick 10-component) ---"

V2_PCA=$V2/../pca
mkdir -p "$V2_PCA"

if [ -f "$V2_PCA/UZB_v2_pca.eigenval" ]; then
  log "PCA already computed, skipping"
else
  log "Computing 10 PCs for quick comparison..."
  plink2 --bfile "$V2/UZB_v2_ldpruned" \
    --pca 10 \
    --out "$V2_PCA/UZB_v2_pca" \
    --threads 8
fi

if [ -f "$V2_PCA/UZB_v2_pca.eigenval" ]; then
  log "PCA eigenvalues (V2):"
  cat "$V2_PCA/UZB_v2_pca.eigenval"
fi

###############################################################################
log ""
log "===== VALIDATION COMPLETE ====="
