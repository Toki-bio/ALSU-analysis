#!/bin/bash
###############################################################################
# ALSU Pipeline Fork v2 — Position-Based VCF Filtering
#
# Fixes the bcftools ID=@ dot-matching bug from Dec 25, 2025.
# Original: 10,009,530 HQ variants intended → 10,846,569 actual (837K extra)
# Fix: switch from ID-based to position-based variant selection
#
# Space-efficient design: converts each chromosome VCF→PLINK individually,
# then merges PLINK files. No large VCF intermediates stored. Fits in ~10GB.
#
# Run with: nohup bash /tmp/recalc_v2.sh > ~/v2/recalc_v2.log 2>&1 &
###############################################################################
set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────
THREADS=8
MEM_MB=16000

# Source paths (READ ONLY — never modified)
UNZ=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz
HQ_POS_FILE=$UNZ/UZB_all.HQ.query.chrposrefalt.tsv        # 10,009,530 lines: CHR POS REF ALT ID
HET_EXCLUDE=/staging/ALSU-analysis/winter2025/3_post-imputation/het_outliers_vcf_format.txt
REF_GENOME=/staging/Genomes/Human/chr/GRCh38.fa

# v1 (buggy) paths — for comparison only
V1_DIR=$UNZ/filtered_clean

# v2 (corrected) output — fits in /home/copilot (~22GB free, need ~10GB)
V2_BASE=$HOME/v2
V2_PERCHR=$V2_BASE/perchr
V2_PLINK=$V2_BASE/plink

# Tools
export PATH=/staging/conda/envs/bioinfo/bin:$PATH

# ── Helper ───────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die() { log "FATAL: $*"; exit 1; }

###############################################################################
# PHASE 0: Preflight checks
###############################################################################
log "===== PHASE 0: Preflight checks ====="

command -v bcftools >/dev/null || die "bcftools not found"
command -v plink2   >/dev/null || die "plink2 not found"
log "bcftools: $(bcftools --version | head -1)"
log "plink2:   $(plink2 --version | head -1)"

[ -f "$HQ_POS_FILE" ] || die "HQ position file not found: $HQ_POS_FILE"
[ -f "$HET_EXCLUDE" ] || die "Het outlier file not found: $HET_EXCLUDE"
[ -f "$REF_GENOME"  ] || die "Reference genome not found: $REF_GENOME"

for chr in $(seq 1 22); do
  [ -f "$UNZ/chr${chr}.dose.vcf.gz" ] || die "Missing dose VCF: chr${chr}"
done
log "All 22 dose VCFs present"

hq_count=$(wc -l < "$HQ_POS_FILE")
log "HQ position file: $hq_count variants"

avail_gb=$(df --output=avail "$HOME" | tail -1 | awk '{print int($1/1048576)}')
log "Available disk in HOME: ${avail_gb} GB"
[ "$avail_gb" -gt 8 ] || die "Insufficient disk (need ~10GB, have ${avail_gb}GB)"

mkdir -p "$V2_PERCHR" "$V2_PLINK"
log "Output directory: $V2_BASE"

###############################################################################
# PHASE 1: Create position-based targets file
###############################################################################
log ""
log "===== PHASE 1: Create targets file ====="

TARGETS=$V2_BASE/HQ_targets_posrefalt.tsv

if [ -f "$TARGETS" ]; then
  log "Targets file already exists, skipping creation"
else
  # bcftools --targets-file format: CHROM\tPOS\tREF,ALT (alleles comma-separated in col 3)
  # Source file has: CHROM\tPOS\tREF\tALT\tID (5 tab-separated columns)
  log "Extracting positions from $HQ_POS_FILE ..."
  awk -F'\t' '{print $1"\t"$2"\t"$3","$4}' "$HQ_POS_FILE" > "$TARGETS"

  target_count=$(wc -l < "$TARGETS")
  log "Targets file created: $target_count entries (3-column: CHROM POS REF,ALT)"
  [ "$target_count" -eq "$hq_count" ] || die "Targets count mismatch: $target_count vs $hq_count"
fi

log "Targets file sample:"
head -3 "$TARGETS"
echo "..."
tail -3 "$TARGETS"

###############################################################################
# PHASE 2: Per-chromosome VCF filter → PLINK conversion (space-efficient)
#
# For each chromosome:
#   1. bcftools view (filter) → temp VCF
#   2. plink2 --vcf → per-chr .pgen/.pvar/.psam
#   3. Delete temp VCF immediately
# Peak disk: one temp VCF (~2-3GB) + accumulated .pgen files (~3GB total)
###############################################################################
log ""
log "===== PHASE 2: Per-chromosome filter + convert ====="

total_variants=0
MERGE_LIST=$V2_PERCHR/merge_list.txt
> "$MERGE_LIST"

for chr in $(seq 1 22); do
  pgen=$V2_PERCHR/chr${chr}
  tmpvcf=$V2_PERCHR/chr${chr}_tmp.vcf.gz

  if [ -f "${pgen}.pgen" ]; then
    log "chr${chr}: pgen exists, skipping"
    n=$(wc -l < "${pgen}.pvar" | awk '{print $1-1}')  # subtract header line
    total_variants=$((total_variants + n))
    echo "$pgen" >> "$MERGE_LIST"
    continue
  fi

  log "chr${chr}: filtering + converting..."
  t0=$(date +%s)

  # Step 1: filter VCF (position-based, exclude het outliers)
  bcftools view \
    "$UNZ/chr${chr}.dose.vcf.gz" \
    --samples-file "^${HET_EXCLUDE}" \
    --targets-file "$TARGETS" \
    --output-type z \
    --output "$tmpvcf" \
    --threads $THREADS

  # Step 2: convert to per-chr pgen
  plink2 --vcf "$tmpvcf" 'dosage=DS' \
    --double-id \
    --vcf-half-call missing \
    --make-pgen \
    --out "$pgen" \
    --threads $THREADS \
    --memory $MEM_MB

  # Step 3: delete temp VCF immediately to save space
  rm -f "$tmpvcf"

  t1=$(date +%s)
  elapsed=$((t1 - t0))

  # Count variants (pvar has 1 header line)
  n_variants=$(tail -n +2 "${pgen}.pvar" | wc -l)
  n_samples=$(wc -l < "${pgen}.psam" | awk '{print $1-1}')
  total_variants=$((total_variants + n_variants))

  log "chr${chr}: ${n_samples} samples, ${n_variants} variants (${elapsed}s)"
  echo "$pgen" >> "$MERGE_LIST"

  # Disk check after each chromosome
  avail_now=$(df --output=avail "$HOME" | tail -1 | awk '{print int($1/1048576)}')
  log "  disk remaining: ${avail_now} GB"
  [ "$avail_now" -gt 3 ] || die "Running low on disk (${avail_now}GB). Stopping."
done

log ""
log "PHASE 2 COMPLETE: total $total_variants variants across 22 chromosomes"
diff=$((total_variants - 10009530))
absdiff=${diff#-}
log "Difference from expected 10,009,530: $diff"
if [ "$absdiff" -gt 1000 ]; then
  log "WARNING: Difference > 1000 — investigate before proceeding!"
fi

###############################################################################
# PHASE 3: Merge chromosomes → single PLINK binary
###############################################################################
log ""
log "===== PHASE 3: Merge chromosomes + ref-from-fa ====="

V2_BED=$V2_PLINK/UZB_v2_clean

if [ -f "${V2_BED}.bed" ]; then
  log "Merged PLINK binary already exists, skipping"
else
  log "Merging 22 chromosomes with --pmerge-list..."
  t0=$(date +%s)

  # Take first chr as base, rest in merge list
  first=$(head -1 "$MERGE_LIST")
  tail -n +2 "$MERGE_LIST" > "${MERGE_LIST}.rest"

  plink2 --pfile "$first" \
    --pmerge-list "${MERGE_LIST}.rest" \
    --fa "$REF_GENOME" --ref-from-fa force \
    --make-bed \
    --out "$V2_BED" \
    --threads $THREADS \
    --memory $MEM_MB

  t1=$(date +%s)
  log "Merge done in $((t1 - t0))s"
fi

v2_variants=$(wc -l < "${V2_BED}.bim")
v2_samples=$(wc -l < "${V2_BED}.fam")
log "PLINK binary: $v2_samples samples, $v2_variants variants"

v2_dots=$(awk '$2=="."' "${V2_BED}.bim" | wc -l)
log "Dot-ID variants in v2 BIM: $v2_dots"

# Clean up per-chr pgen files to save space
log "Cleaning per-chromosome files..."
rm -rf "$V2_PERCHR"
log "Freed space from per-chr files"

###############################################################################
# PHASE 4: Standard QC
###############################################################################
log ""
log "===== PHASE 4: Standard QC ====="

V2_QC=$V2_PLINK/UZB_v2_qc

if [ -f "${V2_QC}.bed" ]; then
  log "QC dataset already exists, skipping"
else
  log "Applying --geno 0.05, --mind 0.05, --maf 0.01..."
  t0=$(date +%s)

  plink2 --bfile "$V2_BED" \
    --geno 0.05 \
    --mind 0.05 \
    --maf 0.01 \
    --make-bed \
    --out "$V2_QC" \
    --threads $THREADS \
    --memory $MEM_MB

  t1=$(date +%s)
  log "QC done in $((t1-t0))s"
fi

v2_qc_samples=$(wc -l < "${V2_QC}.fam")
v2_qc_variants=$(wc -l < "${V2_QC}.bim")
log "After QC: $v2_qc_samples samples, $v2_qc_variants variants"

# KEY comparison with v1
log ""
log "=== V1 vs V2 comparison ==="
v1_samples=$(wc -l < "$V1_DIR/UZB_imputed_HQ_qc.fam" 2>/dev/null || echo "?")
v1_variants=$(wc -l < "$V1_DIR/UZB_imputed_HQ_qc.bim" 2>/dev/null || echo "?")
log "V1 (buggy):    $v1_samples samples, $v1_variants variants"
log "V2 (corrected): $v2_qc_samples samples, $v2_qc_variants variants"

if [ "$v2_qc_samples" -gt "${v1_samples:-0}" ] 2>/dev/null; then
  rescued=$((v2_qc_samples - v1_samples))
  log "** $rescued samples RESCUED (no longer removed by --mind 0.05) **"
fi

###############################################################################
# PHASE 5: Unique variant IDs + LD pruning
###############################################################################
log ""
log "===== PHASE 5: Variant IDs + LD pruning ====="

# 6a. Assign CHR:POS:REF:ALT IDs
V2_UNIQUE=$V2_PLINK/UZB_v2_unique

if [ -f "${V2_UNIQUE}.bed" ]; then
  log "Unique-ID dataset already exists, skipping"
else
  log "Assigning CHR:POS:REF:ALT variant IDs..."
  plink2 --bfile "$V2_QC" \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 50 \
    --make-bed \
    --out "$V2_UNIQUE" \
    --threads $THREADS
fi

v2_unique_variants=$(wc -l < "${V2_UNIQUE}.bim")
log "Unique IDs assigned: $v2_unique_variants variants"

# 6b. LD pruning (same parameters as original Step 7)
V2_PRUNE=$V2_PLINK/UZB_v2_pruned

if [ -f "${V2_PRUNE}.prune.in" ]; then
  log "LD pruning already done, skipping"
else
  log "LD pruning (window=1000kb, step=1, r²=0.05)..."
  t0=$(date +%s)

  plink2 --bfile "$V2_UNIQUE" \
    --indep-pairwise 1000kb 1 0.05 \
    --out "$V2_PRUNE" \
    --threads $THREADS

  t1=$(date +%s)
  log "LD pruning done in $((t1-t0))s"
fi

prune_in=$(wc -l < "${V2_PRUNE}.prune.in")
prune_out=$(wc -l < "${V2_PRUNE}.prune.out")
log "LD-independent variants: $prune_in (removed: $prune_out)"

# 6c. Extract pruned set
V2_LDP=$V2_PLINK/UZB_v2_ldpruned

if [ -f "${V2_LDP}.bed" ]; then
  log "LD-pruned dataset already exists, skipping"
else
  log "Extracting LD-pruned variant set..."
  plink2 --bfile "$V2_UNIQUE" \
    --extract "${V2_PRUNE}.prune.in" \
    --make-bed \
    --out "$V2_LDP" \
    --threads $THREADS
fi

ldp_samples=$(wc -l < "${V2_LDP}.fam")
ldp_variants=$(wc -l < "${V2_LDP}.bim")
log "LD-pruned set: $ldp_samples samples, $ldp_variants variants"

###############################################################################
# PHASE 6: Validation & Summary
###############################################################################
log ""
log "================================================================="
log "               RECALCULATION COMPLETE — SUMMARY"
log "================================================================="
log ""
log "Output directory: $V2_BASE"
log ""
log "Files created:"
log "  PLINK (raw):    ${V2_BED}.{bed,bim,fam}"
log "  PLINK (QC):     ${V2_QC}.{bed,bim,fam}"
log "  PLINK (unique): ${V2_UNIQUE}.{bed,bim,fam}"
log "  PLINK (LD):     ${V2_LDP}.{bed,bim,fam}"
log ""
log "Pipeline metrics:"
log "  Pre-filter variants:    $v2_variants (should be ~10,009,530)"
log "  Pre-filter dot-IDs:     $v2_dots"
log "  Post-QC samples:        $v2_qc_samples (v1: ${v1_samples:-?})"
log "  Post-QC variants:       $v2_qc_variants (v1: ${v1_variants:-?})"
log "  LD-pruned variants:     $ldp_variants"
log ""
log "V1 (buggy) had:  10,846,569 pre-QC → 5,383,832 post-QC (1,062 samples)"
log "V2 (corrected):  $v2_variants pre-QC → $v2_qc_variants post-QC ($v2_qc_samples samples)"
log ""

# Validate: no dot IDs should be in the unique-ID set
remaining_dots=$(awk '$2=="."' "${V2_UNIQUE}.bim" | wc -l)
if [ "$remaining_dots" -eq 0 ]; then
  log "VALIDATION: No dot-IDs in final dataset ✓"
else
  log "VALIDATION WARNING: $remaining_dots dot-IDs remain in unique-ID dataset"
fi

log ""
log "NEXT STEPS (manual):"
log "  1. Review this log for errors/warnings"
log "  2. Compare V2 PCA with V1 PCA (local PCA with UZB_v2_ldpruned)"
log "  3. Re-run ADMIXTURE with UZB_v2_ldpruned"
log "  4. Re-run ROH/IBD with UZB_v2_qc"
log "  5. Update step documentation HTML files"
log ""
log "===== ALL DONE: $(date) ====="
