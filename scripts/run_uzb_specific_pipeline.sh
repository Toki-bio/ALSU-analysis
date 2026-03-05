#!/bin/bash
# ============================================================================
# Master script: Run the complete Uzbek-specific SNP identification pipeline
#
# This script:
#   1. Runs 01_extract_multipop.sh (extract EAS/SAS/AFR, compute Fst & freqs)
#   2. Creates symlink to existing UZB-EUR Fst
#   3. Runs 02_calculate_pbs.py (PBS, delta-AF, Uzbek-specific identification)
#
# Usage:
#   nohup bash run_uzb_specific_pipeline.sh > pipeline.log 2>&1 &
# ============================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="/staging/ALSU-analysis/admixture_analysis/pop_specific"
FST_DIR="/staging/ALSU-analysis/Fst_analysis"

echo "===== Uzbek-Specific SNP Pipeline ====="
echo "Start time: $(date)"
echo "Script directory: ${SCRIPT_DIR}"
echo "Output directory: ${OUTDIR}"
echo ""

# ---- Phase 1: Extract populations and compute pairwise Fst ----
echo "=========================================="
echo "PHASE 1: Multi-Population Extraction & Fst"
echo "=========================================="
bash "${SCRIPT_DIR}/01_extract_multipop.sh"

# ---- Symlink existing UZB-EUR Fst to output directory ----
echo ""
echo "Creating symlink to existing UZB-EUR Fst..."
if [ ! -f "${OUTDIR}/fst_UZB_vs_EUR.fst" ]; then
    ln -s "${FST_DIR}/genomewide_uzbek_vs_eur_fst.fst" "${OUTDIR}/fst_UZB_vs_EUR.fst"
fi
echo "  ✓ fst_UZB_vs_EUR.fst → genomewide_uzbek_vs_eur_fst.fst"

# ---- Phase 2: PBS & Uzbek-specific SNP analysis ----
echo ""
echo "=========================================="
echo "PHASE 2: PBS & Uzbek-Specific Analysis"
echo "=========================================="
python3 "${SCRIPT_DIR}/02_calculate_pbs.py" \
    --outdir "${OUTDIR}" \
    --pbs-threshold 0.3 \
    --delta-af-threshold 0.3

echo ""
echo "===== Pipeline Complete! ====="
echo "End time: $(date)"
echo ""
echo "Results in: ${OUTDIR}/"
echo "  pbs_results.tsv"
echo "  uzbek_specific_snps.tsv"
echo "  delta_af_all.tsv"
echo "  near_private_variants.tsv"
echo "  summary_stats.txt"
echo "  freq_*.frq (per-population allele frequencies)"
echo "  fst_*_vs_*.fst (pairwise Fst)"
