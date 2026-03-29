#!/bin/bash
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH

# ================================================================
# V2 FST Recalculation + Tool Version Capture
# ================================================================
# Uses the V2 global merged dataset (global_v2_admix) with 3,595 samples
# (1,047 UZB + 2,548 1000G) × 77,111 LD-pruned SNPs.
# Computes all 10 pairwise weighted Fst values.
# ================================================================

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

GLOBAL=~/v2/global_admixture/global_v2_admix
POP_MAP=~/v2/pca/pop_mapping.txt
OUT=~/v2/fst_v2
T=8

log "===== V2 FST RECALCULATION START ====="

# ── 0. Preflight ──
command -v plink >/dev/null || { echo "FATAL: plink not found"; exit 1; }
[ -f "${GLOBAL}.bed" ] || { echo "FATAL: ${GLOBAL}.bed not found"; exit 1; }
[ -f "$POP_MAP" ] || { echo "FATAL: pop_mapping not found at $POP_MAP"; exit 1; }

mkdir -p $OUT

# ── 1. Capture tool versions ──
log "=== TOOL VERSIONS ==="
log "plink1.9: $(plink --version 2>&1 | head -1 || echo 'N/A')"
log "plink2:   $(plink2 --version 2>&1 | head -1 || echo 'N/A')"
log "bcftools: $(bcftools --version 2>&1 | head -1 || echo 'N/A')"
log "vcftools: $(vcftools --version 2>&1 | head -1 || echo 'N/A')"
log "admixture: $(admixture 2>&1 | grep -i version | head -1 || echo 'N/A')"
log "R version: $(R --version 2>&1 | head -1 || echo 'N/A')"
# LEA package version
log "LEA version: $(Rscript -e 'cat(as.character(packageVersion("LEA")))' 2>/dev/null || echo 'N/A')"
log "=== END VERSIONS ==="

# ── 2. Build population labels from .fam + pop_mapping ──
# pop_mapping.txt format: SampleID\tSuperPop (for 1000G samples)
# UZB samples in .fam are NOT in pop_mapping → assign "UZB"
log "Building population assignments..."

awk -v pm="$POP_MAP" '
BEGIN { while ((getline line < pm) > 0) { split(line, a, "\t"); pop[a[1]] = a[2] } }
{
    fid = $1; iid = $2
    if (iid in pop) p = pop[iid]
    else if (fid in pop) p = pop[fid]
    else p = "UZB"
    print fid, iid, p
}
' ${GLOBAL}.fam > $OUT/all_pop_labels.txt

# Verify counts
log "Population counts:"
awk '{print $3}' $OUT/all_pop_labels.txt | sort | uniq -c | sort -rn

N_TOTAL=$(wc -l < $OUT/all_pop_labels.txt)
N_FAM=$(wc -l < ${GLOBAL}.fam)
log "Label file: $N_TOTAL lines, FAM file: $N_FAM lines"
[ "$N_TOTAL" -eq "$N_FAM" ] || { echo "FATAL: line count mismatch"; exit 1; }

# ── 3. Compute all 10 pairwise Fst ──
POPS="UZB EUR EAS SAS AFR"
log "Computing pairwise Fst for: $POPS"

# Results header
echo -e "Pop1\tPop2\tMeanFst\tWeightedFst\tValidMarkers" > $OUT/fst_v2_results.tsv

for pop1 in $POPS; do
  for pop2 in $POPS; do
    # Only compute upper triangle
    [ "$pop1" \< "$pop2" ] || [ "$pop1" = "$pop2" ] && continue

    log "  Computing Fst: $pop1 vs $pop2"

    # Create --keep file (FID IID for both pops)
    awk -v p1="$pop1" -v p2="$pop2" '$3 == p1 || $3 == p2 {print $1, $2}' \
      $OUT/all_pop_labels.txt > $OUT/keep_${pop1}_${pop2}.txt

    # Create --within file (FID IID CLUSTER)
    awk -v p1="$pop1" -v p2="$pop2" '$3 == p1 || $3 == p2 {print $1, $2, $3}' \
      $OUT/all_pop_labels.txt > $OUT/within_${pop1}_${pop2}.txt

    n1=$(awk -v p="$pop1" '$3==p' $OUT/within_${pop1}_${pop2}.txt | wc -l)
    n2=$(awk -v p="$pop2" '$3==p' $OUT/within_${pop1}_${pop2}.txt | wc -l)
    log "    $pop1: $n1 samples, $pop2: $n2 samples"

    # Run Fst
    plink --bfile $GLOBAL \
      --keep $OUT/keep_${pop1}_${pop2}.txt \
      --within $OUT/within_${pop1}_${pop2}.txt \
      --fst \
      --out $OUT/fst_${pop1}_vs_${pop2} \
      --silent

    # Parse results from .log
    mean_fst=$(grep "Mean Fst estimate:" $OUT/fst_${pop1}_vs_${pop2}.log | awk '{print $NF}')
    wt_fst=$(grep "Weighted Fst estimate:" $OUT/fst_${pop1}_vs_${pop2}.log | awk '{print $NF}')
    # Count valid markers from .fst file
    valid=$(awk 'NR>1 && $5 != "nan" && $5 != "NA"' $OUT/fst_${pop1}_vs_${pop2}.fst 2>/dev/null | wc -l)

    log "    Mean Fst=$mean_fst, Weighted Fst=$wt_fst, Valid markers=$valid"
    echo -e "${pop1}\t${pop2}\t${mean_fst}\t${wt_fst}\t${valid}" >> $OUT/fst_v2_results.tsv
  done
done

log "=== FST V2 RESULTS ==="
cat $OUT/fst_v2_results.tsv
log "=== END FST V2 ==="

# ── 4. Package results for download ──
log "Packaging results..."
cd $OUT
# Base64-encode the summary for easy download
base64 fst_v2_results.tsv > fst_v2_results.b64
# Also encode the full per-SNP Fst for UZB-EUR (largest interest)
base64 fst_UZB_vs_EUR.fst > fst_UZB_vs_EUR.b64 2>/dev/null || true

log "===== V2 FST RECALCULATION COMPLETE ====="
log "Output: $OUT/fst_v2_results.tsv"
log "Files:"
ls -lh $OUT/fst_v2_results.* $OUT/fst_*_vs_*.log 2>/dev/null
