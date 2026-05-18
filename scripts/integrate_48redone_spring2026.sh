#!/usr/bin/env bash
# =============================================================================
# integrate_48redone_spring2026.sh
# Integrate the 48 rescanned samples into the spring2026 analysis
# Scope: QC + PCA + ADMIXTURE (imputation/GWAS deferred)
#
# Prerequisites:
#   - run_48redone_qc.sh completed successfully
#   - gwas48redone_sampleqc.{bed,bim,fam} exists
#
# Strategy:
#   The 48 samples use the same GSA-24v3 DRAGEN platform as the original
#   95 GWAS2026 samples.  At common variant positions (471K sites shared
#   between ConvSK and DRAGEN), we can extend gwas2026_common_nomiss with
#   the 48 new samples, then re-merge with alsu_common_nomiss.
#
# Run as root on the DRAGEN server
# Usage:  bash integrate_48redone_spring2026.sh
# =============================================================================
set -euo pipefail

BASE=/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513
CROSS=${BASE}/cross_gwas2026_alsu_winter
REDONE=${BASE}/gwas48redone
OUTDIR=${BASE}/cross_gwas48_alsu_winter   # new integration output dir
SPRING=/staging/ALSU-analysis/spring2026
PLINK2=/usr/local/bin/plink2
PLINK1=/usr/local/bin/plink
THREADS=16

LOG() { echo "[$(date '+%Y-%m-%dT%H:%M:%S')] $*"; }
VERIFY_FAM() { echo "  -> $(wc -l < $1) samples in $1"; }

mkdir -p "${OUTDIR}"/{pca,admixture,logs}

# ---------- Verify inputs exist -----------------------------------------------
LOG "=== Verifying inputs ==="
for f in \
    "${REDONE}/plink/gwas48redone_sampleqc.bed" \
    "${CROSS}/gwas2026_sampleqc.bed" \
    "${CROSS}/common_pos.ids" \
    "${CROSS}/alsu_common_nomiss.bed"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Missing required file: $f"
    echo "Check:"
    echo "  - run_48redone_qc.sh completed (for gwas48redone_sampleqc)"
    echo "  - cross_gwas2026_alsu_winter step was already run (for alsu_common_nomiss)"
    exit 1
  fi
done
LOG "All inputs present"

# Check how many 48redone samples passed QC
N48=$(wc -l < "${REDONE}/plink/gwas48redone_sampleqc.fam")
LOG "gwas48redone after sample QC: ${N48} samples"
if [[ ${N48} -eq 0 ]]; then
  echo "ERROR: No samples passed QC in gwas48redone. Check run_48redone_qc.sh output."
  exit 1
fi

# ---------- Step 1: Extract common positions from gwas48redone ---------------
LOG "=== Step 1: Extracting common positions from gwas48redone_sampleqc ==="
# common_pos.ids uses CHR:POS format matching gwas2026 variant IDs
${PLINK2} \
    --bfile "${REDONE}/plink/gwas48redone_sampleqc" \
    --extract "${CROSS}/common_pos.ids" \
    --make-bed \
    --out "${OUTDIR}/gwas48redone_common" \
    --threads ${THREADS}
VERIFY_FAM "${OUTDIR}/gwas48redone_common.fam"

# ---------- Step 2: Apply missingness filter (same as original gwas2026) -----
LOG "=== Step 2: Missingness filter on gwas48redone_common ==="
# Original gwas2026_common_nomiss used geno 0.05 (max 5% per-variant missing)
${PLINK2} \
    --bfile "${OUTDIR}/gwas48redone_common" \
    --geno 0.05 \
    --make-bed \
    --out "${OUTDIR}/gwas48redone_common_nomiss" \
    --threads ${THREADS}
VERIFY_FAM "${OUTDIR}/gwas48redone_common_nomiss.fam"

# Also filter variants to match those in gwas2026_common_nomiss
# (take the intersection to ensure same variant set for merging)
cut -f2 "${CROSS}/gwas2026_common_nomiss.bim" > "${OUTDIR}/gwas2026_common_nomiss_vars.txt"
${PLINK2} \
    --bfile "${OUTDIR}/gwas48redone_common_nomiss" \
    --extract "${OUTDIR}/gwas2026_common_nomiss_vars.txt" \
    --make-bed \
    --out "${OUTDIR}/gwas48redone_for_merge" \
    --threads ${THREADS}
VERIFY_FAM "${OUTDIR}/gwas48redone_for_merge.fam"

# ---------- Step 3: Combine gwas2026 + gwas48redone (DRAGEN combined) --------
LOG "=== Step 3: Merging gwas2026_common_nomiss + gwas48redone_for_merge ==="
# Update IDs in gwas48redone to use GWAS2026_XXXX scheme (continuing from 0095)
# Determine max existing GWAS2026 ID number
python3 - <<'PYEOF'
import sys, os
base = "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513"
cross = f"{base}/cross_gwas2026_alsu_winter"
out = f"{base}/cross_gwas48_alsu_winter"

# Read existing GWAS2026 ID mapping from sample_map.tsv
max_num = 0
with open(f"{cross}/sample_map.tsv") as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if parts[0] == 'GWAS2026':
            try:
                n = int(parts[2].replace('GWAS2026_',''))
                if n > max_num: max_num = n
            except: pass

# Read 48redone FAM and create ID update file
fam48_path = f"{base}/gwas48redone/plink/gwas48redone_sampleqc.fam"
update_lines = []
map_lines = []
n = max_num
with open(fam48_path) as f:
    for line in f:
        parts = line.strip().split()
        old_fid, old_iid = parts[0], parts[1]
        n += 1
        new_id = f"GWAS2026_{n:04d}"
        update_lines.append(f"{old_fid}\t{old_iid}\t{new_id}\t{new_id}")
        map_lines.append(f"GWAS2026\tGWAS2026\t{new_id}\t{old_fid}\t{old_iid}")

with open(f"{out}/gwas48redone_update_ids.txt", 'w') as f:
    f.write('\n'.join(update_lines) + '\n')
with open(f"{out}/gwas48redone_sample_map.tsv", 'w') as f:
    f.write('\n'.join(map_lines) + '\n')
print(f"Created ID update for {len(update_lines)} samples (starting at GWAS2026_{max_num+1:04d})")
PYEOF

${PLINK2} \
    --bfile "${OUTDIR}/gwas48redone_for_merge" \
    --update-ids "${OUTDIR}/gwas48redone_update_ids.txt" \
    --make-bed \
    --out "${OUTDIR}/gwas48redone_renamed" \
    --threads ${THREADS}

# Merge GWAS2026 (95 samples at common nomiss) + gwas48redone (N renamed samples)
# Use plink1 merge for BED merging (plink2 merge requires different syntax)
${PLINK1} \
    --bfile "${CROSS}/gwas2026_common_nomiss" \
    --bmerge "${OUTDIR}/gwas48redone_renamed" \
    --allow-no-sex \
    --out "${OUTDIR}/gwas_dragen_combined" \
    2>&1 | tee "${OUTDIR}/logs/merge_dragen.log"
VERIFY_FAM "${OUTDIR}/gwas_dragen_combined.fam"

# ---------- Step 4: Re-merge ALSU (ConvSK) + combined DRAGEN -----------------
LOG "=== Step 4: Re-merging alsu_common_nomiss + gwas_dragen_combined ==="
${PLINK1} \
    --bfile "${CROSS}/alsu_common_nomiss" \
    --bmerge "${OUTDIR}/gwas_dragen_combined" \
    --allow-no-sex \
    --out "${OUTDIR}/merged_with_48" \
    2>&1 | tee "${OUTDIR}/logs/merge_final.log"

MERGED_N=$(wc -l < "${OUTDIR}/merged_with_48.fam")
MERGED_V=$(wc -l < "${OUTDIR}/merged_with_48.bim")
LOG "Merged dataset: ${MERGED_N} samples, ${MERGED_V} variants"

# Append the 48redone sample entries to the full sample map
cp "${CROSS}/sample_map.tsv" "${OUTDIR}/sample_map.tsv"
cat "${OUTDIR}/gwas48redone_sample_map.tsv" >> "${OUTDIR}/sample_map.tsv"
LOG "Updated sample_map.tsv: $(wc -l < ${OUTDIR}/sample_map.tsv) entries (incl header)"

# ---------- Step 5: PCA prep (LD prune, MAF filter) --------------------------
LOG "=== Step 5: PCA — LD pruning and MAF filter ==="
${PLINK2} \
    --bfile "${OUTDIR}/merged_with_48" \
    --maf 0.01 \
    --geno 0.05 \
    --chr 1-22 \
    --indep-pairwise 200 50 0.25 \
    --out "${OUTDIR}/pca/prune" \
    --threads ${THREADS}

${PLINK2} \
    --bfile "${OUTDIR}/merged_with_48" \
    --extract "${OUTDIR}/pca/prune.prune.in" \
    --maf 0.01 \
    --geno 0.05 \
    --make-bed \
    --out "${OUTDIR}/pca/merged_ldpruned" \
    --threads ${THREADS}
VERIFY_FAM "${OUTDIR}/pca/merged_ldpruned.fam"

# ---------- Step 6: Run PCA --------------------------------------------------
LOG "=== Step 6: PCA (20 PCs) ==="
${PLINK2} \
    --bfile "${OUTDIR}/pca/merged_ldpruned" \
    --pca 20 \
    --out "${OUTDIR}/pca/merged_pca" \
    --threads ${THREADS} \
    2>&1 | tee "${OUTDIR}/logs/pca.log"
LOG "PCA eigenvec: ${OUTDIR}/pca/merged_pca.eigenvec"

# ---------- Step 7: ADMIXTURE prep -------------------------------------------
LOG "=== Step 7: ADMIXTURE prep ==="
# Use LD-pruned dataset with MAF filter
${PLINK2} \
    --bfile "${OUTDIR}/pca/merged_ldpruned" \
    --make-bed \
    --out "${OUTDIR}/admixture/uzb_for_admixture" \
    --threads ${THREADS}
VERIFY_FAM "${OUTDIR}/admixture/uzb_for_admixture.fam"

LOG ""
LOG "=== INTEGRATION COMPLETE ==="
LOG "New merged dataset: ${OUTDIR}/merged_with_48.{bed,bim,fam}"
LOG "  -> ${MERGED_N} samples (previous: 1193; added ~${N48} rescanned)"
LOG "PCA:       ${OUTDIR}/pca/merged_pca.eigenvec"
LOG "ADMIXTURE: ${OUTDIR}/admixture/uzb_for_admixture.{bed,bim,fam}"
LOG ""
LOG "NEXT STEPS:"
LOG "  1. Run ADMIXTURE (K=2..8) on uzb_for_admixture"
LOG "  2. Update spring2026 HTML pages with new sample counts and PCA plots"
LOG "  3. Decide on imputation scope (adding 48 to Michigan submission)"
