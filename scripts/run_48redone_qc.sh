#!/usr/bin/env bash
# =============================================================================
# run_48redone_qc.sh
# DRAGEN Array genotyping + PLINK QC for 48 rescanned samples
# (chips 209422280005 and 209422280032; GSA-24v3 platform)
#
# Mirrors the pipeline used for the original 95 GWAS2026 samples
# (dragen_array_test_20260513/gwas2026)
#
# Run as root on the DRAGEN server
# Usage:  bash run_48redone_qc.sh  [--skip-dragena]
# =============================================================================
set -euo pipefail

# ---------- Paths ------------------------------------------------------------
BASE=/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513
OUTDIR=${BASE}/gwas48redone
IDAT_DIR="/staging/GWAS2026/48redone/ALSU rescan"
DRAGENA=${BASE}/software/dragena-linux-x64-DAv1.4.0-rc2-sha.9b051a4be2c6ade41bfa804c3bccc41cdbdc1c5b/dragena/dragena
BPM=/staging/ALSU-analysis/Conversion/GSA-24v3-0_A2.bpm
CSV_MANIFEST=/staging/ALSU-analysis/Conversion/GSA-24v3-0_A2.csv
EGT=${BASE}/GSA-24v3-0_A1_ClusterFile.egt
GENOME_FA=/staging/Genomes/Human/ensembl/unpacked/Homo_sapiens.GRCh38.dna.primary_assembly.fa
PLINK2=/usr/local/bin/plink2
THREADS=16

SKIP_DRAGENA=0
[[ "${1:-}" == "--skip-dragena" ]] && SKIP_DRAGENA=1

LOG() { echo "[$(date '+%Y-%m-%dT%H:%M:%S')] $*"; }

# ---------- 1. Setup directories ---------------------------------------------
LOG "=== Creating output directories ==="
mkdir -p "${OUTDIR}"/{gtc,vcf,plink,qc,logs}

SAMPLE_SHEET=${OUTDIR}/GWAS48redone_dragen_sample_sheet.csv

# ---------- 2. Generate sample sheet -----------------------------------------
LOG "=== Generating DRAGEN sample sheet ==="
cat > "${SAMPLE_SHEET}" << 'EOF'
[Header],,,,,,,,,,
Institute Name,Illumina,Inc.,,,,,,,,
Investigator Name,Abdushukur Rakhmatullaev,,,,,,,,,
Project Name,GSA-24v3-ALSU_rescan_48_samples-18052026,,,,,,,,,
Date,2026-05-18,,,,,,,,,
,,,,,,,,,,
[Manifests],,,,,,,,,,
A,GSA-24v3-0_A2,,,,,,,,,
,,,,,,,,,,
[Data]
Sample_ID,SentrixBarcode_A,SentrixPosition_A,Sample_Plate,Sample_Well,Sample_Group,Gender,Sample_Name,Replicate,Parent1,Parent2
03-328,209422280032,R01C01,WG7120892-DNA,A1,,,,,,
03-72,209422280032,R02C01,WG7120892-DNA,A2,,,,,,
03-73,209422280032,R03C01,WG7120892-DNA,A3,,,,,,
03-85,209422280032,R04C01,WG7120892-DNA,A4,,,,,,
07-42,209422280032,R05C01,WG7120892-DNA,A5,,,,,,
07-79,209422280032,R06C01,WG7120892-DNA,A6,,,,,,
07-96,209422280032,R07C01,WG7120892-DNA,A7,,,,,,
08-290,209422280032,R08C01,WG7120892-DNA,A8,,,,,,
08-291,209422280032,R09C01,WG7120892-DNA,A9,,,,,,
08-292,209422280032,R10C01,WG7120892-DNA,A10,,,,,,
08-296,209422280032,R11C01,WG7120892-DNA,A11,,,,,,
08-297,209422280032,R12C01,WG7120892-DNA,A12,,,,,,
08-298,209422280032,R01C02,WG7120892-DNA,B1,,,,,,
08-305,209422280032,R02C02,WG7120892-DNA,B2,,,,,,
08-306,209422280032,R03C02,WG7120892-DNA,B3,,,,,,
08-309,209422280032,R04C02,WG7120892-DNA,B4,,,,,,
08-310,209422280032,R05C02,WG7120892-DNA,B5,,,,,,
08-311,209422280032,R06C02,WG7120892-DNA,B6,,,,,,
08-314,209422280032,R07C02,WG7120892-DNA,B7,,,,,,
08-317,209422280032,R08C02,WG7120892-DNA,B8,,,,,,
08-321,209422280032,R09C02,WG7120892-DNA,B9,,,,,,
08-322,209422280032,R10C02,WG7120892-DNA,B10,,,,,,
08-323,209422280032,R11C02,WG7120892-DNA,B11,,,,,,
08-326,209422280032,R12C02,WG7120892-DNA,B12,,,,,,
08-331,209422280005,R01C01,WG7120892-DNA,C1,,,,,,
08-332,209422280005,R02C01,WG7120892-DNA,C2,,,,,,
08-336,209422280005,R03C01,WG7120892-DNA,C3,,,,,,
08-341,209422280005,R04C01,WG7120892-DNA,C4,,,,,,
08-352,209422280005,R05C01,WG7120892-DNA,C5,,,,,,
08-356,209422280005,R06C01,WG7120892-DNA,C6,,,,,,
08-360,209422280005,R07C01,WG7120892-DNA,C7,,,,,,
08-374,209422280005,R08C01,WG7120892-DNA,C8,,,,,,
08-386,209422280005,R09C01,WG7120892-DNA,C9,,,,,,
08-390,209422280005,R10C01,WG7120892-DNA,C10,,,,,,
08-391,209422280005,R11C01,WG7120892-DNA,C11,,,,,,
08-392,209422280005,R12C01,WG7120892-DNA,C12,,,,,,
08-399,209422280005,R01C02,WG7120892-DNA,D1,,,,,,
08-407,209422280005,R02C02,WG7120892-DNA,D2,,,,,,
08-410,209422280005,R03C02,WG7120892-DNA,D3,,,,,,
08-411,209422280005,R04C02,WG7120892-DNA,D4,,,,,,
08-416,209422280005,R05C02,WG7120892-DNA,D5,,,,,,
08-419,209422280005,R06C02,WG7120892-DNA,D6,,,,,,
08-421,209422280005,R07C02,WG7120892-DNA,D7,,,,,,
08-430,209422280005,R08C02,WG7120892-DNA,D8,,,,,,
08-433,209422280005,R09C02,WG7120892-DNA,D9,,,,,,
08-487,209422280005,R10C02,WG7120892-DNA,D10,,,,,,
08-490,209422280005,R11C02,WG7120892-DNA,D11,,,,,,
08-505,209422280005,R12C02,WG7120892-DNA,D12,,,,,,
EOF
LOG "Sample sheet written to ${SAMPLE_SHEET}"

if [[ $SKIP_DRAGENA -eq 0 ]]; then

# ---------- 3. DRAGEN Array: IDAT → GTC --------------------------------------
LOG "=== Step 3: DRAGEN Array genotype calling (IDAT → GTC) ==="
"${DRAGENA}" genotype-call \
    --bpm-manifest-file "${BPM}" \
    --cluster-file "${EGT}" \
    --idat-folder "${IDAT_DIR}" \
    --idat-sample-sheet "${SAMPLE_SHEET}" \
    --output-folder "${OUTDIR}/gtc" \
    --gencall-score-cut-off 0.15 \
    --num-threads ${THREADS} \
    2>&1 | tee "${OUTDIR}/logs/genotype_call.log"
LOG "GTC files written to ${OUTDIR}/gtc"
LOG "GTC count: $(ls ${OUTDIR}/gtc/*.gtc 2>/dev/null | wc -l)"

# ---------- 4. DRAGEN Array: GTC → per-sample VCF ----------------------------
LOG "=== Step 4: DRAGEN Array GTC → VCF ==="
"${DRAGENA}" gtc-to-vcf \
    --bpm-manifest-file "${BPM}" \
    --csv-manifest-file "${CSV_MANIFEST}" \
    --genome-fasta-file "${GENOME_FA}" \
    --gtc-folder "${OUTDIR}/gtc" \
    --gtc-sample-sheet "${SAMPLE_SHEET}" \
    --output-folder "${OUTDIR}/vcf" \
    --bgzip \
    2>&1 | tee "${OUTDIR}/logs/gtc_to_vcf.log"
LOG "VCF files written to ${OUTDIR}/vcf"
LOG "VCF count: $(ls ${OUTDIR}/vcf/*.vcf.gz 2>/dev/null | wc -l)"

else
    LOG "Skipping DRAGEN Array steps (--skip-dragena)"
fi

# ---------- 5. Merge per-sample VCFs -----------------------------------------
LOG "=== Step 5: Merging per-sample VCFs ==="
VCF_LIST=${OUTDIR}/vcf_list.txt
ls "${OUTDIR}/vcf/"*.snv.vcf.gz > "${VCF_LIST}"
VCF_COUNT=$(wc -l < "${VCF_LIST}")
LOG "Found ${VCF_COUNT} VCFs to merge"

# bcftools merge requires tabix-indexed inputs — they should already be
# produced by DRAGEN Array (*.snv.vcf.gz.tbi).  If not, index them.
while IFS= read -r vcf; do
    [[ -f "${vcf}.tbi" ]] || tabix -p vcf "${vcf}"
done < "${VCF_LIST}"

bcftools merge \
    --file-list "${VCF_LIST}" \
    --output-type z \
    --output "${OUTDIR}/gwas48redone_merged.vcf.gz" \
    --threads ${THREADS}
tabix -p vcf "${OUTDIR}/gwas48redone_merged.vcf.gz"
LOG "Merged VCF: ${OUTDIR}/gwas48redone_merged.vcf.gz"

# ---------- 6. VCF → PLINK2 BED (raw) ----------------------------------------
LOG "=== Step 6: VCF → PLINK2 BED (raw) ==="
${PLINK2} \
    --vcf "${OUTDIR}/gwas48redone_merged.vcf.gz" \
    --double-id \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 2000 missing \
    --max-alleles 2 \
    --fa "${GENOME_FA}" \
    --ref-from-fa force \
    --not-chr X Y MT \
    --vcf-half-call missing \
    --make-bed \
    --out "${OUTDIR}/plink/gwas48redone_raw" \
    --threads ${THREADS} \
    2>&1 | tee "${OUTDIR}/logs/plink_import.log"
LOG "Raw PLINK: $(wc -l < ${OUTDIR}/plink/gwas48redone_raw.fam) samples"

# ---------- 7. Sample QC: missingness ----------------------------------------
LOG "=== Step 7: Sample QC — missingness ==="
${PLINK2} \
    --bfile "${OUTDIR}/plink/gwas48redone_raw" \
    --missing \
    --out "${OUTDIR}/qc/gwas48redone_miss" \
    --threads ${THREADS}

# Flag samples with F_MISS > 0.20
python3 - <<'PYEOF'
import sys
miss_file = "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/gwas48redone/qc/gwas48redone_miss.smiss"
remove = []
with open(miss_file) as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split()
        if float(parts[4]) > 0.20:
            remove.append(f"{parts[0]}\t{parts[1]}")
print(f"Samples with F_MISS > 0.20: {len(remove)}")
out = "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/gwas48redone/qc/remove_miss.txt"
with open(out, 'w') as f:
    if remove:
        f.write('\n'.join(remove) + '\n')
    else:
        f.write('#FID\tIID\n')
PYEOF
LOG "Missingness > 20%: $(grep -v '^#' ${OUTDIR}/qc/remove_miss.txt | wc -l) samples"

# ---------- 8. Sample QC: heterozygosity -------------------------------------
LOG "=== Step 8: Sample QC — heterozygosity ==="
${PLINK2} \
    --bfile "${OUTDIR}/plink/gwas48redone_raw" \
    --indep-pairwise 50 5 0.1 \
    --chr 1-22 \
    --out "${OUTDIR}/qc/prune_ibd" \
    --threads ${THREADS}

${PLINK2} \
    --bfile "${OUTDIR}/plink/gwas48redone_raw" \
    --extract "${OUTDIR}/qc/prune_ibd.prune.in" \
    --het \
    --out "${OUTDIR}/qc/gwas48redone_het" \
    --threads ${THREADS}

python3 - <<'PYEOF'
import statistics, sys
het_file = "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/gwas48redone/qc/gwas48redone_het.het"
remove = []
vals = []
rows = []
with open(het_file) as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split()
        # PLINK2 het columns: #FID IID O(HOM) E(HOM) OBS_CT F
        het_f = float(parts[5])
        vals.append(het_f)
        rows.append((parts[0], parts[1], het_f))
if vals:
    mn = statistics.mean(vals)
    sd = statistics.stdev(vals)
    lo, hi = mn - 3*sd, mn + 3*sd
    for fid, iid, f in rows:
        if f < lo or f > hi:
            remove.append(f"{fid}\t{iid}")
    print(f"Het F: mean={mn:.4f} sd={sd:.4f} cutoffs=[{lo:.4f},{hi:.4f}]")
    print(f"Samples failing het ±3SD: {len(remove)}")
out = "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/gwas48redone/qc/remove_het.txt"
with open(out, 'w') as f:
    if remove:
        f.write('\n'.join(remove) + '\n')
    else:
        f.write('#FID\tIID\n')
PYEOF

# ---------- 9. Sample QC: KING relatedness -----------------------------------
LOG "=== Step 9: Sample QC — KING-robust relatedness ==="
${PLINK2} \
    --bfile "${OUTDIR}/plink/gwas48redone_raw" \
    --extract "${OUTDIR}/qc/prune_ibd.prune.in" \
    --king-cutoff 0.185 \
    --out "${OUTDIR}/qc/gwas48redone_king" \
    --threads ${THREADS}

KING_OUT=$(wc -l < "${OUTDIR}/qc/gwas48redone_king.king.cutoff.out.id" 2>/dev/null || echo 0)
LOG "KING: ${KING_OUT} samples removed at cutoff 0.185"

# ---------- 10. Combine exclusion lists and apply ----------------------------
LOG "=== Step 10: Combining removal lists ==="
python3 - <<'PYEOF'
base = "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/gwas48redone/qc"
seen = set()
all_remove = []
for fname in ["remove_miss.txt", "remove_het.txt"]:
    try:
        with open(f"{base}/{fname}") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'): continue
                key = tuple(line.split())
                if key not in seen:
                    seen.add(key)
                    all_remove.append(line)
    except FileNotFoundError:
        pass
# KING exclusion (plink2 format: sample-id-only tsv with header)
king_out = f"{base}/gwas48redone_king.king.cutoff.out.id"
try:
    with open(king_out) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                key = (parts[0], parts[1])
                if key not in seen:
                    seen.add(key)
                    all_remove.append(f"{parts[0]}\t{parts[1]}")
except FileNotFoundError:
    pass
with open(f"{base}/remove_samples_all.txt", 'w') as f:
    f.write('#FID\tIID\n')
    for r in all_remove:
        f.write(r + '\n')
print(f"Total samples to remove: {len(all_remove)}")
PYEOF

${PLINK2} \
    --bfile "${OUTDIR}/plink/gwas48redone_raw" \
    --remove <(grep -v '^#' "${OUTDIR}/qc/remove_samples_all.txt" || true) \
    --make-bed \
    --out "${OUTDIR}/plink/gwas48redone_sampleqc" \
    --threads ${THREADS}
LOG "After sample QC: $(wc -l < ${OUTDIR}/plink/gwas48redone_sampleqc.fam) samples"

# ---------- 11. Variant QC ---------------------------------------------------
LOG "=== Step 11: Variant QC (geno 10%, HWE 1e-6, MAF 0.001) ==="
${PLINK2} \
    --bfile "${OUTDIR}/plink/gwas48redone_sampleqc" \
    --geno 0.10 \
    --hwe 1e-6 \
    --maf 0.001 \
    --make-bed \
    --out "${OUTDIR}/plink/gwas48redone_varqc" \
    --threads ${THREADS} \
    2>&1 | tee "${OUTDIR}/logs/varqc.log"

FINAL_VARS=$(wc -l < "${OUTDIR}/plink/gwas48redone_varqc.bim")
FINAL_SAMP=$(wc -l < "${OUTDIR}/plink/gwas48redone_varqc.fam")
LOG "=== QC COMPLETE ==="
LOG "Final: ${FINAL_SAMP} samples, ${FINAL_VARS} variants"
LOG "Output PLINK: ${OUTDIR}/plink/gwas48redone_varqc.{bed,bim,fam}"
