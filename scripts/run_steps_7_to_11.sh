#!/bin/bash
# Spring 2026: Steps 7–11 (Local PCA → Global PCA → Fst → PBS → ADMIXTURE)
# Prereq: Step 6 completed — UZB_imputed_HQ_clean.{bed,bim,fam} in POSTDIR
# Run on DRAGEN: nohup bash /tmp/run_steps_7_to_11.sh > /tmp/step7_11_out.txt 2>&1 &
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH

# ============================
# Configuration
# ============================
WORKDIR=/staging/ALSU-analysis/spring2026
POSTDIR=${WORKDIR}/post_imputation
PCADIR=${WORKDIR}/pca
GLOBALDIR=${WORKDIR}/global_pca
FSTDIR=${WORKDIR}/fst
PBSDIR=${WORKDIR}/pbs
ADMIXDIR=${WORKDIR}/admixture

# Pre-built 1000G reference (ALL populations: AFR=671, EUR=522, EAS=515, SAS=492, AMR=348 = 2548 total)
KG_REF=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/1000g_grch38/KG_reference_final

T=8; M=200000  # threads, memory (KB)
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ============================
# Pre-flight checks
# ============================
log "========== SPRING 2026: STEPS 7–11 PIPELINE =========="
log "Checking prerequisites..."

for f in "${POSTDIR}/UZB_imputed_HQ_clean.bed" \
         "${POSTDIR}/UZB_imputed_HQ_clean.bim" \
         "${POSTDIR}/UZB_imputed_HQ_clean.fam" \
         "${KG_REF}.bed" "${KG_REF}.bim" "${KG_REF}.fam"; do
    if [ ! -f "$f" ]; then
        log "FATAL: Missing required file: $f"
        exit 1
    fi
done

NSAMPLES_IN=$(wc -l < ${POSTDIR}/UZB_imputed_HQ_clean.fam)
NVARS_IN=$(wc -l < ${POSTDIR}/UZB_imputed_HQ_clean.bim)
NKG=$(wc -l < ${KG_REF}.fam)
log "Step 6 output: ${NSAMPLES_IN} samples × ${NVARS_IN} variants"
log "1000G reference: ${NKG} samples"
log "Disk free: $(df -h /staging/ | awk 'NR==2{print $4}')"

mkdir -p "$PCADIR" "$GLOBALDIR" "$FSTDIR" "$PBSDIR" "$ADMIXDIR"

# ================================================================
# STEP 7: Local PCA (Uzbek-only)
# ================================================================
log "========== STEP 7: Local PCA =========="
t0=$(date +%s)

# 7a: QC filters
log "  7a: QC (geno 0.05, mind 0.05, maf 0.01)"
plink2 --bfile ${POSTDIR}/UZB_imputed_HQ_clean \
  --geno 0.05 --mind 0.05 --maf 0.01 \
  --make-bed --out ${PCADIR}/UZB_qc \
  --threads $T --memory $M
log "  After QC: $(wc -l < ${PCADIR}/UZB_qc.fam) samples × $(wc -l < ${PCADIR}/UZB_qc.bim) variants"

# 7b: Unique variant IDs (CHR:POS:REF:ALT) — use --output-chr 26 to ensure numeric chromosomes
log "  7b: Assign unique variant IDs"
plink2 --bfile ${PCADIR}/UZB_qc \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 50 \
  --output-chr 26 \
  --make-bed --out ${PCADIR}/UZB_unique \
  --threads $T

# 7c: LD pruning (aggressive for PCA)
log "  7c: LD pruning (1000kb, step 1, r²<0.05)"
plink2 --bfile ${PCADIR}/UZB_unique \
  --indep-pairwise 1000kb 1 0.05 \
  --out ${PCADIR}/UZB_pruned --threads $T

NPRUNED=$(wc -l < ${PCADIR}/UZB_pruned.prune.in)
log "  Independent SNPs: ${NPRUNED}"

# 7d: Create LD-pruned BED for downstream use
plink2 --bfile ${PCADIR}/UZB_unique \
  --extract ${PCADIR}/UZB_pruned.prune.in \
  --make-bed --out ${PCADIR}/UZB_ldpruned \
  --threads $T

# 7e: PCA
log "  7e: PCA (10 components)"
plink2 --bfile ${PCADIR}/UZB_ldpruned \
  --pca 10 --out ${PCADIR}/UZB_final_pca --threads $T

log "  Eigenvalues:"
cat ${PCADIR}/UZB_final_pca.eigenval
log "STEP 7 done in $(($(date +%s)-t0))s"

# ================================================================
# STEP 8: Global PCA with 1000 Genomes
# ================================================================
log "========== STEP 8: Global PCA =========="
t0=$(date +%s)

# 8a: Find common variants (by ID) between UZB and KG_reference
log "  8a: Finding common SNPs"
awk '{print $2}' ${PCADIR}/UZB_unique.bim | sort > ${GLOBALDIR}/uzb_snps.txt
awk '{print $2}' ${KG_REF}.bim | sort > ${GLOBALDIR}/kg_snps.txt
comm -12 ${GLOBALDIR}/uzb_snps.txt ${GLOBALDIR}/kg_snps.txt > ${GLOBALDIR}/common_snps.txt
N_COMMON=$(wc -l < ${GLOBALDIR}/common_snps.txt)
log "  Common SNPs: ${N_COMMON}"

# Sanity check: expect at least 10,000 common SNPs for meaningful PCA
if [ "$N_COMMON" -lt 10000 ]; then
    log "FATAL: Only ${N_COMMON} common SNPs found — likely a variant ID format mismatch!"
    log "  UZB IDs (first 3): $(head -3 ${GLOBALDIR}/uzb_snps.txt | tr '\n' ' ')"
    log "  KG  IDs (first 3): $(head -3 ${GLOBALDIR}/kg_snps.txt | tr '\n' ' ')"
    exit 1
fi

# 8b: Extract common SNPs from both datasets
plink2 --bfile ${PCADIR}/UZB_unique \
  --extract ${GLOBALDIR}/common_snps.txt \
  --make-bed --out ${GLOBALDIR}/UZB_common \
  --threads $T --memory $M

plink2 --bfile ${KG_REF} \
  --extract ${GLOBALDIR}/common_snps.txt \
  --make-bed --out ${GLOBALDIR}/KG_common \
  --threads $T --memory $M

# 8c: Merge UZB + 1000G
log "  8c: Merging datasets"
plink2 --bfile ${GLOBALDIR}/UZB_common \
  --pmerge ${GLOBALDIR}/KG_common \
  --make-bed --out ${GLOBALDIR}/UZB_1kG_merged \
  --threads $T --memory $M || {
  log "  plink2 pmerge failed, trying plink1.9 with strand handling"
  plink --bfile ${GLOBALDIR}/UZB_common \
    --bmerge ${GLOBALDIR}/KG_common.bed ${GLOBALDIR}/KG_common.bim ${GLOBALDIR}/KG_common.fam \
    --make-bed --out ${GLOBALDIR}/UZB_1kG_merged --allow-no-sex || {
    log "  Excluding $(wc -l < ${GLOBALDIR}/UZB_1kG_merged-merge.missnp) mismatched SNPs"
    plink --bfile ${GLOBALDIR}/UZB_common \
      --exclude ${GLOBALDIR}/UZB_1kG_merged-merge.missnp \
      --make-bed --out ${GLOBALDIR}/UZB_common_clean --allow-no-sex
    plink --bfile ${GLOBALDIR}/KG_common \
      --exclude ${GLOBALDIR}/UZB_1kG_merged-merge.missnp \
      --make-bed --out ${GLOBALDIR}/KG_common_clean --allow-no-sex
    plink --bfile ${GLOBALDIR}/UZB_common_clean \
      --bmerge ${GLOBALDIR}/KG_common_clean.bed ${GLOBALDIR}/KG_common_clean.bim ${GLOBALDIR}/KG_common_clean.fam \
      --make-bed --out ${GLOBALDIR}/UZB_1kG_merged --allow-no-sex
  }
}

N_MERGED_S=$(wc -l < ${GLOBALDIR}/UZB_1kG_merged.fam)
N_MERGED_V=$(wc -l < ${GLOBALDIR}/UZB_1kG_merged.bim)
log "  Merged: ${N_MERGED_S} samples × ${N_MERGED_V} variants"

# 8d: Global PCA
log "  8d: Global PCA (10 components)"
plink2 --bfile ${GLOBALDIR}/UZB_1kG_merged \
  --pca 10 --out ${GLOBALDIR}/GLOBAL_PCA --threads $T

log "  Global eigenvalues:"
cat ${GLOBALDIR}/GLOBAL_PCA.eigenval

# 8e: Population mapping (needed for Fst, PBS, ADMIXTURE)
log "  8e: Building population mapping"
# Strategy: Use a known-good 3-column clusters file, or build from 1000G panel file.
# CRITICAL: ~/v2/pca/pop_mapping.txt is BROKEN (single concatenated column).
# Validate any candidate has ≥3 tab/space-separated fields.

KG_POP_MAP=""
# Priority 1: v2 clusters file (known correct 3-column: FID IID POP)
# Priority 2: 1000G panel file (has sample, pop, super_pop, gender)
for p in ~/v2/pbs_v2/clusters.txt \
         /staging/ALSU-analysis/Fst_analysis/1000G_data/integrated_call_samples_v3.20130502.ALL.panel \
         /staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/1000g_grch38/1000g_panel.txt; do
    if [ -f "$p" ]; then
        # Validate: must have ≥3 fields on data lines
        NFIELDS=$(awk 'NR>1 && NF>=3 {print NF; exit}' "$p")
        if [ -n "$NFIELDS" ] && [ "$NFIELDS" -ge 3 ]; then
            KG_POP_MAP="$p"
            log "  Pop mapping source: ${KG_POP_MAP} (${NFIELDS} fields)"
            break
        else
            log "  SKIP (bad format): $p"
        fi
    fi
done

if [ -n "$KG_POP_MAP" ]; then
    # Detect format: clusters.txt has FID IID POP; panel has sample pop super_pop gender
    IS_PANEL=false
    if head -1 "$KG_POP_MAP" | grep -qi "super_pop"; then
        IS_PANEL=true
        log "  Detected 1000G panel format"
    fi

    if [ "$IS_PANEL" = true ]; then
        # Panel file: col1=sample, col3=super_pop
        awk 'NR==FNR && FNR>1 {
            gsub(/\r/,"")
            pop[$1] = $3
            next
        }
        {
            iid = $2
            if (iid in pop && pop[iid] ~ /^(EUR|EAS|SAS|AFR|AMR)$/) {
                print $1, iid, pop[iid]
            } else {
                print $1, iid, "UZB"
            }
        }' "$KG_POP_MAP" ${GLOBALDIR}/UZB_1kG_merged.fam > ${GLOBALDIR}/clusters.txt
    else
        # 3-column clusters file: col2=IID, col3=POP
        awk 'NR==FNR {
            gsub(/\r/,"")
            pop[$2] = $3
            next
        }
        {
            iid = $2
            if (iid in pop && pop[iid] ~ /^(EUR|EAS|SAS|AFR|AMR|UZB)$/) {
                print $1, iid, pop[iid]
            } else {
                print $1, iid, "UZB"
            }
        }' "$KG_POP_MAP" ${GLOBALDIR}/UZB_1kG_merged.fam > ${GLOBALDIR}/clusters.txt
    fi
else
    log "  WARNING: No pop mapping found, inferring from sample naming"
    awk '{
        iid = $2
        if (iid ~ /^(HG|NA)/) print $1, iid, "UNKNOWN_1KG"
        else print $1, iid, "UZB"
    }' ${GLOBALDIR}/UZB_1kG_merged.fam > ${GLOBALDIR}/clusters.txt
fi

log "  Population counts:"
awk '{print $3}' ${GLOBALDIR}/clusters.txt | sort | uniq -c | sort -rn | while read n p; do
    log "    $p: $n"
done

# Validate pop mapping: 1000G should have ~2548 samples across 5 superpops
N_1KG_MAPPED=$(awk '$3 ~ /^(EUR|EAS|SAS|AFR|AMR)$/' ${GLOBALDIR}/clusters.txt | wc -l)
N_UZB_MAPPED=$(awk '$3=="UZB"' ${GLOBALDIR}/clusters.txt | wc -l)
if [ "$N_1KG_MAPPED" -lt 2000 ]; then
    log "FATAL: Only ${N_1KG_MAPPED} 1000G samples mapped to superpops — pop mapping is broken!"
    log "  First 5 lines of clusters.txt:"
    head -5 ${GLOBALDIR}/clusters.txt | while read line; do log "    $line"; done
    exit 1
fi
log "  Validated: ${N_1KG_MAPPED} 1000G + ${N_UZB_MAPPED} UZB samples"

cp ${GLOBALDIR}/clusters.txt ${FSTDIR}/clusters.txt
cp ${GLOBALDIR}/clusters.txt ${PBSDIR}/clusters.txt
log "STEP 8 done in $(($(date +%s)-t0))s"

# ================================================================
# STEP 9: Genome-wide Fst Analysis
# ================================================================
log "========== STEP 9: Fst Analysis =========="
t0=$(date +%s)
cd "$FSTDIR"

MERGED=${GLOBALDIR}/UZB_1kG_merged

# 9a: Create per-population sample lists
log "  9a: Per-population sample lists"
for POP in UZB EUR EAS SAS AFR; do
    awk -v p="$POP" '$3==p {print $1, $2}' clusters.txt > keep_${POP}.txt
    log "    $POP: $(wc -l < keep_${POP}.txt) samples"
done

# 9b: Population frequencies
log "  9b: Allele frequencies"
for POP in UZB EUR EAS SAS AFR; do
    plink --bfile "$MERGED" \
          --keep keep_${POP}.txt \
          --freq \
          --out freq_${POP} \
          --allow-no-sex --silent
done
log "  Frequencies computed"

# 9c: Pairwise Fst (5 pairs)
log "  9c: Pairwise Fst"
compute_fst() {
    local P1=$1 P2=$2
    cat keep_${P1}.txt keep_${P2}.txt > keep_${P1}_${P2}.txt
    awk -v p="$P1" '{print $1, $2, p}' keep_${P1}.txt > within_${P1}_${P2}.txt
    awk -v p="$P2" '{print $1, $2, p}' keep_${P2}.txt >> within_${P1}_${P2}.txt
    plink --bfile "$MERGED" \
          --keep keep_${P1}_${P2}.txt \
          --fst --within within_${P1}_${P2}.txt \
          --out fst_${P1}_${P2} \
          --allow-no-sex --silent
    local wfst=$(grep 'Weighted' fst_${P1}_${P2}.log | awk '{print $NF}')
    log "    ${P1} vs ${P2}: weighted Fst = ${wfst}"
}

compute_fst UZB EUR
compute_fst UZB EAS
compute_fst EUR EAS
compute_fst UZB SAS
compute_fst UZB AFR

log "STEP 9 done in $(($(date +%s)-t0))s"

# ================================================================
# STEP 10: PBS Analysis
# ================================================================
log "========== STEP 10: PBS Analysis =========="
t0=$(date +%s)
cd "$PBSDIR"

# Copy Fst and freq data from step 9
cp ${FSTDIR}/fst_*.fst . 2>/dev/null || true
cp ${FSTDIR}/freq_*.frq . 2>/dev/null || true

# PBS calculation via embedded Python
log "  Computing PBS..."
python3 << 'PYEOF'
import math, json, os, sys

def read_frq(path):
    freqs = {}
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.split()
            if len(p)>=5:
                freqs[p[1]] = float(p[4]) if p[4]!='NA' else None
    return freqs

def read_fst(path):
    fst = {}
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.split()
            if len(p)>=5:
                fst[p[1]] = float(p[4]) if p[4] not in ('nan','NA') else 0.0
    return fst

def fst_to_T(f):
    return -math.log(1.0 - max(0.0, min(f, 0.999)))

# Read data
freq = {p: read_frq(f"freq_{p}.frq") for p in ['UZB','EUR','EAS','SAS','AFR']}
fst_data = {
    'UZB_EUR': read_fst("fst_UZB_EUR.fst"),
    'UZB_EAS': read_fst("fst_UZB_EAS.fst"),
    'EUR_EAS': read_fst("fst_EUR_EAS.fst"),
    'UZB_SAS': read_fst("fst_UZB_SAS.fst"),
    'UZB_AFR': read_fst("fst_UZB_AFR.fst"),
}

# BIM for chromosome/position
merged_bim = os.environ.get('MERGED_BIM',
    '/staging/ALSU-analysis/spring2026/global_pca/UZB_1kG_merged.bim')
bim = {}
with open(merged_bim) as f:
    for line in f:
        p = line.split()
        bim[p[1]] = (p[0], int(p[3]))

# Intersection of all datasets
common = set(fst_data['UZB_EUR']) & set(fst_data['UZB_EAS']) & set(fst_data['EUR_EAS'])
for p in ['UZB','EUR','EAS','SAS','AFR']:
    common &= set(freq[p])
print(f"  Common SNPs for PBS: {len(common)}")

results = []
for snp in common:
    T_ue = fst_to_T(fst_data['UZB_EUR'][snp])
    T_ua = fst_to_T(fst_data['UZB_EAS'][snp])
    T_ea = fst_to_T(fst_data['EUR_EAS'][snp])

    pbs_uzb = (T_ue + T_ua - T_ea) / 2.0
    pbs_eur = (T_ue + T_ea - T_ua) / 2.0
    pbs_eas = (T_ua + T_ea - T_ue) / 2.0

    mafs = {p: freq[p].get(snp) or 0 for p in ['UZB','EUR','EAS','SAS','AFR']}
    delta_af = max(abs(mafs['UZB']-mafs[p]) for p in ['EUR','EAS','SAS','AFR'])
    near_priv = mafs['UZB'] >= 0.05 and all(mafs[p] <= 0.01 for p in ['EUR','EAS','SAS','AFR'])

    chrom, bp = bim.get(snp, ('?', 0))
    results.append({
        'snp': snp, 'chr': chrom, 'bp': bp,
        'pbs_uzb': round(pbs_uzb, 6), 'pbs_eur': round(pbs_eur, 6), 'pbs_eas': round(pbs_eas, 6),
        'maf_uzb': round(mafs['UZB'], 6), 'maf_eur': round(mafs['EUR'], 6),
        'maf_eas': round(mafs['EAS'], 6), 'maf_sas': round(mafs['SAS'], 6),
        'maf_afr': round(mafs['AFR'], 6),
        'delta_af': round(delta_af, 4), 'near_private': near_priv,
        'tier1': pbs_uzb >= 0.3, 'tier2': delta_af >= 0.3, 'tier3': near_priv,
    })

# Sort by PBS descending
results.sort(key=lambda x: -x['pbs_uzb'])

# Summary stats
pbs_vals = [r['pbs_uzb'] for r in results]
import statistics
print(f"  PBS stats: mean={statistics.mean(pbs_vals):.5f}, median={statistics.median(pbs_vals):.5f}, stdev={statistics.stdev(pbs_vals):.5f}")
print(f"  PBS>=0.3 (Tier 1): {sum(1 for r in results if r['tier1'])}")
print(f"  PBS>=0.1: {sum(1 for r in results if r['pbs_uzb']>=0.1)}")
print(f"  delta_AF>=0.3 (Tier 2): {sum(1 for r in results if r['tier2'])}")
print(f"  Near-private (Tier 3): {sum(1 for r in results if r['tier3'])}")

# Save full results
with open("pbs_all.tsv", 'w') as f:
    cols = ['snp','chr','bp','pbs_uzb','pbs_eur','pbs_eas','maf_uzb','maf_eur','maf_eas','maf_sas','maf_afr','delta_af','near_private','tier1','tier2','tier3']
    f.write('\t'.join(cols) + '\n')
    for r in results:
        f.write('\t'.join(str(r[c]) for c in cols) + '\n')

# Candidates JSON
candidates = [r for r in results if r['tier1'] or r['tier2'] or r['tier3']]
with open("pbs_candidates.json", 'w') as f:
    json.dump(candidates, f, indent=2)
print(f"  Candidates (any tier): {len(candidates)}")

# Stats JSON
stats = {
    'total_snps': len(results),
    'mean_pbs': round(statistics.mean(pbs_vals), 5),
    'median_pbs': round(statistics.median(pbs_vals), 5),
    'stdev_pbs': round(statistics.stdev(pbs_vals), 5),
    'tier1_count': sum(1 for r in results if r['tier1']),
    'tier2_count': sum(1 for r in results if r['tier2']),
    'tier3_count': sum(1 for r in results if r['tier3']),
}
with open("pbs_stats.json", 'w') as f:
    json.dump(stats, f, indent=2)

print("  Top 10 PBS SNPs:")
for r in results[:10]:
    print(f"    {r['snp']}: PBS={r['pbs_uzb']:.4f} MAF_UZB={r['maf_uzb']:.4f} chr{r['chr']}:{r['bp']}")
PYEOF

log "STEP 10 done in $(($(date +%s)-t0))s"

# ================================================================
# STEP 11: ADMIXTURE
# ================================================================
log "========== STEP 11: ADMIXTURE =========="
t0=$(date +%s)
cd "$ADMIXDIR"

# 11a: Prepare LD-pruned dataset for ADMIXTURE
log "  11a: QC + LD prune for ADMIXTURE"
plink2 --bfile ${GLOBALDIR}/UZB_1kG_merged \
  --maf 0.01 --geno 0.02 --hwe 1e-6 \
  --make-bed --out ${ADMIXDIR}/global_qc \
  --threads $T --memory $M

plink2 --bfile ${ADMIXDIR}/global_qc \
  --indep-pairwise 50 10 0.1 \
  --out ${ADMIXDIR}/global_ld --threads $T

plink2 --bfile ${ADMIXDIR}/global_qc \
  --extract ${ADMIXDIR}/global_ld.prune.in \
  --make-bed --out ${ADMIXDIR}/global_for_admixture \
  --threads $T

N_ADM_S=$(wc -l < ${ADMIXDIR}/global_for_admixture.fam)
N_ADM_V=$(wc -l < ${ADMIXDIR}/global_for_admixture.bim)
log "  ADMIXTURE input: ${N_ADM_S} samples × ${N_ADM_V} variants"

# 11b: Run ADMIXTURE K=2..8 with 5-fold CV (32 threads)
log "  11b: ADMIXTURE K=2..8"
for K in $(seq 2 8); do
    t1=$(date +%s)
    log "  Starting K=${K}..."
    admixture --cv=5 -j32 ${ADMIXDIR}/global_for_admixture.bed $K \
        2>&1 | tee ${ADMIXDIR}/admixture_K${K}.log
    cv=$(grep -i 'CV error' ${ADMIXDIR}/admixture_K${K}.log | tail -1)
    log "  K=${K} done in $(($(date +%s)-t1))s — ${cv}"
done

# 11c: CV Error Summary
log ""
log "=== CV ERROR SUMMARY ==="
for K in $(seq 2 8); do
    grep -i 'CV error' ${ADMIXDIR}/admixture_K${K}.log 2>/dev/null || echo "K=${K}: (not found)"
done

# 11d: Also run UZB-only ADMIXTURE K=2..8
log ""
log "  11d: UZB-only ADMIXTURE"
# LD-prune UZB-only data
plink2 --bfile ${PCADIR}/UZB_unique \
  --maf 0.01 --geno 0.02 --hwe 1e-6 \
  --make-bed --out ${ADMIXDIR}/uzb_qc --threads $T

plink2 --bfile ${ADMIXDIR}/uzb_qc \
  --indep-pairwise 50 10 0.1 \
  --out ${ADMIXDIR}/uzb_ld --threads $T

plink2 --bfile ${ADMIXDIR}/uzb_qc \
  --extract ${ADMIXDIR}/uzb_ld.prune.in \
  --make-bed --out ${ADMIXDIR}/uzb_for_admixture --threads $T

log "  UZB ADMIXTURE input: $(wc -l < ${ADMIXDIR}/uzb_for_admixture.fam) samples × $(wc -l < ${ADMIXDIR}/uzb_for_admixture.bim) variants"

for K in $(seq 2 8); do
    t1=$(date +%s)
    log "  UZB K=${K}..."
    admixture --cv=5 -j32 ${ADMIXDIR}/uzb_for_admixture.bed $K \
        2>&1 | tee ${ADMIXDIR}/uzb_admixture_K${K}.log
    cv=$(grep -i 'CV error' ${ADMIXDIR}/uzb_admixture_K${K}.log | tail -1)
    log "  UZB K=${K} done in $(($(date +%s)-t1))s — ${cv}"
done

log ""
log "=== UZB CV ERROR SUMMARY ==="
for K in $(seq 2 8); do
    grep -i 'CV error' ${ADMIXDIR}/uzb_admixture_K${K}.log 2>/dev/null || echo "K=${K}: (not found)"
done

log "STEP 11 done in $(($(date +%s)-t0))s"

# ================================================================
# FINAL SUMMARY
# ================================================================
log ""
log "========================================="
log "PIPELINE COMPLETE: Steps 7–11"
log "========================================="
log "Step 7 (Local PCA): $(wc -l < ${PCADIR}/UZB_final_pca.eigenvec) samples, ${NPRUNED} LD-pruned SNPs"
log "Step 8 (Global PCA): ${N_MERGED_S} samples × ${N_MERGED_V} variants"
log "Step 9 (Fst): 5 pairwise Fst computed"
log "Step 10 (PBS): $(wc -l < ${PBSDIR}/pbs_all.tsv) rows"
log "Step 11 (ADMIXTURE): K=2..8 global + UZB-only"
log ""
log "Output locations:"
log "  ${PCADIR}/"
log "  ${GLOBALDIR}/"
log "  ${FSTDIR}/"
log "  ${PBSDIR}/"
log "  ${ADMIXDIR}/"
log ""
log "Disk usage:"
du -sh ${PCADIR} ${GLOBALDIR} ${FSTDIR} ${PBSDIR} ${ADMIXDIR}
df -h /staging/
log "=== ALL DONE ==="
