#!/bin/bash
# =============================================================================
# ALSU RPL GWAS — DRAGEN execution script
# Run on DRAGEN server: /staging/ALSU-analysis/spring2026/gwas/
#
# Prerequisites:
#   - PLINK2 binary at /staging/conda/envs/bioinfo/bin/plink2
#   - UZB_imputed_HQ_clean.{bed,bim,fam} in spring2026/post_imputation/
#   - gwas_pheno.txt and gwas_covar.txt uploaded to spring2026/gwas/
#
# Usage:
#   nohup bash /staging/ALSU-analysis/spring2026/gwas/run_gwas.sh > gwas.log 2>&1 &
# =============================================================================

set -euo pipefail
export PATH=/usr/local/bin:/usr/bin:/bin:/staging/conda/envs/bioinfo/bin

WORKDIR=/staging/ALSU-analysis/spring2026/gwas
BFILE=/staging/ALSU-analysis/spring2026/post_imputation/UZB_imputed_HQ_clean
PHENO=${WORKDIR}/gwas_pheno.txt
COVAR=${WORKDIR}/gwas_covar.txt
THREADS=16
MEMORY=200000

cd ${WORKDIR}
echo "$(date): GWAS pipeline started"
echo "BFILE: ${BFILE}"
echo "PHENO: ${PHENO}"
echo "COVAR: ${COVAR}"

# Verify inputs exist
for f in ${BFILE}.bed ${BFILE}.bim ${BFILE}.fam ${PHENO} ${COVAR}; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Missing file: $f"
        exit 1
    fi
done

NSNP=$(wc -l < ${BFILE}.bim)
NSAMP=$(wc -l < ${BFILE}.fam)
echo "Input: ${NSAMP} samples, ${NSNP} SNPs"

# =============================================================================
# TIER 1: Targeted test — 8 PBS candidates
# =============================================================================
echo ""
echo "$(date): ========== TIER 1: PBS CANDIDATE TESTING =========="

# Create PBS candidate position list (CHR:POS from Step 10/12)
# We'll find matching SNP IDs from the BIM file by position
cat > ${WORKDIR}/pbs_positions.txt << 'EOF'
12	22967890
12	125520190
12	5664803
11	20665570
5	53879140
3	133749168
11	207698
10	8512594
EOF

# Extract actual SNP IDs from BIM by matching chromosomal position
echo "Looking up PBS candidate SNP IDs from BIM..."
> ${WORKDIR}/pbs_candidates.txt
while IFS=$'\t' read -r CHR POS; do
    awk -v chr="$CHR" -v pos="$POS" '$1==chr && $4==pos {print $2}' ${BFILE}.bim >> ${WORKDIR}/pbs_candidates.txt
done < ${WORKDIR}/pbs_positions.txt
N_FOUND=$(wc -l < ${WORKDIR}/pbs_candidates.txt)
echo "  Found ${N_FOUND}/8 PBS candidates in BIM"
cat ${WORKDIR}/pbs_candidates.txt

if [ "${N_FOUND}" -gt 0 ]; then
    echo "Testing ${N_FOUND} PBS candidates..."
    plink2 --bfile ${BFILE} \
        --pheno ${PHENO} \
        --pheno-name RPL \
        --1 \
        --covar ${COVAR} \
        --covar-name AGE,BMI,PC1,PC2,PC3 \
        --covar-variance-standardize \
        --extract ${WORKDIR}/pbs_candidates.txt \
        --glm hide-covar cols=chrom,pos,ref,alt,a1freq,test,nobs,orbeta,se,ci,tz,p \
        --ci 0.95 \
        --out ${WORKDIR}/tier1_pbs \
        --threads ${THREADS} \
        --memory ${MEMORY}
else
    echo "WARNING: No PBS candidates found in BIM — skipping Tier 1"
fi

echo "$(date): Tier 1 complete"
if [ -f ${WORKDIR}/tier1_pbs.RPL.glm.logistic.hybrid ]; then
    echo "--- PBS CANDIDATE RESULTS ---"
    cat ${WORKDIR}/tier1_pbs.RPL.glm.logistic.hybrid
    echo ""
    echo "Bonferroni threshold for ${N_FOUND} tests: p < $(python3 -c "print(f'{0.05/max(1,${N_FOUND}):.6f}')")"
fi

# =============================================================================
# TIER 2: Known RPL variant screen
# =============================================================================
echo ""
echo "$(date): ========== TIER 2: KNOWN RPL VARIANT SCREEN =========="

# Known RPL-associated variants
cat > ${WORKDIR}/rpl_known_variants.txt << 'EOF'
rs6025
rs1799963
rs1801133
rs1799889
EOF

# Extract allele frequencies for known variants
plink2 --bfile ${BFILE} \
    --freq \
    --pheno ${PHENO} \
    --pheno-name RPL \
    --1 \
    --loop-cats RPL \
    --out ${WORKDIR}/rpl_known_freq \
    --threads ${THREADS} \
    --memory ${MEMORY} || echo "WARN: --loop-cats may not work; running overall freq"

plink2 --bfile ${BFILE} \
    --freq \
    --out ${WORKDIR}/allele_freq \
    --threads ${THREADS}

# Test known variants if they exist in the data
echo "Checking known RPL variants in dataset..."
for RS in rs6025 rs1799963 rs1801133 rs1799889; do
    FOUND=$(grep -c "${RS}" ${BFILE}.bim || true)
    echo "  ${RS}: ${FOUND} found in .bim"
done

# If rsIDs present, test them
plink2 --bfile ${BFILE} \
    --pheno ${PHENO} \
    --pheno-name RPL \
    --1 \
    --covar ${COVAR} \
    --covar-name AGE,BMI,PC1,PC2,PC3 \
    --covar-variance-standardize \
    --extract ${WORKDIR}/rpl_known_variants.txt \
    --glm hide-covar cols=chrom,pos,ref,alt,a1freq,test,nobs,orbeta,se,ci,tz,p \
    --ci 0.95 \
    --out ${WORKDIR}/tier2_rpl_known \
    --threads ${THREADS} \
    --memory ${MEMORY} || echo "WARN: Some known variants may not be in dataset"

echo "$(date): Tier 2 complete"

# =============================================================================
# TIER 3: Genome-wide scan
# =============================================================================
echo ""
echo "$(date): ========== TIER 3: GENOME-WIDE ASSOCIATION =========="

plink2 --bfile ${BFILE} \
    --pheno ${PHENO} \
    --pheno-name RPL \
    --1 \
    --covar ${COVAR} \
    --covar-name AGE,BMI,PC1,PC2,PC3 \
    --covar-variance-standardize \
    --glm hide-covar cols=chrom,pos,ref,alt,a1freq,test,nobs,orbeta,se,ci,tz,p \
    --ci 0.95 \
    --out ${WORKDIR}/tier3_gwas \
    --threads ${THREADS} \
    --memory ${MEMORY}

echo "$(date): Tier 3 complete"

# Count significant results
if [ -f ${WORKDIR}/tier3_gwas.RPL.glm.logistic.hybrid ]; then
    TOTAL=$(tail -n +2 ${WORKDIR}/tier3_gwas.RPL.glm.logistic.hybrid | wc -l)
    GW_SIG=$(awk 'NR>1 && $NF+0 < 5e-8' ${WORKDIR}/tier3_gwas.RPL.glm.logistic.hybrid | wc -l)
    SUGGESTIVE=$(awk 'NR>1 && $NF+0 < 1e-5' ${WORKDIR}/tier3_gwas.RPL.glm.logistic.hybrid | wc -l)
    echo "--- GENOME-WIDE RESULTS ---"
    echo "Total SNPs tested: ${TOTAL}"
    echo "Genome-wide significant (p < 5e-8): ${GW_SIG}"
    echo "Suggestive (p < 1e-5): ${SUGGESTIVE}"
    echo ""
    echo "Top 20 hits:"
    head -1 ${WORKDIR}/tier3_gwas.RPL.glm.logistic.hybrid
    sort -k$(head -1 ${WORKDIR}/tier3_gwas.RPL.glm.logistic.hybrid | awk '{for(i=1;i<=NF;i++) if($i=="P") print i}') -g \
        ${WORKDIR}/tier3_gwas.RPL.glm.logistic.hybrid | head -20
fi

# =============================================================================
# QC: Genomic inflation (lambda)
# =============================================================================
echo ""
echo "$(date): ========== QC: GENOMIC INFLATION =========="

python3 -c "
import math
results = []
with open('${WORKDIR}/tier3_gwas.RPL.glm.logistic.hybrid') as f:
    header = f.readline().strip().split('\t')
    p_idx = header.index('P')
    for line in f:
        parts = line.strip().split('\t')
        try:
            p = float(parts[p_idx])
            if p > 0 and p < 1:
                # Convert to chi-squared
                from scipy.stats import norm, chi2
                z = norm.ppf(1 - p/2)
                results.append(z**2)
        except (ValueError, IndexError):
            pass

if results:
    results.sort()
    median_chi2 = results[len(results)//2]
    lambda_gc = median_chi2 / 0.4549364
    print(f'Lambda_GC = {lambda_gc:.4f}')
    print(f'  (Ideally ~1.0; values > 1.05 suggest residual population stratification)')
    if lambda_gc > 1.10:
        print('  WARNING: Lambda > 1.10 — consider adding more PCs or using mixed models')
    elif lambda_gc < 0.95:
        print('  NOTE: Lambda < 0.95 — test may be conservative (small sample size effect)')
    else:
        print('  OK: Lambda within acceptable range')
" 2>/dev/null || echo "WARN: scipy not available; compute lambda manually"

echo ""
echo "$(date): ========== ALL GWAS ANALYSES COMPLETE =========="
echo ""
echo "Output files:"
echo "  Tier 1 (PBS):     ${WORKDIR}/tier1_pbs.RPL.glm.logistic.hybrid"
echo "  Tier 2 (Known):   ${WORKDIR}/tier2_rpl_known.RPL.glm.logistic.hybrid"
echo "  Tier 3 (GW):      ${WORKDIR}/tier3_gwas.RPL.glm.logistic.hybrid"
echo "  Frequencies:      ${WORKDIR}/allele_freq.afreq"
