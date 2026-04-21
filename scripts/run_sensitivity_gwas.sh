#!/bin/bash
# Run 4 sensitivity GWAS analyses against UZB_imputed_HQ_clean.
# Each uses the existing covariate file gwas_covar.txt with PCs+AGE+BMI.
set -euo pipefail
export PATH=/usr/local/bin:/usr/bin:/bin:/staging/conda/envs/bioinfo/bin

WORKDIR=/staging/ALSU-analysis/spring2026/gwas
SENSDIR=${WORKDIR}/sensitivity
BFILE=/staging/ALSU-analysis/spring2026/post_imputation/UZB_imputed_HQ_clean
COVAR=${WORKDIR}/gwas_covar.txt
THREADS=16
MEMORY=200000

mkdir -p ${SENSDIR}/results
cd ${SENSDIR}

run_one () {
    local NAME=$1
    local PHENO=${SENSDIR}/gwas_pheno_${NAME}.txt
    local OUT=${SENSDIR}/results/${NAME}
    echo ""
    echo "$(date): ===== ${NAME} ====="
    echo "Pheno: ${PHENO}"
    plink2 --bfile ${BFILE} \
        --pheno ${PHENO} --pheno-name RPL --1 \
        --covar ${COVAR} \
        --covar-name AGE,BMI,PC1,PC2,PC3 \
        --covar-variance-standardize \
        --glm hide-covar cols=chrom,pos,ref,alt,a1freq,test,nobs,orbeta,se,ci,tz,p \
        --ci 0.95 \
        --out ${OUT} \
        --threads ${THREADS} --memory ${MEMORY}

    local RES=${OUT}.RPL.glm.logistic.hybrid
    if [ -f ${RES} ]; then
        local TOTAL GW SUG
        TOTAL=$(tail -n +2 ${RES} | wc -l)
        GW=$(awk 'NR>1 && $NF+0 < 5e-8 && $NF != "NA"' ${RES} | wc -l)
        SUG=$(awk 'NR>1 && $NF+0 < 1e-5 && $NF != "NA"' ${RES} | wc -l)
        echo "  ${NAME}: tested=${TOTAL}  GW(<5e-8)=${GW}  suggestive(<1e-5)=${SUG}"

        # lambda_GC
        python3 -c "
import math
from statistics import median
ps=[]
with open('${RES}') as f:
    h=f.readline().rstrip().split('\t')
    pi=h.index('P')
    for line in f:
        p=line.rstrip().split('\t')[pi]
        try:
            v=float(p)
            if 0<v<1: ps.append(v)
        except: pass
ps.sort()
mp=ps[len(ps)//2]
# chi2 inverse with 1 df via Wilson-Hilferty approx
import scipy.stats as st
chi=st.chi2.ppf(1-mp,1)
print(f'  lambda_GC = {chi/0.4549364:.4f}  (n={len(ps)})')
" 2>/dev/null || echo "  (lambda calc skipped - scipy missing)"

        # Top 5 hits
        echo "  Top 5 hits:"
        head -1 ${RES} | awk -F'\t' '{for(i=1;i<=NF;i++) printf "  %d:%s\n",i,$i}'
        head -1 ${RES} > ${OUT}_top5.tsv
        sort -k$(head -1 ${RES} | tr '\t' '\n' | grep -n '^P$' | cut -d: -f1) -g ${RES} | head -5 >> ${OUT}_top5.tsv
        cat ${OUT}_top5.tsv
    fi
}

run_one S1_strict3
run_one S2_review_as_case
run_one S3_gray_as_case
run_one S4_probable

echo ""
echo "$(date): ALL SENSITIVITY ANALYSES COMPLETE"
ls -la ${SENSDIR}/results/
