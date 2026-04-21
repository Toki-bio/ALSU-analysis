#!/bin/bash
set -uo pipefail
export PATH=/usr/local/bin:/usr/bin:/bin:/staging/conda/envs/bioinfo/bin

WORKDIR=/staging/ALSU-analysis/spring2026/gwas
SENSDIR=${WORKDIR}/sensitivity
BFILE=/staging/ALSU-analysis/spring2026/post_imputation/UZB_imputed_HQ_clean
COVAR=${WORKDIR}/gwas_covar.txt

mkdir -p ${SENSDIR}/results
cd ${SENSDIR}

run_one () {
    local NAME=$1
    local PHENO=${SENSDIR}/gwas_pheno_${NAME}.txt
    local OUT=${SENSDIR}/results/${NAME}
    echo ""
    echo "$(date): ===== ${NAME} ====="
    plink2 --bfile ${BFILE} \
        --pheno ${PHENO} --pheno-name RPL --1 \
        --covar ${COVAR} --covar-name AGE,BMI,PC1,PC2,PC3 \
        --covar-variance-standardize \
        --glm hide-covar cols=chrom,pos,ref,alt,a1freq,test,nobs,orbeta,se,ci,tz,p \
        --ci 0.95 --out ${OUT} --threads 16 --memory 200000 \
        2>&1 | tail -8
    local RES=${OUT}.RPL.glm.logistic.hybrid
    if [ ! -f ${RES} ]; then
        echo "  ERROR: result file missing"
        return 1
    fi
    # Single-pass awk: counts + median p + top 10 (without sorting full file)
    awk -F'\t' -v OFS='\t' -v OUT="${OUT}" '
    NR==1 { for(i=1;i<=NF;i++) if($i=="P") pi=i; print > (OUT "_top10.tsv"); next }
    $pi == "NA" || $pi+0 <= 0 || $pi+0 >= 1 { next }
    {
        n++
        if($pi+0 < 5e-8) gw++
        if($pi+0 < 1e-5) sug++
        # reservoir of top 10 lowest p
        if(k<10) { rows[++k]=$0; pv[k]=$pi+0; if(pv[k]>maxp){maxp=pv[k]; maxi=k} }
        else if($pi+0 < maxp) {
            rows[maxi]=$0; pv[maxi]=$pi+0
            maxp=pv[1]; maxi=1
            for(j=2;j<=10;j++) if(pv[j]>maxp){maxp=pv[j]; maxi=j}
        }
        # also store all p for median (memory ok ~10M floats)
        ap[n]=$pi+0
    }
    END {
        printf "  tested=%d  GW(<5e-8)=%d  suggestive(<1e-5)=%d\n", n, gw+0, sug+0
        # median
        asort(ap)
        mp = ap[int(n/2)]
        printf "  median_p=%.6g\n", mp
        # write top-10 sorted
        for(j=1;j<=k;j++) for(m=j+1;m<=k;m++) if(pv[m]<pv[j]) {
            t=pv[m]; pv[m]=pv[j]; pv[j]=t
            t=rows[m]; rows[m]=rows[j]; rows[j]=t
        }
        for(j=1;j<=k;j++) print rows[j] >> (OUT "_top10.tsv")
    }
    ' ${RES}
    echo "  -- top 10 --"
    cat ${OUT}_top10.tsv
}

# S1 already done; skip
for n in S2_review_as_case S3_gray_as_case S4_probable; do
    run_one $n
done

echo ""
echo "$(date): ALL DONE"
