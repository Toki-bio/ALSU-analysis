#!/bin/bash
# Independent check of reference-pair FST: all common biallelic SNVs on chr22 (NYGC 30x), PLINK weighted/mean FST.
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH
export TMPDIR=/staging/tmp/scratch
OUT=/staging/tmp/scratch/fst22_check
mkdir -p "$OUT"
cd "$OUT" || exit 1
F=/staging/ALSU-analysis/spring2026/full_expanded_cohort
K=$F/pbs_refqc_2026-10/nygc30x_1256
echo "start $(date)"
bcftools view -m2 -M2 -v snps $F/nygc30x/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -Oz -o c22.vcf.gz 2> /dev/null
plink --vcf c22.vcf.gz --double-id --maf 0.01 --make-bed --out c22 --threads 8 > c22.make.log 2>&1
grep -i "variants and\|pass filters" c22.make.log || true
for P in EUR EAS SAS AFR; do awk -v p=$P '{print $1,$2,p}' $K/keep_$P.txt > w_$P.txt; done
for PAIR in EUR_EAS EUR_SAS EAS_SAS AFR_EUR AFR_EAS AFR_SAS; do
  A=${PAIR%_*}; B=${PAIR#*_}
  cat w_$A.txt w_$B.txt > within_$PAIR.txt
  cut -d' ' -f1,2 within_$PAIR.txt > keep_$PAIR.txt
  plink --bfile c22 --keep keep_$PAIR.txt --within within_$PAIR.txt --fst --out fst_$PAIR --threads 8 > /dev/null 2>&1
  echo "$PAIR $(grep -E 'Mean Fst|Weighted Fst' fst_$PAIR.log | tr '\n' ' ') n_snps=$(awk 'NR>1 && $5!="nan"' fst_$PAIR.fst | wc -l)"
done
echo "finished $(date)"
