#!/bin/bash
# Project the 1,193 relatedness-dataset samples (ALSU + GWAS2026 arrays) onto the fixed NYGC K=4 allele frequencies
# (admixture -P with adm_in.4.P), using the SNPs present in both data sets with identical allele sets.
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH
export TMPDIR=/staging/tmp/scratch
OUT=/staging/tmp/scratch/proj_k4
mkdir -p "$OUT"
cd "$OUT" || exit 1
B=/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/cross_gwas2026_alsu_winter
A=/staging/ALSU-analysis/spring2026/full_expanded_cohort/nygc30x_1256/adm
echo "start $(date)"
# overlap by chr:pos with the same unordered allele pair; keep merged.bim order
awk 'NR==FNR{a1[$1":"$4]=$5; a2[$1":"$4]=$6; row[$1":"$4]=FNR; next}
     ($1":"$4) in a1 {k=$1":"$4; x=$5; y=$6;
       if ((x==a1[k] && y==a2[k]) || (x==a2[k] && y==a1[k])) print $2"\t"k"\t"a1[k]"\t"row[k]}' $A/adm_in.bim $B/merged.bim > keep_map.tsv
wc -l keep_map.tsv
cut -f1 keep_map.tsv > keep_ids.txt
plink --bfile $B/merged --extract keep_ids.txt --a1-allele keep_map.tsv 3 1 --make-bed --out sub --threads 8 > sub.plink.log 2>&1
grep -i "variants and\|pass filters" sub.plink.log || true
# check orientation: A1/A2 in sub.bim must equal adm_in.bim for every kept SNP
awk 'NR==FNR{a1[$1":"$4]=$5; a2[$1":"$4]=$6; next} {k=$1":"$4; if ($5!=a1[k] || $6!=a2[k]) bad++} END{print "allele-orientation mismatches vs adm_in.bim:", bad+0}' $A/adm_in.bim sub.bim
# P rows in the same order as sub.bim
cut -f4 keep_map.tsv > rows.txt
awk 'NR==FNR{r[NR]=$1; n=NR; next} {p[FNR]=$0} END{for(i=1;i<=n;i++) print p[r[i]]}' rows.txt $A/s11/adm_in.4.P > sub.4.P.in
wc -l sub.4.P.in
admixture -P -j8 --seed=7 sub.bed 4 > admixture_proj.log 2>&1
tail -3 admixture_proj.log
wc -l sub.4.Q
echo "finished $(date)"
