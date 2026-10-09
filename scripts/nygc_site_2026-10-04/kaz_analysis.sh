#!/bin/bash
# PCA, ADMIXTURE and FST on the merged Uzbek + NYGC + Kazakh panel (low parallelism: DRAGEN had a CPU fault 2026-10-07).
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH
export TMPDIR=/staging/tmp/scratch
R=/staging/ALSU-analysis/admixture_analysis/reference_panels/kazakh_gsa_prjeb89820/merge_full
K=/staging/ALSU-analysis/spring2026/full_expanded_cohort/pbs_refqc_2026-10/nygc30x_1256
cd "$R" || exit 1
echo "start $(date)"
wc -l panel_kaz.fam panel_kaz.bim
# population labels: cohort FID 0 = UZB; NYGC samples from the keep lists; Kazakh = rest
cut -d' ' -f2 kaz_final.fam > kaz_ids.txt
awk -v d="$K" 'BEGIN{split("EUR EAS SAS AFR",P," "); for(i in P){f=d"/keep_"P[i]".txt"; while((getline l < f)>0){split(l,x," "); lab[x[2]]=P[i]}}; while((getline l < "kaz_ids.txt")>0) kz[l]=1}
     {if($2 in kz) print $1,$2,"KAZ"; else if($2 in lab) print $1,$2,lab[$2]; else print $1,$2,"UZB"}' panel_kaz.fam > labels.txt
cut -d' ' -f3 labels.txt | sort | uniq -c
# PCA
plink --bfile panel_kaz --pca 20 --out pca_kaz --threads 4 > pca_kaz.plink.log 2>&1
echo "pca done $(date)"
# FST, all pairs among six populations (same SNP set)
pops="UZB KAZ EUR EAS SAS AFR"
echo -e "pair\tn1\tn2\tmean_fst\tweighted_fst" > fst_pairs.tsv
for a in $pops; do for b in $pops; do
  [[ "$a" < "$b" ]] || continue
  awk -v a=$a -v b=$b '$3==a||$3==b' labels.txt > within_${a}_${b}.txt
  cut -d' ' -f1,2 within_${a}_${b}.txt > keep_${a}_${b}.txt
  plink --bfile panel_kaz --keep keep_${a}_${b}.txt --within within_${a}_${b}.txt --fst --out fst_${a}_${b} --threads 2 > fst_${a}_${b}.out 2>&1
  m=$(grep "Mean Fst" fst_${a}_${b}.log | awk '{print $NF}'); w=$(grep "Weighted Fst" fst_${a}_${b}.log | awk '{print $NF}')
  echo -e "${a}_${b}\t$(awk -v x=$a '$3==x' labels.txt | wc -l)\t$(awk -v x=$b '$3==x' labels.txt | wc -l)\t$m\t$w" >> fst_pairs.tsv
done; done
cat fst_pairs.tsv
echo "fst done $(date)"
# ADMIXTURE K=2..7, 5-fold CV, one seed
for k in 2 3 4 5 6 7; do
  admixture --cv=5 -j4 -s 11 panel_kaz.bed $k > adm_K$k.log 2>&1
  echo "K=$k $(grep 'CV error' adm_K$k.log)"
done
echo "finished $(date)"
