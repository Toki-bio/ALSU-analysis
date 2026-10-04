#!/bin/bash
# Add Kazakh GSA genotypes (EVA PRJEB89820, GRCh37) to the NYGC + Uzbek panel (adm_in, GRCh38) at the shared SNPs.
# usage: kaz_pipeline.sh <outname> <max number of vcf files, 0 = all>
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH
export TMPDIR=/staging/tmp/scratch
R=/staging/ALSU-analysis/admixture_analysis/reference_panels/kazakh_gsa_prjeb89820
OUT=$R/merge_$1
NMAX=$2
mkdir -p "$OUT"
cd "$OUT" || exit 1
A=/staging/ALSU-analysis/spring2026/full_expanded_cohort/nygc30x_1256/adm
echo "start $(date) out=$OUT nmax=$NMAX"
ls $R/vcf_all/sample_*.snv.vcf.gz > vcf_all.list
if [ "$NMAX" -gt 0 ]; then head -n "$NMAX" vcf_all.list > vcf.list; else cp vcf_all.list vcf.list; fi
wc -l vcf.list
# 1. index + merge (missing where a sample lacks a site)
xargs -P 8 -I{} bash -c '[ -s {}.csi ] || bcftools index -f {}' < vcf.list
bcftools merge -m none --threads 8 -l vcf.list -Oz -o kaz_all.vcf.gz 2> merge.err || { tail -3 merge.err; exit 1; }
echo "merged: $(bcftools query -l kaz_all.vcf.gz | wc -l) samples, $(bcftools view -H kaz_all.vcf.gz | wc -l) sites"
# 2. plink (GRCh37), unique IDs from position and alleles
plink2 --vcf kaz_all.vcf.gz --chr 1-22 --snps-only just-acgt --max-alleles 2 --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 20 --make-bed --out kaz37 --threads 8 > kaz37.log 2>&1
tail -2 kaz37.log
# 3. liftOver GRCh37 -> GRCh38
CH=/staging/tmp/scratch/kaz_overlap/hg19ToHg38.over.chain.gz
awk '{print "chr"$1"\t"$4-1"\t"$4"\t"$2}' kaz37.bim > k37.bed
liftOver k37.bed $CH k38.bed unmapped.bed 2>/dev/null
awk '{c=$1; sub("chr","",c); print $4"\t"c"\t"$3}' k38.bed > map38.txt        # id chr pos38
# keep only variants that stay on the same chromosome
awk 'NR==FNR{c[$2]=$1; next} ($1 in c) && ($2==c[$1]) {print $1"\t"$3}' kaz37.bim map38.txt > upd_pos.txt
cut -f1 upd_pos.txt > keep_ids.txt
wc -l keep_ids.txt
plink --bfile kaz37 --extract keep_ids.txt --update-map upd_pos.txt 2 1 --make-bed --out kaz38 --threads 8 > kaz38.log 2>&1
# 4. harmonise to adm_in.bim by chr:pos and allele set; ID rename to the adm_in ID; adm_in allele orientation
awk 'NR==FNR{k=$1":"$4; a1[k]=$5; a2[k]=$6; id[k]=$2; next}
     {k=$1":"$4; if (!(k in a1)) next; x=$5; y=$6;
      if ((x==a1[k] && y==a2[k]) || (x==a2[k] && y==a1[k])) print $2"\t"id[k]"\t"a1[k]"\tkeep";
      else {cx=(x=="A"?"T":x=="T"?"A":x=="C"?"G":"C"); cy=(y=="A"?"T":y=="T"?"A":y=="C"?"G":"C");
            if ((cx==a1[k] && cy==a2[k]) || (cx==a2[k] && cy==a1[k])) print $2"\t"id[k]"\t"a1[k]"\tflip"}}' $A/adm_in.bim kaz38.bim > harm.tsv
echo "harmonised sites: $(wc -l < harm.tsv)  (flip: $(awk '$4=="flip"' harm.tsv | wc -l))"
awk '$4=="flip"{print $1}' harm.tsv > flip_ids.txt
cut -f1 harm.tsv > h_ids.txt
plink --bfile kaz38 --extract h_ids.txt --flip flip_ids.txt --make-bed --out kaz_h --threads 8 > kaz_h.log 2>&1
cut -f1,2 harm.tsv > rename.txt
cut -f2,3 harm.tsv > a1.txt
plink --bfile kaz_h --update-name rename.txt --a1-allele a1.txt 2 1 --make-bed --out kaz_h2 --threads 8 > kaz_h2.log 2>&1
# 5. Kazakh QC: SNP call rate >= 95%, sample call rate >= 95%, drop one of each pair PI_HAT > 0.2
plink --bfile kaz_h2 --geno 0.05 --mind 0.05 --make-bed --out kaz_qc1 --threads 8 > kaz_qc1.log 2>&1
grep -i "removed\|pass filters" kaz_qc1.log || true
plink --bfile kaz_qc1 --genome --min 0.2 --out kaz_rel --threads 8 > /dev/null 2>&1
awk 'NR>1{print $1,$2}' kaz_rel.genome | sort -u > kaz_rel_drop.txt
echo "relatives to drop (one per pair, first listed): $(wc -l < kaz_rel_drop.txt)"
plink --bfile kaz_qc1 --remove kaz_rel_drop.txt --make-bed --out kaz_final --threads 8 > kaz_final.log 2>&1
tail -2 kaz_final.log
# 6. merge with the existing panel restricted to the shared SNPs
cut -f2 kaz_final.bim > shared_ids.txt
plink --bfile $A/adm_in --extract shared_ids.txt --keep-allele-order --make-bed --out base_sub --threads 8 > base_sub.log 2>&1
plink --bfile base_sub --bmerge kaz_final --keep-allele-order --make-bed --out panel_kaz --threads 8 > panel_kaz.log 2>&1 || { tail -5 panel_kaz.log; exit 1; }
tail -3 panel_kaz.log
echo "panel_kaz: $(wc -l < panel_kaz.fam) samples, $(wc -l < panel_kaz.bim) SNPs"
echo "finished $(date)"
