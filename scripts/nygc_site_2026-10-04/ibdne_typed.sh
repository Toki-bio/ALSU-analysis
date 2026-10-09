#!/bin/bash
# Recent Ne from IBD on TYPED array SNPs (Beagle phasing, hap-ibd, IBDNe). The imputed-data run failed sanity checks.
# usage: ibdne_typed.sh <outname> "<chr list>"
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH
export TMPDIR=/staging/tmp/scratch
OUT=/staging/ALSU-analysis/admixture_analysis/pophistory_typed_$1
CHRS="$2"
mkdir -p "$OUT"
cd "$OUT" || exit 1
echo "start $(date) out=$OUT chrs=[$CHRS]"
T=/staging/tmp/scratch/ibdne_tools
[ -s $T/beagle.jar ] || curl -sS -fL -m 300 -o $T/beagle.jar https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar
M=$T/no_chr_in_chrom_field
F=/staging/ALSU-analysis/spring2026/full_expanded_cohort
# final cohort minus second-degree relatives (same list as the imputed-data run)
KEEP=/staging/ALSU-analysis/admixture_analysis/pophistory_full/keep_ids.txt
awk 'NR==FNR{k[$1]=1; next} ($2 in k){print $1,$2}' $KEEP $F/FULL_QC_FINAL.fam > keep_fidiid.txt
echo "typed samples kept: $(wc -l < keep_fidiid.txt) of $(wc -l < $KEEP)"
run_chr() {
  C=$1
  plink2 --bfile $F/FULL_QC_FINAL --keep keep_fidiid.txt --chr $C --maf 0.01 --geno 0.05 --export vcf bgz id-paste=iid --out typed_chr$C --threads 2 > typed_chr$C.plink.log 2>&1
  java -Xmx16g -jar $T/beagle.jar gt=typed_chr$C.vcf.gz map=$M/plink.chr$C.GRCh38.map out=ph_chr$C nthreads=4 impute=false > beagle_chr$C.log 2>&1
  java -Xmx16g -jar $T/hap-ibd.jar gt=ph_chr$C.vcf.gz map=$M/plink.chr$C.GRCh38.map out=hap_chr$C min-seed=2 min-output=2 nthreads=4 > hap_chr$C.log 2>&1
  echo "chr$C: markers $(bcftools index -n ph_chr$C.vcf.gz 2>/dev/null || echo ?) ; segments $(zcat hap_chr$C.ibd.gz | wc -l)"
  rm -f typed_chr$C.vcf.gz
}
export -f run_chr
export T M F
printf '%s\n' $CHRS | xargs -P 2 -I{} bash -c 'run_chr {}'
cat hap_chr*.ibd.gz > all.ibd.gz
for C in $CHRS; do cat $M/plink.chr$C.GRCh38.map; done > all.map
NSEG=$(zcat all.ibd.gz | wc -l); NS=$(wc -l < keep_fidiid.txt)
echo "ibd segments total: $NSEG ; per pair: $(awk -v n=$NSEG -v s=$NS 'BEGIN{printf "%.3f", n/(s*(s-1)/2)}')"
if [ "$(echo $CHRS | wc -w)" -lt 10 ]; then echo "toy run: stop before IBDNe"; echo "finished $(date)"; exit 0; fi
zcat all.ibd.gz | java -Xmx40g -jar $T/ibdne.jar map=all.map out=uzb mincm=2 nthreads=8 > ibdne.log 2>&1 || { tail -5 ibdne.log; exit 1; }
tail -3 ibdne.log
echo "finished $(date)"
