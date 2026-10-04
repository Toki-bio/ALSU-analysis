#!/bin/bash
# Recent effective population size of the Uzbek cohort from IBD segments (hap-ibd + IBDNe, Browning & Browning).
# usage: ibdne.sh <outname> "<chr list>"
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH
export TMPDIR=/staging/tmp/scratch
OUT=/staging/ALSU-analysis/admixture_analysis/pophistory_$1
CHRS="$2"
mkdir -p "$OUT"
cd "$OUT" || exit 1
echo "start $(date) out=$OUT chrs=[$CHRS]"
T=/staging/tmp/scratch/ibdne_tools
mkdir -p $T
[ -s $T/hap-ibd.jar ] || curl -sS -fL -m 300 -o $T/hap-ibd.jar https://faculty.washington.edu/browning/hap-ibd.jar
[ -s $T/ibdne.jar ] || curl -sS -fL -m 300 -o $T/ibdne.jar https://faculty.washington.edu/browning/ibdne/ibdne.23Apr20.ae9.jar
if [ ! -s $T/chr_in_chrom_field/plink.chrchr1.GRCh38.map ]; then
  curl -sS -fL -m 300 -o $T/maps.zip https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
  (cd $T && unzip -oq maps.zip)
fi
ls $T | head -8
V=/staging/ALSU-analysis/spring2026/full_expanded_cohort/imputation_results/hq_filtered
# samples: final cohort (1,256) minus one of each pair with PI_HAT > 0.177 (second degree or closer)
G=/staging/ALSU-analysis/spring2026/full_expanded_cohort/step15_roh_ibd/UZB_expanded_IBD.genome
awk '{print $2}' $V/UZB_imputed_HQ_qc.fam | sort -u > cohort_ids.txt
awk 'NR>1 && $10>0.177 {print $4}' $G | sort -u > drop_rel.txt
comm -23 cohort_ids.txt drop_rel.txt > keep_ids.txt
echo "cohort $(wc -l < cohort_ids.txt), relatives dropped $(wc -l < drop_rel.txt), kept $(wc -l < keep_ids.txt)"
run_chr() {
  C=$1
  bcftools view -S keep_ids.txt --force-samples -m2 -M2 -v snps -i 'INFO/MAF>=0.05' $V/chr$C.HQ.vcf.gz -Oz -o chr$C.vcf.gz 2> chr$C.view.err
  java -Xmx24g -jar $T/hap-ibd.jar gt=chr$C.vcf.gz map=$T/chr_in_chrom_field/plink.chrchr$C.GRCh38.map out=hap_chr$C min-seed=2 min-output=2 nthreads=6 > hap_chr$C.log 2>&1
  echo "chr$C: $(bcftools index -n chr$C.vcf.gz 2>/dev/null || true) markers; segments $(zcat hap_chr$C.ibd.gz | wc -l)"
  rm -f chr$C.vcf.gz chr$C.vcf.gz.csi
}
export -f run_chr
export T V
printf '%s\n' $CHRS | xargs -P 5 -I{} bash -c 'run_chr {}'
cat hap_chr*.ibd.gz > all.ibd.gz
for C in $CHRS; do cat $T/chr_in_chrom_field/plink.chrchr$C.GRCh38.map; done > all.map
echo "ibd segments total: $(zcat all.ibd.gz | wc -l)"
zcat all.ibd.gz | java -Xmx40g -jar $T/ibdne.jar map=all.map out=uzb mincm=2 nthreads=16 > ibdne.log 2>&1 || { tail -5 ibdne.log; exit 1; }
tail -5 ibdne.log
ls
echo "finished $(date)"
