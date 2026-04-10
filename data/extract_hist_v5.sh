#!/bin/bash
set -euo pipefail
WD="/staging/ALSU-analysis/winter2025/PLINK_301125_0312"

echo "### DUP_IDS ###"
cat "$WD/remove_dups_pihat098.txt"

echo "### F_ROH ###"
ROHFILE="/staging/ALSU-analysis/admixture_analysis/roh_analysis/UZB_ROH.hom.indiv"
echo "HEAD:"
head -2 "$ROHFILE"
echo "N:"
awk 'NR>1' "$ROHFILE" | wc -l
echo "DATA:"
awk 'NR>1{
  froh=$5/2881033.0
  if(froh<0.005) b=1
  else if(froh<0.010) b=2
  else if(froh<0.015) b=3
  else if(froh<0.020) b=4
  else if(froh<0.025) b=5
  else if(froh<0.030) b=6
  else if(froh<0.035) b=7
  else if(froh<0.040) b=8
  else if(froh<0.050) b=9
  else if(froh<0.060) b=10
  else if(froh<0.080) b=11
  else if(froh<0.100) b=12
  else if(froh<0.150) b=13
  else b=14
  a[b]++
} END{
  for(i=1;i<=14;i++) printf "%d\n", a[i]+0
}' "$ROHFILE"

echo "### INFO_FORMAT ###"
VCFDIR="$WD/michigan_ready_chr/imputation_results/unz"
# Check one VCF INFO field format
zcat "$VCFDIR/chr22.info.gz" | grep -v "^#" | head -1 | cut -f8 | tr ';' '\n'

echo "### INFO_DATA ###"
# Use bcftools if available, otherwise awk on smallest chr first
if command -v bcftools &>/dev/null; then
  echo "METHOD:bcftools"
  for chr in $(seq 1 22); do
    bcftools query -f '%INFO/R2\n' "$VCFDIR/chr${chr}.info.gz" 2>/dev/null || \
    bcftools query -f '%INFO/DR2\n' "$VCFDIR/chr${chr}.info.gz" 2>/dev/null || true
  done | awk '{
    v=$1+0
    if(v<0.10) b=1
    else if(v<0.20) b=2
    else if(v<0.30) b=3
    else if(v<0.40) b=4
    else if(v<0.50) b=5
    else if(v<0.60) b=6
    else if(v<0.70) b=7
    else if(v<0.80) b=8
    else if(v<0.90) b=9
    else b=10
    a[b]++
  } END{
    t=0; for(i=1;i<=10;i++) t+=a[i]+0
    printf "TOTAL:%d\n",t
    for(i=1;i<=10;i++) printf "%d\n", a[i]+0
  }'
else
  echo "METHOD:awk"
  for chr in $(seq 1 22); do
    zcat "$VCFDIR/chr${chr}.info.gz" 2>/dev/null | grep -v "^#" | awk '{
      n=split($8,info,";")
      for(i=1;i<=n;i++){
        if(info[i]~/^R2=/){sub("R2=","",info[i]); print info[i]}
        else if(info[i]~/^DR2=/){sub("DR2=","",info[i]); print info[i]}
      }
    }'
  done | awk '{
    v=$1+0
    if(v<0.10) b=1
    else if(v<0.20) b=2
    else if(v<0.30) b=3
    else if(v<0.40) b=4
    else if(v<0.50) b=5
    else if(v<0.60) b=6
    else if(v<0.70) b=7
    else if(v<0.80) b=8
    else if(v<0.90) b=9
    else b=10
    a[b]++
  } END{
    t=0; for(i=1;i<=10;i++) t+=a[i]+0
    printf "TOTAL:%d\n",t
    for(i=1;i<=10;i++) printf "%d\n", a[i]+0
  }'
fi

echo "### DONE ###"
