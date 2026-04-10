#!/bin/bash
set -euo pipefail
VCFDIR="/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz"

echo "### BCFTOOLS_CHECK ###"
which bcftools 2>/dev/null && bcftools --version | head -1 || echo "NO_BCFTOOLS"

echo "### INFO_BINS ###"
# Extract R2 from VCF INFO field and bin - use zgrep+sed for speed
for chr in $(seq 1 22); do
  zcat "$VCFDIR/chr${chr}.info.gz" 2>/dev/null | awk '!/^#/{
    n=split($8,a,";")
    for(i=1;i<=n;i++) if(a[i]~/^R2=/) {sub("R2=","",a[i]); print a[i]; break}
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

echo "### DONE ###"
