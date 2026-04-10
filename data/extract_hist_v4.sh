#!/bin/bash
set -euo pipefail
WD="/staging/ALSU-analysis/winter2025/PLINK_301125_0312"

echo "### FIND_INFO ###"
# Michigan Imputation Server creates .info.gz files - find them
find "$WD/michigan_ready_chr/imputation_results" -maxdepth 3 -type f -name "*.info*" 2>/dev/null | head -20
echo "---"
# Also check for dose files
find "$WD/michigan_ready_chr/imputation_results" -maxdepth 3 -type f 2>/dev/null | head -40
echo "---"
# Check unz directory specifically
ls -la "$WD/michigan_ready_chr/imputation_results/unz/" 2>/dev/null | head -30

echo "### DUP_IDS ###"
cat "$WD/remove_dups_pihat098.txt"

echo "### F_ROH ###"
ROHFILE="/staging/ALSU-analysis/admixture_analysis/roh_analysis/UZB_ROH.hom.indiv"
head -2 "$ROHFILE"
echo "COLS:"
head -1 "$ROHFILE" | awk '{for(i=1;i<=NF;i++) print i, $i}'
echo "N:"
awk 'NR>1' "$ROHFILE" | wc -l
echo "DATA:"
awk 'NR>1{
  froh = $5 / 2881033.0
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

echo "### INFO_SEARCH ###"
# The .info.gz in unz/ is actually VCF. Look for proper info files
find "$WD/michigan_ready_chr" -maxdepth 4 -type f -name "*.info" 2>/dev/null | head -10
echo "---"
# Check if there are info files alongside dose vcfs in the original download
find "$WD" -maxdepth 5 -type f -name "chr*.info.gz" 2>/dev/null | while read f; do
  echo "FILE: $f"
  zcat "$f" | head -1
  echo "---"
done 2>/dev/null | head -40

echo "### INFO_ALT ###"
# Try reading Rsq from VCF INFO field instead
# Michigan imputation puts R2/DR2 in the VCF INFO field
VCFDIR="$WD/michigan_ready_chr/imputation_results/unz"
FIRST_VCF=$(ls "$VCFDIR"/chr1.*.gz 2>/dev/null | head -1)
if [ -z "$FIRST_VCF" ]; then
  # Maybe the .info.gz ARE the VCFs
  FIRST_VCF="$VCFDIR/chr1.info.gz"
fi
echo "VCF: $FIRST_VCF"
zcat "$FIRST_VCF" | grep -v "^##" | head -3 | cut -f1-9

echo "### INFO_R2_EXTRACT ###"
# Extract R2 from VCF INFO field for all chromosomes
for f in $(ls "$VCFDIR"/chr*.info.gz 2>/dev/null | sort -V); do
  zcat "$f" | grep -v "^#" | awk '{
    split($8, info, ";")
    for(i in info) {
      if(info[i] ~ /^R2=/) {
        split(info[i], kv, "=")
        print kv[2]
      }
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
  printf "TOTAL:%d\n", t
  for(i=1;i<=10;i++) printf "%d\n", a[i]+0
}'

echo "### DONE ###"
