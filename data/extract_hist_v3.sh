#!/bin/bash
set -euo pipefail
WD="/staging/ALSU-analysis/winter2025/PLINK_301125_0312"
cd "$WD"

echo "### HWE_HEAD ###"
head -3 /tmp/alsu_hwe.hwe

echo "### HWE_COLS ###"
awk 'NR==1{for(i=1;i<=NF;i++) print i, $i}' /tmp/alsu_hwe.hwe

echo "### HWE_DATA ###"
# Try without TEST filter, use last column as P
awk 'NR>1{
  # Find the p-value - last column that looks numeric
  p=$NF+0
  if(p<=0) p=1e-300
  lp=-log(p)/log(10)
  if(lp<1) b=1
  else if(lp<2) b=2
  else if(lp<3) b=3
  else if(lp<4) b=4
  else if(lp<5) b=5
  else if(lp<6) b=6
  else if(lp<10) b=7
  else b=8
  a[b]++
} END{
  for(i=1;i<=8;i++) printf "%d\n", a[i]+0
}' /tmp/alsu_hwe.hwe

echo "### INFO_HEAD ###"
IMPDIR="$WD/michigan_ready_chr/imputation_results/unz"
FIRST_INFO=$(find "$IMPDIR" -name "chr1.info.gz" -type f 2>/dev/null | head -1)
if [ -z "$FIRST_INFO" ]; then
  FIRST_INFO=$(find "$IMPDIR" -name "*.info*" -type f 2>/dev/null | sort | head -1)
fi
echo "FILE: $FIRST_INFO"
zcat "$FIRST_INFO" | head -3

echo "### INFO_COLS ###"
zcat "$FIRST_INFO" | head -1 | awk '{for(i=1;i<=NF;i++) print i, $i}'

echo "### INFO_COUNT ###"
# Count total lines across all info files (minus headers)
TOTAL=0
for f in $(find "$IMPDIR" -name "chr*.info.gz" -type f 2>/dev/null | sort); do
  N=$(zcat "$f" | tail -n+2 | wc -l)
  echo "$(basename $f): $N"
  TOTAL=$((TOTAL + N))
done
echo "TOTAL: $TOTAL"

echo "### INFO_DATA ###"
# Determine Rsq column from header, then compute bins
RSQCOL=$(zcat "$FIRST_INFO" | head -1 | awk '{for(i=1;i<=NF;i++) if($i=="Rsq") print i}')
echo "Rsq column: $RSQCOL"
for f in $(find "$IMPDIR" -name "chr*.info.gz" -type f 2>/dev/null | sort); do
  zcat "$f" | tail -n+2
done | awk -v col="$RSQCOL" '{
  v=$col+0
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
  for(i=1;i<=10;i++) printf "%d\n", a[i]+0
}'

echo "### DUP_IDS ###"
cat "$WD/remove_dups_pihat098.txt"

echo "### F_ROH ###"
ROHFILE="/staging/ALSU-analysis/admixture_analysis/roh_analysis/UZB_ROH.hom.indiv"
echo "HEAD:"
head -3 "$ROHFILE"
echo "COLS:"
awk 'NR==1{for(i=1;i<=NF;i++) print i, $i}' "$ROHFILE"
echo "NSAMPLES:"
awk 'NR>1' "$ROHFILE" | wc -l
echo "DATA:"
# F_ROH = KB / 2881033 (autosomal genome length in KB for hg19)
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

echo "### DONE ###"
