#!/bin/bash
# Extract ALL histogram data - fixed version using numeric bin indices
set -euo pipefail

WD="/staging/ALSU-analysis/winter2025/PLINK_301125_0312"
cd "$WD"

echo "### STEP1_FMISS ###"
awk 'NR>1{
  v=$6+0
  if(v<0.005) b=1
  else if(v<0.01) b=2
  else if(v<0.02) b=3
  else if(v<0.05) b=4
  else if(v<0.10) b=5
  else if(v<0.15) b=6
  else if(v<0.20) b=7
  else if(v<0.30) b=8
  else if(v<0.50) b=9
  else b=10
  a[b]++
} END{
  for(i=1;i<=10;i++) printf "%d\n", a[i]+0
}' ConvSK_raw_miss.imiss

echo "### STEP1_IDS ###"
cat remove_miss20.txt

echo "### STEP2_GENOME_FIND ###"
find /staging/ALSU-analysis -maxdepth 4 -name "*.genome" -type f 2>/dev/null

echo "### STEP2_DUP_IDS_FIND ###"
find /staging/ALSU-analysis -maxdepth 4 -name "duplicates*" -type f 2>/dev/null
find /staging/ALSU-analysis -maxdepth 4 -name "remove_dup*" -type f 2>/dev/null

echo "### STEP2_CLUSTERS ###"
cat dup_clusters_summary.tsv

echo "### STEP2_PIHAT ###"
# Search for the genome file
GFILE=$(find /staging/ALSU-analysis -maxdepth 5 -name "ConvSK_mind20.genome" -type f 2>/dev/null | head -1)
if [ -z "$GFILE" ]; then
  # Try alternative names
  GFILE=$(find "$WD" -maxdepth 2 -name "*.genome" -type f 2>/dev/null | head -1)
fi
if [ -n "$GFILE" ] && [ -f "$GFILE" ]; then
  echo "FILE:$GFILE"
  awk 'NR>1{
    v=$10+0
    if(v<0.05) b=1
    else if(v<0.10) b=2
    else if(v<0.20) b=3
    else if(v<0.35) b=4
    else if(v<0.50) b=5
    else if(v<0.75) b=6
    else if(v<0.98) b=7
    else b=8
    a[b]++
  } END{
    for(i=1;i<=8;i++) printf "%d\n", a[i]+0
  }' "$GFILE"
else
  echo "NOT_FOUND"
  # Try step2 subdirectory
  echo "LISTING:"
  ls -la "$WD"/step2_*/ 2>/dev/null || echo "no step2 dir"
  find "$WD" -maxdepth 2 -name "*ibd*" -o -name "*genome*" -o -name "*pihat*" 2>/dev/null
fi

echo "### STEP3_MAF ###"
awk 'NR>1{
  v=$5+0; if(v<0) v=-v; if(v>0.5) v=1-v
  if(v<0.01) b=1
  else if(v<0.05) b=2
  else if(v<0.10) b=3
  else if(v<0.15) b=4
  else if(v<0.20) b=5
  else if(v<0.25) b=6
  else if(v<0.30) b=7
  else if(v<0.35) b=8
  else if(v<0.40) b=9
  else if(v<0.45) b=10
  else b=11
  a[b]++
} END{
  for(i=1;i<=11;i++) printf "%d\n", a[i]+0
}' /tmp/alsu_frq.frq

echo "### STEP3_HWE ###"
awk 'NR>1 && $3=="ALL"{
  p=$9+0
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

echo "### STEP3_LMISS ###"
awk 'NR>1{
  v=$5+0
  if(v<0.01) b=1
  else if(v<0.02) b=2
  else if(v<0.05) b=3
  else if(v<0.10) b=4
  else if(v<0.15) b=5
  else if(v<0.20) b=6
  else if(v<0.30) b=7
  else b=8
  a[b]++
} END{
  for(i=1;i<=8;i++) printf "%d\n", a[i]+0
}' ConvSK_raw_miss.lmiss

echo "### STEP4_INFO ###"
IMPDIR="$WD/michigan_ready_chr/imputation_results/unz"
for f in $(find "$IMPDIR" -name "chr*.info.gz" -type f 2>/dev/null | sort); do
  zcat "$f"
done | awk '!/^SNP/ && NR>1{
  v=$7+0
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

echo "### STEP15_ROH_FIND ###"
find /staging/ALSU-analysis -maxdepth 5 -name "*.hom.indiv" -type f 2>/dev/null

echo "### STEP15_PIHAT_POST ###"
# Post-dedup IBD — use the v2 file
POSTG=$(find /staging/ALSU-analysis -maxdepth 5 -name "UZB_v2_IBD.genome" -type f 2>/dev/null | head -1)
if [ -n "$POSTG" ] && [ -f "$POSTG" ]; then
  echo "FILE:$POSTG"
  echo "PAIRS:"
  awk 'NR>1' "$POSTG" | wc -l
  awk 'NR>1{
    v=$10+0
    if(v<0.05) b=1
    else if(v<0.10) b=2
    else if(v<0.20) b=3
    else if(v<0.35) b=4
    else if(v<0.50) b=5
    else b=6
    a[b]++
  } END{
    for(i=1;i<=6;i++) printf "%d\n", a[i]+0
  }' "$POSTG"
else
  echo "NOT_FOUND"
fi

echo "### DONE ###"
