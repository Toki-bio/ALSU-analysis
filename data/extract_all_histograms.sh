#!/bin/bash
# Extract ALL histogram/distribution data from PLINK QC files
# Output: tab-separated bin data for each chart
set -euo pipefail

WD="/staging/ALSU-analysis/winter2025/PLINK_301125_0312"
cd "$WD"

echo "=== STEP1: F_MISS distribution (raw, 1247 samples) ==="
echo "# bins: lo hi count"
awk 'NR>1{
  v=$6+0
  if(v<0.005) b="0.000\t0.005"
  else if(v<0.01) b="0.005\t0.010"
  else if(v<0.02) b="0.010\t0.020"
  else if(v<0.05) b="0.020\t0.050"
  else if(v<0.10) b="0.050\t0.100"
  else if(v<0.15) b="0.100\t0.150"
  else if(v<0.20) b="0.150\t0.200"
  else if(v<0.30) b="0.200\t0.300"
  else if(v<0.50) b="0.300\t0.500"
  else b="0.500\t1.000"
  a[b]++
} END{
  split("0.000\t0.005 0.005\t0.010 0.010\t0.020 0.020\t0.050 0.050\t0.100 0.100\t0.150 0.150\t0.200 0.200\t0.300 0.300\t0.500 0.500\t1.000",bins," ")
  for(i=1;i<=10;i++) print bins[i] "\t" (a[bins[i]]+0)
}' ConvSK_raw_miss.imiss

echo ""
echo "=== STEP1: removed sample IDs (first 10 + total) ==="
echo "# total:"
wc -l < remove_miss20.txt
echo "# first 10:"
head -10 remove_miss20.txt

echo ""
echo "=== STEP2: PI_HAT distribution (all pairs from IBD) ==="
echo "# bins: lo hi count"
# Use the genome file - it has all pairs
GFILE="ConvSK_mind20.genome"
if [ ! -f "$GFILE" ]; then
  echo "ERROR: $GFILE not found"
else
  awk 'NR>1{
    v=$10+0
    if(v<0.05) b="0.00\t0.05"
    else if(v<0.10) b="0.05\t0.10"
    else if(v<0.20) b="0.10\t0.20"
    else if(v<0.35) b="0.20\t0.35"
    else if(v<0.50) b="0.35\t0.50"
    else if(v<0.75) b="0.50\t0.75"
    else if(v<0.98) b="0.75\t0.98"
    else b="0.98\t1.00"
    a[b]++
  } END{
    split("0.00\t0.05 0.05\t0.10 0.10\t0.20 0.20\t0.35 0.35\t0.50 0.50\t0.75 0.75\t0.98 0.98\t1.00",bins," ")
    for(i=1;i<=8;i++) print bins[i] "\t" (a[bins[i]]+0)
  }' "$GFILE"
fi

echo ""
echo "=== STEP2: duplicate cluster summary ==="
if [ -f "dup_clusters_summary.tsv" ]; then
  echo "# First 5 clusters:"
  head -6 dup_clusters_summary.tsv
  echo "# Total clusters:"
  tail -n+2 dup_clusters_summary.tsv | wc -l
else
  echo "ERROR: dup_clusters_summary.tsv not found"
fi

echo ""
echo "=== STEP2: removed duplicate IDs (all) ==="
if [ -f "duplicates_pihat098.txt" ]; then
  echo "# total:"
  wc -l < duplicates_pihat098.txt
  echo "# all IDs:"
  cat duplicates_pihat098.txt
else
  echo "ERROR: duplicates_pihat098.txt not found"
fi

echo ""
echo "=== STEP3: MAF spectrum (raw 654k variants) ==="
echo "# bins: lo hi count"
# Use .frq file if available, otherwise compute
if [ -f "ConvSK_raw_miss.frq" ]; then
  FRQFILE="ConvSK_raw_miss.frq"
elif [ -f "ConvSK_raw.frq" ]; then
  FRQFILE="ConvSK_raw.frq"
else
  echo "# Computing MAF from raw..."
  plink --bfile ConvSK_raw --freq --out /tmp/alsu_frq 2>/dev/null
  FRQFILE="/tmp/alsu_frq.frq"
fi
echo "# using: $FRQFILE"
awk 'NR>1{
  v=$5+0; if(v<0) v=-v; if(v>0.5) v=1-v
  if(v<0.01) b="0.00\t0.01"
  else if(v<0.05) b="0.01\t0.05"
  else if(v<0.10) b="0.05\t0.10"
  else if(v<0.15) b="0.10\t0.15"
  else if(v<0.20) b="0.15\t0.20"
  else if(v<0.25) b="0.20\t0.25"
  else if(v<0.30) b="0.25\t0.30"
  else if(v<0.35) b="0.30\t0.35"
  else if(v<0.40) b="0.35\t0.40"
  else if(v<0.45) b="0.40\t0.45"
  else b="0.45\t0.50"
  a[b]++
} END{
  split("0.00\t0.01 0.01\t0.05 0.05\t0.10 0.10\t0.15 0.15\t0.20 0.20\t0.25 0.25\t0.30 0.30\t0.35 0.35\t0.40 0.40\t0.45 0.45\t0.50",bins," ")
  for(i=1;i<=11;i++) print bins[i] "\t" (a[bins[i]]+0)
}' "$FRQFILE"

echo ""
echo "=== STEP3: HWE p-value distribution (raw variants, -log10 p) ==="
echo "# bins: lo hi count"
if [ -f "ConvSK_raw_miss.hwe" ]; then
  HWEFILE="ConvSK_raw_miss.hwe"
elif [ -f "ConvSK_raw.hwe" ]; then
  HWEFILE="ConvSK_raw.hwe"
else
  echo "# Computing HWE from raw..."
  plink --bfile ConvSK_raw --hardy --out /tmp/alsu_hwe 2>/dev/null
  HWEFILE="/tmp/alsu_hwe.hwe"
fi
echo "# using: $HWEFILE"
awk 'NR>1 && $3=="ALL"{
  p=$9+0
  if(p<=0) p=1e-300
  lp=-log(p)/log(10)
  if(lp<1) b="0\t1"
  else if(lp<2) b="1\t2"
  else if(lp<3) b="2\t3"
  else if(lp<4) b="3\t4"
  else if(lp<5) b="4\t5"
  else if(lp<6) b="5\t6"
  else if(lp<10) b="6\t10"
  else b="10\tInf"
  a[b]++
} END{
  split("0\t1 1\t2 2\t3 3\t4 4\t5 5\t6 6\t10 10\tInf",bins," ")
  for(i=1;i<=8;i++) print bins[i] "\t" (a[bins[i]]+0)
}' "$HWEFILE"

echo ""
echo "=== STEP3: Per-SNP missingness distribution (raw variants) ==="
echo "# bins: lo hi count"
if [ -f "ConvSK_raw_miss.lmiss" ]; then
  LMFILE="ConvSK_raw_miss.lmiss"
else
  echo "# Computing lmiss from raw..."
  plink --bfile ConvSK_raw --missing --out /tmp/alsu_miss 2>/dev/null
  LMFILE="/tmp/alsu_miss.lmiss"
fi
echo "# using: $LMFILE"
awk 'NR>1{
  v=$5+0
  if(v<0.01) b="0.00\t0.01"
  else if(v<0.02) b="0.01\t0.02"
  else if(v<0.05) b="0.02\t0.05"
  else if(v<0.10) b="0.05\t0.10"
  else if(v<0.15) b="0.10\t0.15"
  else if(v<0.20) b="0.15\t0.20"
  else if(v<0.30) b="0.20\t0.30"
  else b="0.30\t1.00"
  a[b]++
} END{
  split("0.00\t0.01 0.01\t0.02 0.02\t0.05 0.05\t0.10 0.10\t0.15 0.15\t0.20 0.20\t0.30 0.30\t1.00",bins," ")
  for(i=1;i<=8;i++) print bins[i] "\t" (a[bins[i]]+0)
}' "$LMFILE"

echo ""
echo "=== STEP4: INFO score distribution (post-imputation) ==="
echo "# bins: lo hi count"
# Check for info files in imputation results
IMPDIR="$WD/michigan_ready_chr/imputation_results/unz"
if [ -d "$IMPDIR" ]; then
  echo "# Scanning info files in $IMPDIR"
  # Info files: chr*.info.gz or chr*.info
  cat_cmd="cat"
  infos=$(find "$IMPDIR" -name "chr*.info*" -type f 2>/dev/null | sort)
  if [ -z "$infos" ]; then
    echo "ERROR: no info files found in $IMPDIR"
  else
    echo "# files found: $(echo "$infos" | wc -l)"
    for f in $infos; do
      if [[ "$f" == *.gz ]]; then
        zcat "$f"
      else
        cat "$f"
      fi
    done | awk 'NR>1 && !/^SNP/{
      v=$7+0
      if(v<0.10) b="0.00\t0.10"
      else if(v<0.20) b="0.10\t0.20"
      else if(v<0.30) b="0.20\t0.30"
      else if(v<0.40) b="0.30\t0.40"
      else if(v<0.50) b="0.40\t0.50"
      else if(v<0.60) b="0.50\t0.60"
      else if(v<0.70) b="0.60\t0.70"
      else if(v<0.80) b="0.70\t0.80"
      else if(v<0.90) b="0.80\t0.90"
      else b="0.90\t1.01"
      a[b]++
    } END{
      split("0.00\t0.10 0.10\t0.20 0.20\t0.30 0.30\t0.40 0.40\t0.50 0.50\t0.60 0.60\t0.70 0.70\t0.80 0.80\t0.90 0.90\t1.01",bins," ")
      for(i=1;i<=10;i++) print bins[i] "\t" (a[bins[i]]+0)
    }'
  fi
else
  echo "ERROR: imputation dir not found: $IMPDIR"
fi

echo ""
echo "=== STEP15: F_ROH distribution (post-dedup, 1098 samples) ==="
# Check for ROH files
ROHDIR="$WD"
if [ -f "$ROHDIR/ConvSK_mind20_dedup.hom.indiv" ]; then
  ROHFILE="$ROHDIR/ConvSK_mind20_dedup.hom.indiv"
elif [ -f "/staging/ALSU-analysis/winter2025/3_post-imputation/UZB_v2_ROH.hom.indiv" ]; then
  ROHFILE="/staging/ALSU-analysis/winter2025/3_post-imputation/UZB_v2_ROH.hom.indiv"
else
  # Search for any .hom.indiv
  ROHFILE=$(find /staging/ALSU-analysis -name "*.hom.indiv" -type f 2>/dev/null | head -1)
fi
if [ -n "$ROHFILE" ] && [ -f "$ROHFILE" ]; then
  echo "# using: $ROHFILE"
  echo "# total samples:"
  awk 'NR>1' "$ROHFILE" | wc -l
  # Need genome length to compute F_ROH = total_KB_roh / genome_length_KB
  # Approximation: use column KBAVG * NSEG or KB column
  # Typical autosome length ~2,881,033 KB (hg19)
  # PLINK .hom.indiv cols: FID IID PHE NSEG KB KBAVG
  echo "# F_ROH = KB / 2881033"
  awk 'NR>1{
    froh = $5 / 2881033.0
    if(froh<0.005) b="0.000\t0.005"
    else if(froh<0.010) b="0.005\t0.010"
    else if(froh<0.015) b="0.010\t0.015"
    else if(froh<0.020) b="0.015\t0.020"
    else if(froh<0.025) b="0.020\t0.025"
    else if(froh<0.030) b="0.025\t0.030"
    else if(froh<0.035) b="0.030\t0.035"
    else if(froh<0.040) b="0.035\t0.040"
    else if(froh<0.050) b="0.040\t0.050"
    else if(froh<0.060) b="0.050\t0.060"
    else if(froh<0.080) b="0.060\t0.080"
    else if(froh<0.100) b="0.080\t0.100"
    else if(froh<0.150) b="0.100\t0.150"
    else b="0.150\t1.000"
    a[b]++
  } END{
    split("0.000\t0.005 0.005\t0.010 0.010\t0.015 0.015\t0.020 0.020\t0.025 0.025\t0.030 0.030\t0.035 0.035\t0.040 0.040\t0.050 0.050\t0.060 0.060\t0.080 0.080\t0.100 0.100\t0.150 0.150\t1.000",bins," ")
    for(i=1;i<=14;i++) print bins[i] "\t" (a[bins[i]]+0)
  }' "$ROHFILE"
else
  echo "ERROR: No .hom.indiv file found"
fi

echo ""
echo "=== STEP15: PI_HAT distribution (post-dedup, 885 unique samples) ==="
# This uses post-dedup IBD
POSTG="/staging/ALSU-analysis/winter2025/3_post-imputation/UZB_v2_IBD.genome"
if [ ! -f "$POSTG" ]; then
  POSTG=$(find /staging/ALSU-analysis -name "UZB_v2_IBD.genome" -type f 2>/dev/null | head -1)
fi
if [ -n "$POSTG" ] && [ -f "$POSTG" ]; then
  echo "# using: $POSTG"
  echo "# total pairs:"
  awk 'NR>1' "$POSTG" | wc -l
  awk 'NR>1{
    v=$10+0
    if(v<0.05) b="0.00\t0.05"
    else if(v<0.10) b="0.05\t0.10"
    else if(v<0.20) b="0.10\t0.20"
    else if(v<0.35) b="0.20\t0.35"
    else if(v<0.50) b="0.35\t0.50"
    else b="0.50\tmax"
    a[b]++
  } END{
    split("0.00\t0.05 0.05\t0.10 0.10\t0.20 0.20\t0.35 0.35\t0.50 0.50\tmax",bins," ")
    for(i=1;i<=6;i++) print bins[i] "\t" (a[bins[i]]+0)
  }' "$POSTG"
else
  echo "ERROR: UZB_v2_IBD.genome not found"
fi

echo ""
echo "=== DONE ==="
