#!/bin/bash
echo "=== LD/Imputation/Manifest check for rs5822325 (chr17:80107730 G>A) ==="
PLINK=/staging/Tools/plink_linux_x86_64_20241022/plink

# --- 1) ILLUMINA MANIFEST ---
echo "=== 1) ILLUMINA MANIFEST ==="
find /staging -maxdepth 5 -type f \( -iname "*.bpm" -o -iname "*manifest*.csv" -o -iname "*GSA*.csv" \) 2>/dev/null | head -20
echo "Also checking ALSU conversion dir:"
ls /staging/ALSU-analysis/Conversion/ 2>/dev/null
ls /staging/ALSU-analysis/Conversion/OUTPUT/ 2>/dev/null
find /staging/ALSU-analysis/Conversion -maxdepth 3 -type f \( -name "*.csv" -o -name "*.bpm" -o -name "*.txt" \) 2>/dev/null | head -20
echo "Check data/Foodomics manifest:"
head -5 /staging/data/Foodomics/GSA_192samples_Foodomics_sample_sheet.csv 2>/dev/null
echo "---"

# --- 2) IMPUTED DATA ---
echo ""
echo "=== 2) IMPUTED DATA ==="
find /staging -maxdepth 4 -type d \( -iname "*imput*" -o -iname "*beagle*" -o -iname "*minimac*" -o -iname "*phased*" \) 2>/dev/null | head -10
find /staging -maxdepth 5 -type f \( -iname "*impute*" -o -iname "*dose*vcf*" -o -iname "*info.gz" \) 2>/dev/null | head -20
echo "---"
# Check for any chr17-specific imputed files
find /staging -maxdepth 6 -type f -name "*chr17*" 2>/dev/null | grep -iv "1000G\|phase3" | head -10
echo "Check for Michigan/TOPMed imputation output:"
find /staging -maxdepth 5 -type f -name "*.dose.vcf.gz" 2>/dev/null | head -10
find /staging -maxdepth 5 -type f -name "*.info.gz" 2>/dev/null | head -10

# --- 3) ARRAY PROBES near rs5822325 ---
echo ""
echo "=== 3) NEARBY ARRAY PROBES (±100kb of 80107730) ==="
RAW="/staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_120625_0107/ConvSK.bim"
echo "Chr  rsID                          cM      BP        A1 A2"
awk '($1==17||$1=="17") && $4>=80007730 && $4<=80207730 {printf "%s  %-30s %s  %10d  %s  %s\n", $1, $2, $3, $4, $5, $6}' "$RAW"
echo ""
echo "Probes in ±100kb:"
awk '($1==17||$1=="17") && $4>=80007730 && $4<=80207730' "$RAW" | wc -l

# --- 4) 1000G VCF: extract our SNP and its neighbors ---
echo ""
echo "=== 4) 1000G chr17 VCF - rs5822325 and neighbors ==="
VCF17="/staging/ALSU-analysis/Fst_analysis/1000G_data/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
if [ -f "${VCF17}.tbi" ]; then
  echo "VCF indexed, querying..."
  echo "--- Exact position 80107730 ---"
  bcftools view "$VCF17" -r 17:80107730-80107730 -H 2>/dev/null | cut -f1-8
  echo "--- All variants 80100000-80115000 (±7.5kb tight window) ---"
  bcftools view "$VCF17" -r 17:80100000-80115000 -H 2>/dev/null | cut -f1-5 | head -30
  echo "(count):"
  bcftools view "$VCF17" -r 17:80100000-80115000 -H 2>/dev/null | wc -l
else
  echo "No tbi index, trying bcftools..."
  bcftools view "$VCF17" -r 17:80107730-80107730 -H 2>/dev/null | cut -f1-8
fi

# --- 5) LD computation using plink1.9 from 1000G VCF ---
echo ""
echo "=== 5) LD COMPUTATION ==="
# Extract region from 1000G VCF to plink format
echo "Extracting ±200kb region from 1000G VCF..."
bcftools view "$VCF17" -r 17:79907730-80307730 -Oz -o /tmp/region_chr17.vcf.gz 2>/dev/null
bcftools index -t /tmp/region_chr17.vcf.gz 2>/dev/null

# Convert to plink format
"$PLINK" --vcf /tmp/region_chr17.vcf.gz --make-bed --out /tmp/region_chr17 --memory 2000 2>/dev/null

echo "Variants in extracted region:"
wc -l /tmp/region_chr17.bim 2>/dev/null
echo "Samples:"
wc -l /tmp/region_chr17.fam 2>/dev/null

# Find our SNP in the region
echo ""
echo "--- Our target in 1000G region data ---"
awk '$4==80107730' /tmp/region_chr17.bim 2>/dev/null
TARGET_ID=$(awk '$4==80107730 {print $2; exit}' /tmp/region_chr17.bim 2>/dev/null)
if [ -n "$TARGET_ID" ]; then
  echo "Target ID: $TARGET_ID"
  
  # Get list of array probe IDs in this region
  awk '($1==17||$1=="17") && $4>=79907730 && $4<=80307730 {print $2}' "$RAW" > /tmp/array_probes.txt
  echo "Array probes in region: $(wc -l < /tmp/array_probes.txt)"
  
  # Compute LD between target and everything
  "$PLINK" --bfile /tmp/region_chr17 --ld-snp "$TARGET_ID" \
    --ld-window-kb 250 --ld-window 99999 --ld-window-r2 0 \
    --r2 --out /tmp/ld_rs5822325 --memory 2000 2>/dev/null
  
  if [ -f /tmp/ld_rs5822325.ld ]; then
    echo ""
    echo "--- Top LD partners (r² > 0.3) ---"
    awk 'NR>1 && $7>=0.3 {printf "%-25s pos=%s  r2=%.4f\n", $6, $5, $7}' /tmp/ld_rs5822325.ld | sort -t= -k3 -rn | head -20
    
    echo ""
    echo "--- LD with GSA array probes specifically ---"
    while IFS= read -r probe; do
      hit=$(awk -v p="$probe" '($3==p || $6==p) {printf "r2=%.4f  D=%.4f\n", $7, $8}' /tmp/ld_rs5822325.ld 2>/dev/null)
      if [ -n "$hit" ]; then
        pos=$(awk -v p="$probe" '($3==p) {print $2} ($6==p) {print $5}' /tmp/ld_rs5822325.ld 2>/dev/null)
        echo "  $probe (pos=$pos): $hit"
      fi
    done < /tmp/array_probes.txt
    
    echo ""
    echo "--- Best proxy among array probes ---"
    # Cross-reference LD file with array probes
    awk 'NR==FNR {probes[$1]=1; next} FNR>1 && ($6 in probes) {printf "%-30s pos=%-10s r2=%.4f\n", $6, $5, $7}' \
      /tmp/array_probes.txt /tmp/ld_rs5822325.ld | sort -t= -k3 -rn | head -10
    awk 'NR==FNR {probes[$1]=1; next} FNR>1 && ($3 in probes) {printf "%-30s pos=%-10s r2=%.4f\n", $3, $2, $7}' \
      /tmp/array_probes.txt /tmp/ld_rs5822325.ld | sort -t= -k3 -rn | head -10
  else
    echo "LD computation failed. Checking log:"
    tail -10 /tmp/ld_rs5822325.log 2>/dev/null
  fi
else
  echo "rs5822325 / pos 80107730 NOT in 1000G phase3 VCF for this region!"
  echo "Checking what's near that position:"
  awk '$4>=80107000 && $4<=80108500' /tmp/region_chr17.bim 2>/dev/null
fi

# Cleanup
rm -f /tmp/region_chr17.* /tmp/ld_rs5822325.* /tmp/array_probes.txt 2>/dev/null

echo ""
echo "=== DONE ==="
