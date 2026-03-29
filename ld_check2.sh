#!/bin/bash
echo "=== Phase 2: Manifest + Post-imputation check for rs5822325 ==="

# --- 1) GSA Manifest check ---
echo "=== 1) GSA MANIFEST (GSA-24v3-0_A2.csv) ==="
MANIFEST="/staging/ALSU-analysis/Conversion/GSA-24v3-0_A2.csv"
echo "File size: $(ls -lh "$MANIFEST" | awk '{print $5}')"
echo "--- Check manifest column format ---"
sed -n '1,10p' "$MANIFEST" | head -10
echo "..."
# Find the data header line
DATA_HEADER=$(grep -n "IlmnID\|Name\|Chr\|MapInfo\|SNP" "$MANIFEST" | head -1)
echo "Data header: $DATA_HEADER"
HEADER_LINE=$(echo "$DATA_HEADER" | cut -d: -f1)
echo "Header at line: $HEADER_LINE"
echo "--- Search for rs5822325 ---"
grep -i "rs5822325" "$MANIFEST" | head -5
echo "--- Search for position 80107730 ---"
grep "80107730" "$MANIFEST" | head -5
echo "--- Probes near 80107xxx ---"
grep "8010[67]" "$MANIFEST" | grep -i ",17," | head -10

# --- 2) Post-imputation output ---
echo ""
echo "=== 2) POST-IMPUTATION DATA ==="
echo "Directory structure:"
find /staging/ALSU-analysis/winter2025/3_post-imputation -maxdepth 2 -type f 2>/dev/null | head -40
echo "---"
ls -la /staging/ALSU-analysis/winter2025/3_post-imputation/ 2>/dev/null | head -30
echo "---"
# Also check parent
find /staging/ALSU-analysis/winter2025 -maxdepth 1 -type f -o -type d 2>/dev/null | head -20
echo "---"
# Check for dose/info files anywhere in winter2025
find /staging/ALSU-analysis/winter2025 -maxdepth 4 \( -name "*.dose.*" -o -name "*.info*" -o -name "*imputed*" -o -name "*.vcf.gz" \) -type f 2>/dev/null | head -30

# --- 3) Check all VCFs in winter2025 for rs5822325 ---
echo ""
echo "=== 3) CHECK WINTER2025 VCFs ==="
for vcf in $(find /staging/ALSU-analysis/winter2025 -maxdepth 4 -name "*.vcf.gz" -type f 2>/dev/null | head -15); do
  if [ -f "${vcf}.tbi" ]; then
    hit=$(bcftools view "$vcf" -r 17:80107730-80107730 -H 2>/dev/null)
    if [ -n "$hit" ]; then
      echo "FOUND in $vcf:"
      echo "$hit" | cut -f1-8
    fi
  fi
done

# Also try the main imputation input
echo "--- Imputation input VCF nearby probes ---"
INPUT_VCF="/staging/ALSU-analysis/winter2025/PLINK_301125_0312/impute_in.autosomes.clean.nodup.rehead.vcf.gz"
if [ -f "${INPUT_VCF}.tbi" ]; then
  echo "Querying ±3kb around 80107730..."
  bcftools view "$INPUT_VCF" -r 17:80104730-80110730 -H 2>/dev/null | cut -f1-5
fi

# --- 4) Check all winter2025 .bim files ---
echo ""
echo "=== 4) WINTER2025 PLINK FILES ==="
find /staging/ALSU-analysis/winter2025 -maxdepth 4 -name "*.bim" -type f 2>/dev/null | while read bim; do
  echo "File: $bim ($(wc -l < "$bim") SNPs)"
  hit=$(awk '$2=="rs5822325" || $4==80107730' "$bim" 2>/dev/null)
  if [ -n "$hit" ]; then
    echo "  FOUND: $hit"
  else
    nearby=$(awk '($1==17||$1=="17") && $4>=80105000 && $4<=80110000 {printf "  %s %d %s %s\n", $2, $4, $5, $6}' "$bim" 2>/dev/null)
    [ -n "$nearby" ] && echo "  Nearby:
$nearby"
  fi
done

# --- 5) SNP list ---
echo ""
echo "=== 5) SNP_LIST.TXT ==="
SNP_LIST="/staging/ALSU-analysis/Conversion/snp_list.txt"
if [ -f "$SNP_LIST" ]; then
  echo "$(wc -l < "$SNP_LIST") lines"
  grep -i "rs5822325\|80107730" "$SNP_LIST" | head -5
  echo "Sample:"
  head -3 "$SNP_LIST"
fi

echo ""
echo "=== DONE ==="
