#!/bin/bash
echo "=== Targeted search for rs5822325 / chr17:80107730 ==="

# Known .bim files from previous scan
BIMS=(
  /staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_120625_0107/ConvSK.bim
  /staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_sept2025/ConvSK.bim
  /staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_sept2025_1/ConvSK_qc1_geno_mind.bim
  /staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_sept2025_2/ConvSK.bim
  /staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_sept2025_3/ConvSK_1.bim
  /staging/Sergei/ALSU/ConvSK/PLINK_120625_0107/alsu_qc_00.bim
  /staging/ALSU-analysis/Fst_analysis/1000G_data/1000G_EUR_chr17.bim
)

echo "--- 1) Exact rsID match ---"
for f in "${BIMS[@]}"; do
  [ -f "$f" ] || continue
  hit=$(awk '$2=="rs5822325"' "$f")
  [ -n "$hit" ] && echo "FOUND $f: $hit"
done

echo "--- 2) Exact position chr17:80107730 ---"
for f in "${BIMS[@]}"; do
  [ -f "$f" ] || continue
  hit=$(awk '($1==17||$1=="17") && $4==80107730' "$f")
  [ -n "$hit" ] && echo "FOUND $f: $hit"
done

echo "--- 3) Neighborhood +/-500kb in raw array ---"
RAW=/staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_120625_0107/ConvSK.bim
echo "Total SNPs: $(wc -l < "$RAW")"
echo "Chr17 SNPs in 79.6M-80.6M:"
awk '($1==17||$1=="17") && $4>=79607730 && $4<=80607730' "$RAW"

echo "--- 4) Check 1000G chr17 VCF if available ---"
VCF17=/staging/ALSU-analysis/Fst_analysis/1000G_data/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
if [ -f "$VCF17" ]; then
  echo "Found 1000G chr17 VCF"
  tabix "$VCF17" 17:80107730-80107730 2>/dev/null || \
  zgrep -m5 "80107730" "$VCF17" 2>/dev/null | head -3
else
  echo "No 1000G chr17 VCF found, checking what chr17 files exist:"
  ls /staging/ALSU-analysis/Fst_analysis/1000G_data/*chr17* 2>/dev/null
fi

echo "--- 5) Check ALSU VCFs for chr17:80107730 ---"
for vcf in /staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_sept2025/ConvSK_mind90.vcf.gz \
           /staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_sept2025/ConvSK_mind90_clean.vcf.gz; do
  [ -f "$vcf" ] || continue
  echo "Checking $vcf..."
  tabix "$vcf" 17:80107730-80107730 2>/dev/null || \
  zgrep -m3 "80107730" "$vcf" 2>/dev/null | head -2
done

echo "=== DONE ==="
