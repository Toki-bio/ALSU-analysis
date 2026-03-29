#!/bin/bash
echo "=== LD proxy fix: rename variants to chr:pos, compute properly ==="
PLINK=/staging/Tools/plink_linux_x86_64_20241022/plink
VCF17="/staging/ALSU-analysis/Fst_analysis/1000G_data/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

# Extract ±100kb
bcftools view "$VCF17" -r 17:80007730-80207730 -Oz -o /tmp/r17f.vcf.gz 2>/dev/null
"$PLINK" --vcf /tmp/r17f.vcf.gz --make-bed --out /tmp/r17f --memory 2000 2>/dev/null

# Rename all "." IDs to chr:pos format
awk '{if($2==".") $2=$1":"$4; print}' OFS='\t' /tmp/r17f.bim > /tmp/r17f_renamed.bim
cp /tmp/r17f.bed /tmp/r17f_renamed.bed
cp /tmp/r17f.fam /tmp/r17f_renamed.fam

# Our proxy at 80107725
PROXY="17:80107725"
echo "Proxy variant: $PROXY"
grep "80107725" /tmp/r17f_renamed.bim

# Compute LD
"$PLINK" --bfile /tmp/r17f_renamed --ld-snp "$PROXY" \
  --ld-window-kb 150 --ld-window 99999 --ld-window-r2 0 \
  --r2 --out /tmp/ldfix --memory 2000 2>/dev/null

echo "LD pairs: $(wc -l < /tmp/ldfix.ld)"

# GSA probe positions
RAW="/staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_120625_0107/ConvSK.bim"
FCLEAN="/staging/ALSU-analysis/winter2025/3_post-imputation/ConvSK_final_clean.bim"

# Get all GSA positions in ±100kb (raw array)
awk '($1==17||$1=="17") && $4>=80007730 && $4<=80207730 {print $4}' "$RAW" > /tmp/gsa_allpos.txt
# Get QC'd positions
awk '($1==17||$1=="17") && $4>=80007730 && $4<=80207730 {print $4}' "$FCLEAN" > /tmp/gsa_qcpos.txt

echo ""
echo "=== LD RESULTS: Proxy 17:80107725 (5bp from rs5822325) ==="
echo ""
echo "--- Top 20 LD partners (r2>0.1) ---"
awk 'NR>1 && $7>=0.1 {printf "%-20s pos=%-10s r2=%.4f\n", $6, $5, $7}' /tmp/ldfix.ld | sort -t= -k3 -rn | head -20
echo ""
echo "--- All LD WITH positions matching RAW ARRAY probes ---"
echo "(These are GSA probes the genotyping chip actually measured)"
awk 'NR==FNR {p[$1]=1;next} FNR>1 && ($5 in p) {printf "1kG=%s  arrayPos=%-10s r2=%.4f  D'=%.4f\n",$6,$5,$7,$8}' \
  /tmp/gsa_allpos.txt /tmp/ldfix.ld | sort -t= -k4 -rn | head -30
echo "---reverse direction---"
awk 'NR==FNR {p[$1]=1;next} FNR>1 && ($2 in p) {printf "1kG=%s  arrayPos=%-10s r2=%.4f\n",$3,$2,$7}' \
  /tmp/gsa_allpos.txt /tmp/ldfix.ld | sort -t= -k4 -rn | head -30

echo ""
echo "--- LD WITH positions matching QC'd ARRAY probes ---"
echo "(These survived quality control and are in ConvSK_final_clean)"
awk 'NR==FNR {p[$1]=1;next} FNR>1 && ($5 in p) {printf "1kG=%s  qcPos=%-10s r2=%.4f\n",$6,$5,$7}' \
  /tmp/gsa_qcpos.txt /tmp/ldfix.ld | sort -t= -k4 -rn | head -30
awk 'NR==FNR {p[$1]=1;next} FNR>1 && ($2 in p) {printf "1kG=%s  qcPos=%-10s r2=%.4f\n",$3,$2,$7}' \
  /tmp/gsa_qcpos.txt /tmp/ldfix.ld | sort -t= -k4 -rn | head -30

echo ""
echo "--- Specific closest probes ---"
for pos in 80107553 80107574 80107590 80107599 80107630 80107858 80108537 80109927; do
  r2=$(awk -v p="$pos" 'NR>1 && ($5==p || $2==p) {printf "%.4f", $7}' /tmp/ldfix.ld)
  [ -n "$r2" ] && echo "  pos=$pos r2=$r2" || echo "  pos=$pos NOT_IN_1000G"
done

# Also check the .info.gz by position (bcftools query)
echo ""
echo "=== MICHIGAN IMPUTATION - exact position search ==="
INFO="/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz/chr17.info.gz"
echo "Trying bcftools on info.gz..."
bcftools view "$INFO" -r chr17:80107730-80107730 -H 2>/dev/null | head -3
echo "Trying 17:80107730..."
bcftools view "$INFO" -r 17:80107730-80107730 -H 2>/dev/null | head -3
echo "Nearby in info:"
bcftools view "$INFO" -r chr17:80107000-80108500 -H 2>/dev/null | cut -f1-8 | head -15
bcftools view "$INFO" -r 17:80107000-80108500 -H 2>/dev/null | cut -f1-8 | head -15

rm -f /tmp/r17f.* /tmp/r17f_renamed.* /tmp/ldfix.* /tmp/gsa_allpos.txt /tmp/gsa_qcpos.txt 2>/dev/null

echo ""
echo "=== DONE ==="
