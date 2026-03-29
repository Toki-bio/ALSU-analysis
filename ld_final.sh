#!/bin/bash
echo "=== Final: Michigan dose VCF + info + LD proxy ==="
PLINK=/staging/Tools/plink_linux_x86_64_20241022/plink
IMPDIR="/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz"

# --- 1) Check chr17.info.gz for position 80107730 ---
echo "=== 1) chr17.info.gz ==="
echo "Header:"
zcat "$IMPDIR/chr17.info.gz" 2>/dev/null | head -1
echo "--- Search 80107730 ---"
zgrep "80107730" "$IMPDIR/chr17.info.gz" 2>/dev/null | head -5
echo "--- Nearby 80107xxx ---"
zgrep "8010[67][0-9][0-9][0-9]" "$IMPDIR/chr17.info.gz" 2>/dev/null | head -20

# --- 2) Check chr17.dose.vcf.gz (no index - use zgrep) ---
echo ""
echo "=== 2) chr17.dose.vcf.gz ==="
echo "Size: $(ls -lh "$IMPDIR/chr17.dose.vcf.gz" 2>/dev/null | awk '{print $5}')"
echo "--- Search 80107730 ---"
zgrep -m5 "80107730" "$IMPDIR/chr17.dose.vcf.gz" 2>/dev/null | cut -f1-8 | head -5
echo "--- Nearby 80107500-80108000 ---"
# Use bcftools if possible, else zgrep
bcftools view "$IMPDIR/chr17.dose.vcf.gz" -r 17:80107500-80108000 -H 2>/dev/null | cut -f1-10 | head -10
if [ $? -ne 0 ]; then
  echo "bcftools failed (no index), using zgrep..."
  zgrep -m10 "8010[78]" "$IMPDIR/chr17.dose.vcf.gz" 2>/dev/null | awk -F'\t' '$2>=80107500 && $2<=80108000' | cut -f1-8 | head -10
fi

# --- 3) Check rsq quality for nearby imputed variants ---
echo ""
echo "=== 3) IMPUTATION QUALITY (chr17.rsq.keyed.tsv) ==="
head -1 "$IMPDIR/chr17.rsq.keyed.tsv" 2>/dev/null
grep "8010[67]" "$IMPDIR/chr17.rsq.keyed.tsv" 2>/dev/null | head -20

# --- 4) LD from 1000G with proxy variant ---
echo ""
echo "=== 4) LD WITH PROXY ==="
VCF17="/staging/ALSU-analysis/Fst_analysis/1000G_data/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
# Extract just ±100kb to be faster
bcftools view "$VCF17" -r 17:80007730-80207730 -Oz -o /tmp/r17b.vcf.gz 2>/dev/null
"$PLINK" --vcf /tmp/r17b.vcf.gz --make-bed --out /tmp/r17b --memory 2000 2>&1 | tail -5

echo "Variants near 80107730:"
awk '$4>=80107700 && $4<=80107800' /tmp/r17b.bim

# Pick proxy at 80107725 (5bp away from our target)
PROXY=$(awk '$4>=80107700 && $4<=80107770 {print $2; exit}' /tmp/r17b.bim)
echo "Proxy: $PROXY"

if [ -n "$PROXY" ]; then
  "$PLINK" --bfile /tmp/r17b --ld-snp "$PROXY" \
    --ld-window-kb 150 --ld-window 99999 --ld-window-r2 0 \
    --r2 --out /tmp/ldp --memory 2000 2>&1 | tail -3
  
  if [ -f /tmp/ldp.ld ]; then
    echo ""
    echo "Total LD pairs: $(wc -l < /tmp/ldp.ld)"
    echo ""
    echo "--- Top 15 LD partners (r2 > 0.2) ---"
    awk 'NR>1 && $7>=0.2 {printf "%-20s pos=%-10d r2=%.4f\n", $6, $5, $7}' /tmp/ldp.ld | sort -k3 -t= -rn | head -15
    
    echo ""
    echo "--- Cross-ref: which GSA probes are in 1000G? ---"
    RAW="/staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_120625_0107/ConvSK.bim"
    # GSA probes in ±100kb
    awk '($1==17||$1=="17") && $4>=80007730 && $4<=80207730 {print $2}' "$RAW" > /tmp/gsa.txt
    
    # Also check by position since 1000G uses different IDs
    # Get GSA positions
    awk '($1==17||$1=="17") && $4>=80007730 && $4<=80207730 {print $4}' "$RAW" > /tmp/gsa_pos.txt
    
    echo "GSA probes: $(wc -l < /tmp/gsa.txt)"
    echo ""
    echo "--- LD with GSA probes (by position match) ---"
    awk 'NR==FNR {pos[$1]=1; next} FNR>1 && ($5 in pos) {printf "pos=%-10d  1kG_id=%-20s  r2=%.4f\n", $5, $6, $7}' \
      /tmp/gsa_pos.txt /tmp/ldp.ld | sort -t= -k4 -rn | head -20
    awk 'NR==FNR {pos[$1]=1; next} FNR>1 && ($2 in pos) {printf "pos=%-10d  1kG_id=%-20s  r2=%.4f\n", $2, $3, $7}' \
      /tmp/gsa_pos.txt /tmp/ldp.ld | sort -t= -k4 -rn | head -20
    
    # Also check ConvSK_final_clean (the actual QC'd array)
    echo ""
    echo "--- LD with ConvSK_final_clean probes (post-QC array, the ones we actually have) ---"
    FCLEAN="/staging/ALSU-analysis/winter2025/3_post-imputation/ConvSK_final_clean.bim"
    awk '($1==17||$1=="17") && $4>=80007730 && $4<=80207730 {print $4}' "$FCLEAN" > /tmp/qc_pos.txt
    echo "QC'd probes in region: $(wc -l < /tmp/qc_pos.txt)"
    awk 'NR==FNR {pos[$1]=1; next} FNR>1 && ($5 in pos) {printf "pos=%-10d  1kG_id=%-20s  r2=%.4f\n", $5, $6, $7}' \
      /tmp/qc_pos.txt /tmp/ldp.ld | sort -t= -k4 -rn | head -20
    awk 'NR==FNR {pos[$1]=1; next} FNR>1 && ($2 in pos) {printf "pos=%-10d  1kG_id=%-20s  r2=%.4f\n", $2, $3, $7}' \
      /tmp/qc_pos.txt /tmp/ldp.ld | sort -t= -k4 -rn | head -20
  fi
fi

rm -f /tmp/r17b.* /tmp/ldp.* /tmp/gsa.txt /tmp/gsa_pos.txt /tmp/qc_pos.txt 2>/dev/null

echo ""
echo "=== DONE ==="
