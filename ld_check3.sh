#!/bin/bash
echo "=== Phase 3: Imputed data + LD for rs5822325 ==="
PLINK=/staging/Tools/plink_linux_x86_64_20241022/plink

# --- 1) Check the big 6M SNP merged/imputed file ---
echo "=== 1) UZB_for_pruning (5.9M SNPs - post-imputation merged) ==="
BIG="/staging/ALSU-analysis/winter2025/PLINK_301125_0312/step7_ld_pruning/UZB_for_pruning.bim"
echo "--- Search rs5822325 ---"
awk '$2=="rs5822325"' "$BIG"
echo "--- Search position 80107730 ---"
awk '($1==17||$1=="17") && $4==80107730' "$BIG"
echo "--- Search 17:80107730 format ---"
grep -w "17:80107730" "$BIG" | head -3
echo "--- Nearby (80105000-80110000) ---"
awk '($1==17||$1=="17") && $4>=80105000 && $4<=80110000' "$BIG"

# --- 2) Michigan imputation results ---
echo ""
echo "=== 2) IMPUTATION RESULTS DIR ==="
IMPDIR="/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results"
ls -la "$IMPDIR/" 2>/dev/null | head -30
echo "---"
find "$IMPDIR" -maxdepth 2 -type f 2>/dev/null | head -40
echo "---"
# Check for chr17 specifically
find "$IMPDIR" -maxdepth 2 -name "*chr17*" -o -name "*17*" 2>/dev/null | head -10
# Check any chr17 .bim or VCF
for f in $(find "$IMPDIR" -maxdepth 2 \( -name "*chr17*bim" -o -name "*17*bim" \) -type f 2>/dev/null | head -3); do
  echo "Checking $f..."
  awk '$2=="rs5822325" || $4==80107730' "$f" 2>/dev/null
  awk '($1==17||$1=="17") && $4>=80107000 && $4<=80108500' "$f" 2>/dev/null
done
for f in $(find "$IMPDIR" -maxdepth 2 \( -name "*chr17*vcf*" -o -name "*17*vcf*" \) -type f 2>/dev/null | head -3); do
  echo "Checking VCF: $f"
  if [ -f "${f}.tbi" ]; then
    bcftools view "$f" -r 17:80107730-80107730 -H 2>/dev/null | head -3
  fi
done

# --- 3) Old chr17.vcf.gz (pre-imputation per-chr) ---
echo ""
echo "=== 3) PRE-IMPUTATION chr17.vcf.gz ==="
CHR17_VCF="/staging/ALSU-analysis/winter2025/PLINK_301125_0312/old/chr17.vcf.gz"
if [ -f "$CHR17_VCF" ]; then
  echo "Exists: $(ls -lh "$CHR17_VCF" | awk '{print $5}')"
  if [ -f "${CHR17_VCF}.tbi" ]; then
    echo "Indexed, querying 80107730..."
    bcftools view "$CHR17_VCF" -r 17:80107730-80107730 -H 2>/dev/null | head -3
    echo "Nearby:"
    bcftools view "$CHR17_VCF" -r 17:80105000-80110000 -H 2>/dev/null | cut -f1-5 | head -15
  else
    echo "No index, trying zgrep..."
    zgrep "80107730" "$CHR17_VCF" 2>/dev/null | head -3
  fi
fi

# --- 4) Check michigan_ready input chr17 ---
echo ""
echo "=== 4) MICHIGAN-READY INPUT ==="
find /staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr -maxdepth 1 -type f 2>/dev/null | head -30
echo "---"
MR_CHR17=$(find /staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr -maxdepth 1 -name "*chr17*" -o -name "*17*" 2>/dev/null | head -3)
echo "Chr17 files: $MR_CHR17"
for f in $MR_CHR17; do
  echo "Checking: $f"
  if [[ "$f" == *.vcf.gz ]]; then
    if [ -f "${f}.tbi" ]; then
      bcftools view "$f" -r 17:80107730-80107730 -H 2>/dev/null | head -3
    else
      zgrep -m3 "80107730" "$f" 2>/dev/null | head -3
    fi
  fi
done

# --- 5) ConvSK_final_clean.bim check ---
echo ""
echo "=== 5) POST-QC ARRAY (ConvSK_final_clean) ==="
FCLEAN="/staging/ALSU-analysis/winter2025/3_post-imputation/ConvSK_final_clean.bim"
echo "Total SNPs: $(wc -l < "$FCLEAN")"
echo "--- Search rs5822325 ---"
awk '$2=="rs5822325"' "$FCLEAN"
echo "--- Nearby ---"
awk '($1==17||$1=="17") && $4>=80105000 && $4<=80110000' "$FCLEAN"

# --- 6) LD from 1000G using nearby 1000G variants as proxies ---
echo ""
echo "=== 6) LD CHECK: Closest 1000G variants to 80107730 ==="
# From previous run we know 1000G has variants at 80107725 and 80107759 (5bp and 29bp away)
# Let's compute LD between these close 1000G variants and the array probes
VCF17="/staging/ALSU-analysis/Fst_analysis/1000G_data/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
echo "Extracting ±200kb from 1000G..."
bcftools view "$VCF17" -r 17:79907730-80307730 -Oz -o /tmp/r17.vcf.gz 2>/dev/null
bcftools index -t /tmp/r17.vcf.gz 2>/dev/null
"$PLINK" --vcf /tmp/r17.vcf.gz --make-bed --out /tmp/r17 --memory 2000 2>/dev/null

# Get the closest variant IDs to 80107730
echo "Variants at 80107700-80107770:"
awk '$4>=80107700 && $4<=80107770' /tmp/r17.bim

# Use the 80107725 variant (5bp away, likely same LD block) as proxy
PROXY=$(awk '$4==80107725 {print $2; exit}' /tmp/r17.bim 2>/dev/null)
if [ -z "$PROXY" ]; then
  PROXY=$(awk '$4>=80107700 && $4<=80107770 {print $2; exit}' /tmp/r17.bim 2>/dev/null)
fi
echo "Using proxy: $PROXY"

if [ -n "$PROXY" ]; then
  "$PLINK" --bfile /tmp/r17 --ld-snp "$PROXY" \
    --ld-window-kb 250 --ld-window 99999 --ld-window-r2 0 \
    --r2 --out /tmp/ld_proxy --memory 2000 2>/dev/null
  
  if [ -f /tmp/ld_proxy.ld ]; then
    echo ""
    echo "--- Top 20 LD partners of proxy (r2 > 0.3) ---"
    awk 'NR>1 && $7>=0.3 {printf "%-20s pos=%-10s r2=%.4f\n", $6, $5, $7}' /tmp/ld_proxy.ld | sort -t= -k3 -rn | head -20
    
    echo ""
    echo "--- Cross-ref with GSA array probes ---"
    RAW="/staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_120625_0107/ConvSK.bim"
    awk '($1==17||$1=="17") && $4>=79907730 && $4<=80307730 {print $2}' "$RAW" > /tmp/gsa_probes.txt
    echo "GSA probes in region: $(wc -l < /tmp/gsa_probes.txt)"
    
    echo ""
    echo "--- LD between proxy and GSA probes ---"
    awk 'NR==FNR {p[$1]=1;next} FNR>1 && ($6 in p) {printf "GSA: %-30s pos=%-10s r2=%.4f\n",$6,$5,$7}' \
      /tmp/gsa_probes.txt /tmp/ld_proxy.ld | sort -t= -k3 -rn | head -20
    awk 'NR==FNR {p[$1]=1;next} FNR>1 && ($3 in p) {printf "GSA: %-30s pos=%-10s r2=%.4f\n",$3,$2,$7}' \
      /tmp/gsa_probes.txt /tmp/ld_proxy.ld | sort -t= -k3 -rn | head -20
  else
    echo "LD computation failed"
    tail -5 /tmp/ld_proxy.log 2>/dev/null
  fi
fi

rm -f /tmp/r17.* /tmp/ld_proxy.* /tmp/gsa_probes.txt 2>/dev/null

echo ""
echo "=== DONE ==="
