#!/bin/bash
# Forensic investigation of 10 dot-ID variants from the BIM file.
# For each: look up its quality in the dose VCF to prove whether it
# should or should not have passed the HQ filter (R²≥0.8, MAF≥0.001, IMPUTED==1).

UNZ=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz
BIM=$UNZ/step7_ld_pruning/UZB_imputed_HQ_clean.bim
HQ=$UNZ/HQ_variant_ids.txt
HQTSV=$UNZ/UZB_all.HQ_imputed.R2ge0p8.MAFge0p001.tsv

echo "=========================================="
echo "DOT-ID VARIANT FORENSIC ANALYSIS"
echo "=========================================="

# --- Step A: Pick 10 dot-ID variants from BIM, spread across chromosomes ---
echo ""
echo "=== A. Selecting 10 dot-ID variants from BIM ==="
# Get 10 examples: 1 from chr1, 1 from chr5, 2 from chr10, 2 from chr15, 2 from chr19, 2 from chr22
dotfile=/tmp/dot_examples.tsv
awk '$2=="." {print $1, $4, $5, $6}' $BIM | awk '
  $1==1  && c1<1  {c1++; print}
  $1==5  && c5<1  {c5++; print}
  $1==10 && c10<2 {c10++; print}
  $1==15 && c15<2 {c15++; print}
  $1==19 && c19<2 {c19++; print}
  $1==22 && c22<2 {c22++; print}
' > $dotfile

echo "Selected variants (CHR POS A1 A2):"
cat $dotfile
echo ""

# --- Step B: For each variant, look it up in the dose VCF ---
echo "=== B. Looking up each variant in dose VCF ==="
echo ""

while read chr pos a1 a2; do
  echo "--- chr${chr}:${pos} (BIM alleles: $a1 / $a2) ---"
  
  # Search in dose VCF for this exact position
  line=$(zcat $UNZ/chr${chr}.dose.vcf.gz 2>/dev/null | awk -v p="$pos" -F'\t' '$2==p && $1!~/^#/ {print; exit}')
  
  if [ -z "$line" ]; then
    echo "  NOT FOUND in dose VCF"
  else
    vcf_id=$(echo "$line" | cut -f3)
    vcf_ref=$(echo "$line" | cut -f4)
    vcf_alt=$(echo "$line" | cut -f5)
    vcf_info=$(echo "$line" | cut -f8)
    
    echo "  VCF ID:  $vcf_id"
    echo "  REF/ALT: $vcf_ref / $vcf_alt"
    
    # Parse INFO field
    typed="NO"
    imputed="NO"
    r2="N/A"
    maf="N/A"
    echo "$vcf_info" | grep -q 'TYPED' && typed="YES"
    echo "$vcf_info" | grep -q 'IMPUTED' && imputed="YES"
    r2=$(echo "$vcf_info" | grep -oP 'R2=[0-9.e+-]+' | cut -d= -f2)
    maf=$(echo "$vcf_info" | grep -oP 'MAF=[0-9.e+-]+' | cut -d= -f2)
    
    echo "  TYPED:   $typed"
    echo "  IMPUTED: $imputed"
    echo "  R2:      $r2"
    echo "  MAF:     $maf"
    
    # Check: does this pass HQ criteria?
    pass="FAIL"
    if [ "$imputed" = "YES" ]; then
      r2_ok=$(awk "BEGIN {print ($r2 >= 0.8) ? 1 : 0}")
      maf_ok=$(awk "BEGIN {print ($maf >= 0.001) ? 1 : 0}")
      if [ "$r2_ok" = "1" ] && [ "$maf_ok" = "1" ]; then
        pass="PASS"
      fi
    fi
    echo "  HQ filter (IMPUTED & R2>=0.8 & MAF>=0.001): $pass"
    
    # Check: is this position in the HQ position file?
    in_hq_pos="NO"
    if [ -f "$UNZ/UZB_all.HQ.query.chrposrefalt.tsv" ]; then
      grep -q "^chr${chr}	${pos}	" $UNZ/UZB_all.HQ.query.chrposrefalt.tsv && in_hq_pos="YES"
    fi
    echo "  In HQ chrposrefalt file: $in_hq_pos"
  fi
  echo ""
done < $dotfile

# --- Step C: Also check 3 rsID variants that ARE in HQ as controls ---
echo "=== C. Control: 3 HQ rsID variants (should all PASS) ==="
# Pick 3 rsIDs from the HQ file
control_ids=$(grep -v '^\.' $HQ | head -3)
for rid in $control_ids; do
  # Find it in BIM
  bim_line=$(awk -v id="$rid" '$2==id {print $1, $4; exit}' $BIM)
  chr=$(echo "$bim_line" | awk '{print $1}')
  pos=$(echo "$bim_line" | awk '{print $2}')
  
  if [ -n "$chr" ] && [ -n "$pos" ]; then
    echo "--- $rid (chr${chr}:${pos}) ---"
    line=$(zcat $UNZ/chr${chr}.dose.vcf.gz 2>/dev/null | awk -v p="$pos" -F'\t' '$2==p && $1!~/^#/ {print; exit}')
    vcf_info=$(echo "$line" | cut -f8)
    r2=$(echo "$vcf_info" | grep -oP 'R2=[0-9.e+-]+' | cut -d= -f2)
    maf=$(echo "$vcf_info" | grep -oP 'MAF=[0-9.e+-]+' | cut -d= -f2)
    typed="NO"; imputed="NO"
    echo "$vcf_info" | grep -q 'TYPED' && typed="YES"
    echo "$vcf_info" | grep -q 'IMPUTED' && imputed="YES"
    echo "  TYPED=$typed IMPUTED=$imputed R2=$r2 MAF=$maf"
  fi
done

echo ""
echo "=== D. Summary counts ==="
echo "Total dot-ID variants in BIM:"
awk '$2=="."' $BIM | wc -l

echo "Total dots in HQ_variant_ids.txt:"
grep -c '^\.$' $HQ

echo "Dot-ID variants in BIM per chr (sample):"
awk '$2=="." {c[$1]++} END {for(k in c) print "chr"k": "c[k]}' $BIM | sort -t: -k1.4 -n

echo ""
echo "=== DONE ==="
