#!/bin/bash
# Corrected forensic: BIM is in filtered_clean, dose VCFs have no tabix index
# Strategy: pick dots from BIM, look them up in clean VCFs (indexed) and dose VCFs (grep)
UNZ=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz
FC=$UNZ/filtered_clean
BIM=$FC/UZB_imputed_HQ_clean.bim
HQ=$UNZ/HQ_variant_ids.txt
HQPOS=$UNZ/UZB_all.HQ.query.chrposrefalt.tsv

echo "=========================================="
echo "DOT-ID VARIANT FORENSIC ANALYSIS v3"
echo "=========================================="

# Step A: Pick 10 dot-ID variants across chromosomes
echo ""
echo "=== A. Selecting 10 dot-ID variants from BIM ==="
awk '$2=="." {print $1, $4, $5, $6}' $BIM | awk '
  $1==1  && c1<1  {c1++; print}
  $1==5  && c5<1  {c5++; print}
  $1==10 && c10<2 {c10++; print}
  $1==15 && c15<2 {c15++; print}
  $1==19 && c19<2 {c19++; print}
  $1==22 && c22<2 {c22++; print}
' > /tmp/dot_examples.tsv

echo "Selected (CHR POS A1 A2):"
cat /tmp/dot_examples.tsv
echo ""

# Step B: For each, query the clean VCF (indexed!) and dose VCF (small grep)
echo "=== B. Per-variant analysis ==="
n=0
while read chr pos a1 a2; do
  n=$((n+1))
  echo ""
  echo "--- #$n: chr${chr}:${pos} (BIM alleles: $a1 / $a2) ---"

  # 1) Is it in the clean VCF? (these have tabix indexes)
  clean_line=$(bcftools query -r "chr${chr}:${pos}-${pos}" -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n' $FC/chr${chr}.clean.vcf.gz 2>/dev/null | head -1)
  if [ -n "$clean_line" ]; then
    echo "  [clean VCF] PRESENT"
    echo "  $clean_line" | awk -F'\t' '{print "  ID="$3, "REF="$4, "ALT="$5}'
    # Extract INFO flags
    info=$(echo "$clean_line" | cut -f6)
    typed="NO"; imputed="NO"
    echo "$info" | grep -q 'TYPED' && typed="YES"
    echo "$info" | grep -q 'IMPUTED' && imputed="YES"
    r2=$(echo "$info" | grep -oP 'R2=[0-9.e+-]+' | cut -d= -f2)
    maf=$(echo "$info" | grep -oP 'MAF=[0-9.e+-]+' | cut -d= -f2)
    echo "  TYPED=$typed IMPUTED=$imputed R2=$r2 MAF=$maf"
    
    # HQ filter check
    pass="FAIL"
    if [ "$imputed" = "YES" ] && [ -n "$r2" ] && [ -n "$maf" ]; then
      ok=$(awk "BEGIN {print (($r2+0 >= 0.8) && ($maf+0 >= 0.001)) ? 1 : 0}")
      [ "$ok" = "1" ] && pass="PASS"
    fi
    echo "  Would pass HQ filter (IMPUTED & R2>=0.8 & MAF>=0.001): $pass"
  else
    echo "  [clean VCF] NOT FOUND"
  fi

  # 2) Is this position in the HQ position file?
  in_pos="NO"
  grep -q "^chr${chr}	${pos}	" $HQPOS 2>/dev/null && in_pos="YES"
  echo "  In HQ position file (chrposrefalt): $in_pos"

done < /tmp/dot_examples.tsv

# Step C: Control group - 3 rsID variants from HQ that should pass
echo ""
echo "=== C. Control: 3 HQ rsID variants ==="
grep -v '^\.' $HQ | head -3 | while read rid; do
  bim_line=$(awk -v id="$rid" '$2==id {print $1, $4; exit}' $BIM)
  chr=$(echo "$bim_line" | awk '{print $1}')
  pos=$(echo "$bim_line" | awk '{print $2}')
  if [ -n "$chr" ] && [ -n "$pos" ]; then
    info=$(bcftools query -r "chr${chr}:${pos}-${pos}" -f '%ID\t%INFO\n' $FC/chr${chr}.clean.vcf.gz 2>/dev/null | head -1)
    vid=$(echo "$info" | cut -f1)
    ifield=$(echo "$info" | cut -f2)
    typed="NO"; imputed="NO"
    echo "$ifield" | grep -q 'TYPED' && typed="YES"
    echo "$ifield" | grep -q 'IMPUTED' && imputed="YES"
    r2=$(echo "$ifield" | grep -oP 'R2=[0-9.e+-]+' | cut -d= -f2)
    maf=$(echo "$ifield" | grep -oP 'MAF=[0-9.e+-]+' | cut -d= -f2)
    echo "$rid (chr${chr}:${pos}): ID=$vid R2=$r2 MAF=$maf TYPED=$typed IMPUTED=$imputed"
  else
    echo "$rid: NOT FOUND in BIM"
  fi
done

# Step D: Summary stats
echo ""
echo "=== D. Summary ==="
echo "Dots in BIM: $(awk '$2=="."' $BIM | wc -l)"
echo "Dots in HQ:  $(grep -c '^\.$' $HQ)"
echo "Per-chr dots in BIM:"
awk '$2=="." {c[$1]++} END {for(k in c) printf "  chr%s: %d\n", k, c[k]}' $BIM | sort -t: -k1 -V

echo ""
echo "=== DONE ==="
