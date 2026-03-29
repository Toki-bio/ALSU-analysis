#!/bin/bash
# Fast forensic: use tabix for indexed lookups instead of zcat+awk
UNZ=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz
BIM=$UNZ/step7_ld_pruning/UZB_imputed_HQ_clean.bim
HQ=$UNZ/HQ_variant_ids.txt

echo "=========================================="
echo "DOT-ID VARIANT FORENSIC ANALYSIS (FAST)"
echo "=========================================="

# Check for tabix index
echo "=== 0. Index check ==="
ls -la $UNZ/chr22.dose.vcf.gz.tbi 2>/dev/null || echo "No tabix index for dose VCF"
ls -la $UNZ/filtered_clean/chr22.clean.vcf.gz.tbi 2>/dev/null || echo "No tabix index for clean VCF"

# Step A: Pick 10 dot-ID variants from BIM spread across chromosomes
echo ""
echo "=== A. Selecting 10 dot-ID variants ==="
awk '$2=="." {print $1, $4, $5, $6}' $BIM | awk '
  $1==1  && c1<1  {c1++; print}
  $1==5  && c5<1  {c5++; print}
  $1==10 && c10<2 {c10++; print}
  $1==15 && c15<2 {c15++; print}
  $1==19 && c19<2 {c19++; print}
  $1==22 && c22<2 {c22++; print}
' > /tmp/dot_examples.tsv

echo "Selected variants:"
cat /tmp/dot_examples.tsv

# Step B: For each variant, use bcftools query on the dose VCF (indexed)
echo ""
echo "=== B. Querying dose VCF for each variant ==="
while read chr pos a1 a2; do
  echo ""
  echo "--- chr${chr}:${pos} (BIM: $a1/$a2) ---"
  
  # Use bcftools to query specific position - much faster than zcat
  result=$(bcftools query -r "chr${chr}:${pos}-${pos}" \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/MAF\t%INFO/R2\n' \
    $UNZ/chr${chr}.dose.vcf.gz 2>/dev/null)
  
  if [ -z "$result" ]; then
    echo "  NOT FOUND via bcftools query"
    # Fallback: check first few hundred lines after header
    echo "  (Trying grep on position...)"
    result_grep=$(zcat $UNZ/chr${chr}.dose.vcf.gz 2>/dev/null | grep -m1 "^chr${chr}	${pos}	")
    if [ -n "$result_grep" ]; then
      echo "  Found via grep:"
      echo "  $(echo "$result_grep" | cut -f1-8)"
    fi
  else
    echo "$result" | while IFS=$'\t' read rc rp rid rref ralt raf rmaf rr2; do
      echo "  ID:      $rid"
      echo "  REF/ALT: $rref / $ralt"
      echo "  AF:      $raf"
      echo "  MAF:     $rmaf"
      echo "  R2:      $rr2"
    done
    
    # Also get TYPED/IMPUTED flags from INFO
    info_line=$(bcftools query -r "chr${chr}:${pos}-${pos}" \
      -f '%INFO\n' $UNZ/chr${chr}.dose.vcf.gz 2>/dev/null | head -1)
    
    typed="NO"; imputed="NO"
    echo "$info_line" | grep -q 'TYPED' && typed="YES"
    echo "$info_line" | grep -q 'IMPUTED' && imputed="YES"
    echo "  TYPED:   $typed"
    echo "  IMPUTED: $imputed"
    
    # Check HQ criteria
    r2val=$(bcftools query -r "chr${chr}:${pos}-${pos}" -f '%INFO/R2\n' $UNZ/chr${chr}.dose.vcf.gz 2>/dev/null | head -1)
    mafval=$(bcftools query -r "chr${chr}:${pos}-${pos}" -f '%INFO/MAF\n' $UNZ/chr${chr}.dose.vcf.gz 2>/dev/null | head -1)
    
    pass="FAIL"
    if [ "$imputed" = "YES" ]; then
      r2_ok=$(awk "BEGIN {print ($r2val+0 >= 0.8) ? 1 : 0}" 2>/dev/null)
      maf_ok=$(awk "BEGIN {print ($mafval+0 >= 0.001) ? 1 : 0}" 2>/dev/null)
      [ "$r2_ok" = "1" ] && [ "$maf_ok" = "1" ] && pass="PASS"
    fi
    echo "  HQ filter (IMPUTED & R2>=0.8 & MAF>=0.001): $pass"
  fi
  
  # Check position file
  in_pos="NO"
  grep -q "^chr${chr}	${pos}	" $UNZ/UZB_all.HQ.query.chrposrefalt.tsv 2>/dev/null && in_pos="YES"
  echo "  In HQ position file: $in_pos"
  
done < /tmp/dot_examples.tsv

# Step C: 3 control rsID variants from HQ list
echo ""
echo "=== C. Control: 3 HQ rsID variants ==="
grep -v '^\.' $HQ | head -3 | while read rid; do
  bim_line=$(awk -v id="$rid" '$2==id {print $1, $4; exit}' $BIM)
  chr=$(echo "$bim_line" | awk '{print $1}')
  pos=$(echo "$bim_line" | awk '{print $2}')
  
  if [ -n "$chr" ] && [ -n "$pos" ]; then
    result=$(bcftools query -r "chr${chr}:${pos}-${pos}" \
      -f '%ID\t%INFO/R2\t%INFO/MAF\n' \
      $UNZ/chr${chr}.dose.vcf.gz 2>/dev/null | head -1)
    info=$(bcftools query -r "chr${chr}:${pos}-${pos}" -f '%INFO\n' \
      $UNZ/chr${chr}.dose.vcf.gz 2>/dev/null | head -1)
    typed="NO"; imputed="NO"
    echo "$info" | grep -q 'TYPED' && typed="YES"
    echo "$info" | grep -q 'IMPUTED' && imputed="YES"
    echo "$rid (chr${chr}:${pos}): $result  TYPED=$typed IMPUTED=$imputed"
  fi
done

# Step D: Summary
echo ""
echo "=== D. Summary ==="
echo "Dots in BIM: $(awk '$2=="."' $BIM | wc -l)"
echo "Dots in HQ:  $(grep -c '^\.$' $HQ)"
echo ""
echo "Per-chr dot count in BIM:"
awk '$2=="." {c[$1]++} END {for(k in c) printf "chr%s: %d\n", k, c[k]}' $BIM | sort -t: -k1 -V

echo ""
echo "=== DONE ==="
