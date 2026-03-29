#!/bin/bash
UNZ=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz
FC=$UNZ/filtered_clean

echo "=== 1. Michigan info.gz structure (chr22 as sample) ==="
zcat $UNZ/chr22.info.gz 2>/dev/null | head -1
echo ""
echo "First 10 data lines:"
zcat $UNZ/chr22.info.gz 2>/dev/null | head -11 | tail -10

echo ""
echo "=== 2. Count TYPED vs IMPUTED in chr22.info.gz ==="
zcat $UNZ/chr22.info.gz 2>/dev/null | awk -F'\t' '
  NR==1 {for(i=1;i<=NF;i++) if($i=="Genotyped") gcol=i; next}
  {t[$gcol]++}
  END {for(k in t) print k": "t[k]}
'

echo ""
echo "=== 3. Count ALL typed variants across all chromosomes ==="
total_typed=0
total_imputed=0
for chr in $(seq 1 22); do
  counts=$(zcat $UNZ/chr${chr}.info.gz 2>/dev/null | awk -F'\t' '
    NR==1 {for(i=1;i<=NF;i++) if($i=="Genotyped") gcol=i; next}
    {t[$gcol]++}
    END {
      typed = t["Genotyped"]+0
      imputed = t["Imputed"]+0
      print typed, imputed
    }
  ')
  typed=$(echo "$counts" | awk '{print $1}')
  imputed=$(echo "$counts" | awk '{print $2}')
  total_typed=$((total_typed + typed))
  total_imputed=$((total_imputed + imputed))
  echo "chr$chr: Typed=$typed, Imputed=$imputed"
done
echo ""
echo "TOTAL: Typed=$total_typed, Imputed=$total_imputed"
echo "Grand total: $((total_typed + total_imputed))"

echo ""
echo "=== 4. Check if typed variants have R2 in info file ==="
echo "Typed variants R2 values (chr22 sample):"
zcat $UNZ/chr22.info.gz 2>/dev/null | awk -F'\t' '
  NR==1 {for(i=1;i<=NF;i++) {if($i=="Genotyped") gcol=i; if($i=="Rsq") rcol=i}; next}
  $gcol=="Genotyped" {print $rcol}
' | head -10

echo ""
echo "=== 5. Typed variants that pass R2>=0.8 and MAF>=0.001 (chr22) ==="
zcat $UNZ/chr22.info.gz 2>/dev/null | awk -F'\t' '
  NR==1 {for(i=1;i<=NF;i++) {if($i=="Genotyped") gcol=i; if($i=="Rsq") rcol=i; if($i=="MAF") mcol=i}; next}
  $gcol=="Genotyped" && $rcol+0>=0.8 && $mcol+0>=0.001
' | wc -l
echo "(total typed in chr22:"
zcat $UNZ/chr22.info.gz 2>/dev/null | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Genotyped") gcol=i; next} $gcol=="Genotyped"' | wc -l
echo ")"

echo ""
echo "=== 6. Check bash_history ==="
cat ~/.bash_history 2>/dev/null

echo ""
echo "=== 7. Check for scripts in filtered_clean ==="
ls -la $FC/*.sh $FC/*.py $FC/*.log 2>/dev/null | head -20

echo ""
echo "=== 8. What are chr10_positions.bed and chr10_ref.log? ==="
head -5 $FC/chr10_positions.bed 2>/dev/null
echo "---"
cat $FC/chr10_ref.log 2>/dev/null

echo ""
echo "=== 9. Check the UZB_all.HQ_imputed file for TYPED column values ==="
head -6 $UNZ/UZB_all.HQ_imputed.R2ge0p8.MAFge0p001.tsv | awk -F'\t' '{print "IMPUTED="$10, "TYPED="$11}'

echo ""
echo "=== 10. Count unique TYPED values in HQ file (first 1000 lines) ==="
awk -F'\t' 'NR>1 && NR<=1001 {print $11}' $UNZ/UZB_all.HQ_imputed.R2ge0p8.MAFge0p001.tsv | sort | uniq -c

echo ""
echo "=== 11. How many TYPED variants are in the HQ file? ==="
awk -F'\t' 'NR>1 && $10==1 && $11!="." {c++} END {print "Typed+Imputed: "c+0}' $UNZ/UZB_all.HQ_imputed.R2ge0p8.MAFge0p001.tsv
awk -F'\t' 'NR>1 && $11=="." {c++} END {print "Imputed only: "c+0}' $UNZ/UZB_all.HQ_imputed.R2ge0p8.MAFge0p001.tsv

echo ""
echo "=== DONE ==="
