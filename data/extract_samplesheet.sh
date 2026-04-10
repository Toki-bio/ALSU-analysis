#!/bin/bash
echo "=== SAMPLE SHEET HEAD ==="
head -30 "/staging/ALSU-analysis/Conversion/1248_merged_sample_sheet_12022025.csv.csv"

echo "=== SAMPLE SHEET TAIL ==="
tail -20 "/staging/ALSU-analysis/Conversion/1248_merged_sample_sheet_12022025.csv.csv"

echo "=== TOTAL LINES ==="
wc -l "/staging/ALSU-analysis/Conversion/1248_merged_sample_sheet_12022025.csv.csv"

echo "=== CHIP DIR COUNT ==="
ls -d /staging/ALSU-analysis/Conversion/OUTPUT/*/ | wc -l

echo "=== ALL CHIP DIRS ==="
ls -d /staging/ALSU-analysis/Conversion/OUTPUT/*/ | xargs -I{} basename {}

echo "=== IDATS PER CHIP ==="
for d in /staging/ALSU-analysis/Conversion/OUTPUT/*/; do
  chip=$(basename "$d")
  cnt=$(find "$d" -name "*.idat" | wc -l)
  echo "$chip $cnt"
done

echo "=== SAMPLE XML EXAMPLE ==="
head -50 /staging/ALSU-analysis/Conversion/OUTPUT/207585740001/207585740001_R01C01_1_Green.xml 2>/dev/null || echo "no xml"

echo "=== CONVSK DIR ==="
ls -la /staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/ 2>/dev/null | head -20

echo "=== LOOK FOR SAMPLE MAP ==="
find /staging/ALSU-analysis/Conversion/ -maxdepth 1 -type f | head -20

echo "=== FULL SAMPLE SHEET ==="
cat "/staging/ALSU-analysis/Conversion/1248_merged_sample_sheet_12022025.csv.csv"

echo "=== DONE ==="
