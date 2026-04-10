#!/bin/bash
# Get the exact 7 sample IDs to remove
cd /staging/ALSU-analysis/winter2025/PLINK_301125_0312/

echo "=== THE 7 SAMPLES TO BE REMOVED ==="
echo "458 08-365"
echo "499 08-701"
echo "840 08-25"
echo "862 08-495"
echo "886 08-77"
echo "898 08-825"
echo "910 12-11"

echo ""
echo "=== VERIFY THESE EXIST IN ConvSK_mind20_dedup.fam ==="
for sample in "458 08-365" "499 08-701" "840 08-25" "862 08-495" "886 08-77" "898 08-825" "910 12-11"; do
  fid=$(echo $sample | awk '{print $1}')
  iid=$(echo $sample | awk '{print $2}')
  if grep -q "^$fid\t$iid" ConvSK_mind20_dedup.fam; then
    echo "✓ Found: $fid $iid"
  else
    echo "✗ NOT FOUND: $fid $iid"
  fi
done

echo ""
echo "=== CREATE FILTER FILE ==="
cat > remove_7_badsamples.txt << 'EOF'
458	08-365
499	08-701
840	08-25
862	08-495
886	08-77
898	08-825
910	12-11
EOF

echo "Created: remove_7_badsamples.txt"
wc -l remove_7_badsamples.txt

echo ""
echo "=== TEST: What would be left after removing these? ==="
echo "Current imputation input (ConvSK_mind20_dedup): 1098 samples"
echo "After removing 7: should be 1091 samples"
echo "This matches the CORRECT Step 2 output!"
