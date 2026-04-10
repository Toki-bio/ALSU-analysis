#!/bin/bash
cd /staging/ALSU-analysis/winter2025/PLINK_301125_0312/

echo "=== VERIFY THE 7 MISSED SAMPLES ARE IN STEP 1 OUTPUT ==="
cat > /tmp/7samples.txt << 'EOF'
458	08-365
499	08-701
840	08-25
862	08-495
886	08-77
898	08-825
910	12-11
EOF

echo "Looking for these samples in ConvSK_mind20.fam (Step 1):"
while IFS=$'\t' read fid iid; do
  if grep -q "^$fid[[:space:]]$iid" ConvSK_mind20.fam; then
    echo "✓ FOUND: $fid $iid (in Step 1 output - this is correct)"
  else
    echo "✗ NOT FOUND: $fid $iid"
  fi
done < /tmp/7samples.txt

echo ""
echo "=== STRATEGY FOR FILTERING ==="
echo ""
echo "Option 1: Quick fix for downstream (NO IMPUTATION RE-RUN)"
echo "  Step 1: Keep as-is (1155 samples with error)"
echo "  Step 2: Keep as-is (1098 samples)"
echo "  Step 6-15: Filter OUT the 7 bad samples from imputation results"
echo "Results: Final analyses use 1091 samples, matching correct Step 2 output"
echo ""
echo "Option 2: Proper fix (REQUIRES IMPUTATION RE-RUN)"
echo "  Re-run Step 2 with corrected Step 1"
echo "  Re-run imputation with 1,091 samples"
echo "  Re-run Step 6-15"
echo ""
echo "recommendation: Use Option 1 (much faster)"
echo "Just remove these 7 from the final imputation results!"
