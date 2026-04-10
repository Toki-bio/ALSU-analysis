#!/bin/bash
cd /staging/ALSU-analysis/winter2025/PLINK_301125_0312/

echo "=== CHECK STEP 1 OUTPUT (1155 samples) ==="
for sample in "458 08-365" "499 08-701" "840 08-25" "862 08-495" "886 08-77" "898 08-825" "910 12-11"; do
  fid=$(echo $sample | awk '{print $1}')
  iid=$(echo $sample | awk '{print $2}')
  if grep -q "^$fid\t$iid" ConvSK_mind20.fam; then
    echo "✓ IN STEP 1: $fid $iid"
  else
    echo "✗ NOT IN STEP 1: $fid $iid (removed earlier or doesn't exist)"
  fi
done

echo ""
echo "=== CHECK STEP 2 OUTPUT (1098 samples) ===" 
for sample in "458 08-365" "499 08-701" "840 08-25" "862 08-495" "886 08-77" "898 08-825" "910 12-11"; do
  fid=$(echo $sample | awk '{print $1}')
  iid=$(echo $sample | awk '{print $2}')
  if grep -q "^$fid\t$iid" ConvSK_mind20_dedup.fam; then
    echo "✓ IN STEP 2: $fid $iid"
  else
    echo "✗ REMOVED IN STEP 2: $fid $iid (was in Step 1 but removed via IBD/dedup)"
  fi
done

echo ""
echo "=== WHICH ONES ARE IN STEP 1 BUT NOT IN STEP 2? ==="
echo "These are the ones removed by IBD deduplication (good news!)"
