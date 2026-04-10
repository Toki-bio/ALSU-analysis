#!/bin/bash
cd /staging/ALSU-analysis/winter2025/PLINK_301125_0312

# Count total lines in .genome
echo "=== GENOME FILE STATS ==="
wc -l ConvSK_mind20_ibd.genome

# Extract ALL pairs with PI_HAT >= 0.15 using awk
echo "=== HIGH_PIHAT_PAIRS ==="
awk 'NR==1{print "FID1","IID1","FID2","IID2","Z0","Z1","Z2","PI_HAT"} NR>1 && $10+0>=0.15{print $1,$2,$3,$4,$7,$8,$9,$10}' ConvSK_mind20_ibd.genome

# Also get the PI_HAT distribution summary
echo "=== PIHAT_DISTRIBUTION ==="
awk 'NR>1{
  p=$10+0
  if(p<0.05) b["0.00-0.05"]++
  else if(p<0.10) b["0.05-0.10"]++
  else if(p<0.15) b["0.10-0.15"]++
  else if(p<0.25) b["0.15-0.25"]++
  else if(p<0.50) b["0.25-0.50"]++
  else if(p<0.75) b["0.50-0.75"]++
  else if(p<0.90) b["0.75-0.90"]++
  else if(p<0.98) b["0.90-0.98"]++
  else b["0.98-1.00"]++
} END{for(k in b) print k, b[k]}' ConvSK_mind20_ibd.genome

echo "=== DONE ==="
