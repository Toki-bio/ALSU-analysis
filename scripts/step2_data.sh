#!/bin/bash
set -euo pipefail
cd /staging/ALSU-analysis/spring2026/

echo "=== GENOME_FILE_INFO ==="
wc -l ConvSK_mind20.genome
head -1 ConvSK_mind20.genome

echo "=== ALL_PAIRS ==="
# All pairs (skip header), print FID1 IID1 FID2 IID2 PI_HAT
awk 'NR>1 {print $1, $2, $3, $4, $10}' ConvSK_mind20.genome

echo "=== PAIR_COUNT ==="
awk 'NR>1' ConvSK_mind20.genome | wc -l

echo "=== PIHAT_DISTRIBUTION ==="
# Distribution of PI_HAT values in 0.01 bins from 0.98 to 1.00
awk 'NR>1 {
  p = $10+0
  if (p >= 1.00) bin="1.00"
  else if (p >= 0.99) bin="0.99-1.00"
  else bin="0.98-0.99"
  counts[bin]++
}
END {
  for (b in counts) print b, counts[b]
}' ConvSK_mind20.genome

echo "=== UNIQUE_SAMPLES ==="
# All unique sample IDs involved in high-IBD pairs
awk 'NR>1 {print $2; print $4}' ConvSK_mind20.genome | sort -u

echo "=== UNIQUE_SAMPLE_COUNT ==="
awk 'NR>1 {print $2; print $4}' ConvSK_mind20.genome | sort -u | wc -l

echo "=== SAMPLE_PAIR_COUNTS ==="
# How many pairs each sample is in (to find hubs)
awk 'NR>1 {c[$2]++; c[$4]++} END {for (s in c) print c[s], s}' ConvSK_mind20.genome | sort -rn

echo "=== FAM_SAMPLE_COUNT ==="
wc -l ConvSK_mind20.fam

echo "=== FULL_GENOME_TABLE ==="
# Full .genome file content (all columns)
cat ConvSK_mind20.genome

echo "=== END ==="
