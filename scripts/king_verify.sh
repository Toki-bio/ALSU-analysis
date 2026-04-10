#!/bin/bash
set -e
cd /staging/ALSU-analysis/winter2025/PLINK_301125_0312/

echo "=== PLINK2 VERSION ==="
plink2 --version 2>&1 | head -1

echo "=== RUNNING KING-ROBUST KINSHIP (duplicate threshold >= 0.354) ==="
plink2 --bfile ConvSK_raw \
  --make-king-table \
  --king-table-filter 0.354 \
  --out /tmp/king_duplicates \
  --threads 4 \
  --memory 4000

echo "=== KING DUPLICATE-LEVEL PAIRS ==="
wc -l /tmp/king_duplicates.kin0
cat /tmp/king_duplicates.kin0

echo "=== ALSO RUNNING AT 2ND-DEGREE THRESHOLD (>= 0.177) TO CATCH EDGE CASES ==="
plink2 --bfile ConvSK_raw \
  --make-king-table \
  --king-table-filter 0.177 \
  --out /tmp/king_2nd_degree \
  --threads 4 \
  --memory 4000

echo "=== KING 2ND-DEGREE PAIRS (count) ==="
wc -l /tmp/king_2nd_degree.kin0
echo "=== HEADER ==="
head -1 /tmp/king_2nd_degree.kin0
echo "=== PAIRS BETWEEN 0.177 AND 0.354 (non-duplicate but related) ==="
awk 'NR>1 && $10 < 0.354 && $10 >= 0.177' /tmp/king_2nd_degree.kin0

echo "=== DONE ==="
