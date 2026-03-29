#!/bin/bash
echo "=== V2 FILE LAYOUT EXPLORATION ==="

echo "--- ~/v2/ directories (2 levels) ---"
find ~/v2/ -maxdepth 2 -type d 2>/dev/null
echo ""

echo "--- All .bim files under ~/v2/ (SNP counts) ---"
for f in $(find ~/v2/ -name "*.bim" 2>/dev/null); do
    echo "$(wc -l < "$f") $f"
done
echo ""

echo "--- All .fam files under ~/v2/ (sample counts) ---"
for f in $(find ~/v2/ -name "*.fam" 2>/dev/null); do
    echo "$(wc -l < "$f") $f"
done
echo ""

echo "--- Looking for non-LD-pruned V2 UZB data ---"
find ~/v2/ -name "UZB_v2*" 2>/dev/null | grep -v "ld_pruned" | head -20
echo ""

echo "--- Looking for global merged (pre-LD-pruning) ---"
find ~/v2/ -name "global*" 2>/dev/null | head -30
echo ""

echo "--- Looking for 1000G reference data ---"
ls -d ~/1000G* /data/1000G* /ref/* ~/ref* 2>/dev/null | head -10
find ~/v2/ -name "*1000*" -o -name "*ref_*" -o -name "*1kg*" 2>/dev/null | head -20
echo ""

echo "--- V1 PBS/global data ---"
find ~ -maxdepth 2 -name "global_for_admixture*" -o -name "pbs*" 2>/dev/null | head -20
find ~ -maxdepth 2 -name "*.fst" 2>/dev/null | head -10
echo ""

echo "--- Scripts from V1 PBS run ---"
find ~ -maxdepth 3 -name "02_calculate_pbs*" -o -name "01_extract*" 2>/dev/null | head -10
echo ""

echo "--- Tool versions ---"
which plink plink2 python3 2>/dev/null
plink --version 2>/dev/null | head -1
python3 --version 2>/dev/null
echo ""

echo "--- Disk space ---"
df -h ~ | tail -1
echo ""

echo "=== DONE ==="
