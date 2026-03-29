#!/bin/bash
# Investigation: Why does the final dataset have 10,846,569 variants instead of 10,009,531?
# Look for the ALSU analysis data and understand the filtering logic.

set -e

echo "=== 1. FIND ALSU-RELATED DIRECTORIES ==="
find /home /data /staging /mnt -maxdepth 3 -type d -iname '*alsu*' -o -iname '*uzbek*' -o -iname '*imputation*' 2>/dev/null | head -30

echo ""
echo "=== 2. CHECK HOME DIRECTORY ==="
ls -la ~/ 2>/dev/null | head -20

echo ""
echo "=== 3. FIND HQ VARIANT FILES ==="
find /home /data /staging /mnt -maxdepth 5 -type f \( -name '*HQ_variant*' -o -name '*HQ_imputed*' -o -name '*R2ge*' \) 2>/dev/null | head -20

echo ""
echo "=== 4. FIND PLINK BED/FAM FILES ==="
find /home /data /staging /mnt -maxdepth 5 -type f -name '*.fam' 2>/dev/null | head -30

echo ""
echo "=== 5. FIND VCF FILTERING SCRIPTS ==="
find /home /data /staging /mnt -maxdepth 5 -type f \( -name '*.sh' -o -name '*.py' \) 2>/dev/null | grep -i 'alsu\|uzbek\|imput\|filter\|clean\|step6\|post_imp' | head -20

echo ""
echo "=== 6. FIND LARGE TSV/TXT FILES (frequency tables) ==="
find /home /data /staging /mnt -maxdepth 5 -type f \( -name '*UZB*' -o -name '*uzb*' \) 2>/dev/null | head -30

echo ""
echo "=== 7. CHECK /staging IF EXISTS ==="
ls -la /staging/ 2>/dev/null | head -20
ls -la /staging/ALSU-analysis/ 2>/dev/null | head -20

echo ""
echo "=== 8. FIND ANY LOGS ABOUT VARIANT FILTERING ==="
find /home /data /staging /mnt -maxdepth 4 -type f -name '*.log' 2>/dev/null | grep -i 'alsu\|uzbek\|imput\|filter\|step6' | head -20

echo ""
echo "=== DONE ==="
