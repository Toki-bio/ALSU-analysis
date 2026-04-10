#!/bin/bash
cd /staging/ALSU-analysis/winter2025/PLINK_301125_0312

# 1. Sex check
echo "=== SEX CHECK ==="
plink --bfile ConvSK_raw --check-sex --out /tmp/alsu_sexcheck 2>&1 | tail -5
cat /tmp/alsu_sexcheck.sexcheck

# 2. Heterozygosity
echo "=== HET CHECK ==="
plink --bfile ConvSK_raw --het --out /tmp/alsu_het 2>&1 | tail -5
cat /tmp/alsu_het.het

echo "=== ALL DONE ==="
