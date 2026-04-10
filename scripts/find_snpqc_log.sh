#!/bin/bash
# Find the original PLINK command that created ConvSK_mind20_dedup_snpqc in winter 2025
cd /staging/ALSU-analysis/winter2025/PLINK_301125_0312/

echo "=== PLINK log for snpqc ==="
cat ConvSK_mind20_dedup_snpqc.log 2>/dev/null || echo "No .log file found"

echo ""
echo "=== All .log files mentioning snpqc ==="
grep -l "snpqc" *.log 2>/dev/null || echo "No .log files mention snpqc"

echo ""
echo "=== All .log files in directory ==="
ls -lt *.log 2>/dev/null | head -20

echo ""
echo "=== Search for --geno in all logs ==="
grep -h "\-\-geno" *.log 2>/dev/null | head -20

echo ""
echo "=== Check spring2026 logs too ==="
cd /staging/ALSU-analysis/spring2026/
ls -lt *.log 2>/dev/null | head -20

echo ""
echo "=== spring2026 snpqc log ==="
cat ConvSK_mind20_dedup_snpqc.log 2>/dev/null || echo "No .log file"
