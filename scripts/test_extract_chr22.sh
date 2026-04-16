#!/bin/bash
# Test extraction with chr22 (smallest zip, 1.1GB)
set -euo pipefail
cd /staging/ALSU-analysis/spring2026/imputation

echo "=== DISK BEFORE ==="
df -h /staging/

echo "=== TEST: chr22 extraction ==="
# Password contains & and ; — double-quoted variable is safe
IMPUTE_PASSWORD='d6rFCYsg&q;RZ0'

unzip -P "$IMPUTE_PASSWORD" -o chr_22.zip \
  chr22.dose.vcf.gz \
  chr22.info.gz 2>&1 | tail -5

echo "=== VERIFY ==="
ls -lh chr22.dose.vcf.gz chr22.info.gz 2>/dev/null || echo "EXTRACTION FAILED"

echo "=== SAMPLE COUNT ==="
bcftools query -l chr22.dose.vcf.gz 2>/dev/null | wc -l

echo "=== SAMPLE IDS (first 5) ==="
bcftools query -l chr22.dose.vcf.gz 2>/dev/null | head -5

echo "=== VARIANT COUNT ==="
bcftools index -f chr22.dose.vcf.gz 2>/dev/null
bcftools index -n chr22.dose.vcf.gz 2>/dev/null

echo "=== R2 QUICK CHECK ==="
zcat chr22.info.gz | grep -v '^#' | head -3 | cut -f1-8

echo "=== DISK AFTER ==="
df -h /staging/

echo "=== DONE ==="
