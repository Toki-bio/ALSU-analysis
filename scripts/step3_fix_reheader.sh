#!/bin/bash
# Run on Biotech2024 in /staging/ALSU-analysis/spring2026/
# Finds non-ASCII sample IDs and creates rename_samples.tsv

# Step 1: Find samples with non-ASCII characters
echo "=== Non-ASCII sample IDs ==="
bcftools query -l impute_in.autosomes.clean.nodup.vcf.gz | grep -P '[^\x00-\x7F]'

# Step 2: Show all sample IDs for reference (pipe to less or head)
# bcftools query -l impute_in.autosomes.clean.nodup.vcf.gz | head -20

# Step 3: Once you see the Cyrillic IDs, create rename_samples.tsv manually:
#   Format: old_name<TAB>new_name (one line per sample to rename)
#   Example:
#     Иванов_08-123	Ivanov_08-123
#     Петров_08-456	Petrov_08-456
#
# Then run:
#   bcftools reheader -s rename_samples.tsv \
#     -o impute_in.final.vcf.gz \
#     impute_in.autosomes.clean.nodup.vcf.gz

echo ""
echo "=== Total variant count (before reheader) ==="
bcftools index -n impute_in.autosomes.clean.nodup.vcf.gz
