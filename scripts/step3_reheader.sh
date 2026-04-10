#!/bin/bash
# Run on Biotech2024 in /staging/ALSU-analysis/spring2026/
# Creates rename_samples.tsv and applies reheader

# The two non-ASCII IDs found:
#   808_03-25м    (Cyrillic м = em, U+043C)
#   1038_08-176Х-00006  (Cyrillic Х = Ha, U+0425)
#
# Rename to Latin equivalents: м -> m, Х -> X

# Step 1: Create rename file using printf to guarantee real tabs
# Cyrillic м = UTF-8 bytes \xd0\xbc, Cyrillic Х = UTF-8 bytes \xd0\xa5
printf '808_03-25\xd0\xbc\t808_03-25m\n' > rename_samples.tsv
printf '1038_08-176\xd0\xa5-00006\t1038_08-176X-00006\n' >> rename_samples.tsv

# Verify: should see ^I (tab) between old and new names
cat -A rename_samples.tsv

# Step 2: Apply reheader
bcftools reheader -s rename_samples.tsv \
  -o impute_in.final.vcf.gz \
  impute_in.autosomes.clean.nodup.vcf.gz

# Step 3: Index, then count final variants
bcftools index impute_in.final.vcf.gz
bcftools index -n impute_in.final.vcf.gz
