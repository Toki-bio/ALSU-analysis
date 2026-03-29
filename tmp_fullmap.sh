#!/bin/bash
# Build complete sample->superpop mapping for all 3595 samples in global ADMIXTURE
# Output: one line per sample: IID<tab>SUPERPOP
# For ref: use pca pop_mapping. For Uzbeks: assign "UZB"

echo "=== FULL_MAPPING ==="
# First 2548 are ref, remaining 1047 are Uzbek
head -2548 ~/v2/global_admixture/global_v2_admix.fam | awk '{print $2}' | \
  awk 'NR==FNR{pop[$1]=$2;next} {print $1"\t"($1 in pop ? pop[$1] : "UNK")}' ~/v2/pca/pop_mapping.txt -

tail -1047 ~/v2/global_admixture/global_v2_admix.fam | awk '{print $2"\tUZB"}'
echo "=== END ==="
