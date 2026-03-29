#!/bin/bash
# Get the population mapping used in V2 global ADMIXTURE
echo "=== POP_MAP ==="
# Check common locations for the population mapping file
for f in ~/v2/global_admixture/pop_mapping.txt ~/v2/global_admixture/populations.txt ~/v2/global_admixture/ref_populations.txt ~/v2/global_admixture/*.pop ~/v2/plink/*.pop; do
  if [ -f "$f" ]; then
    echo "--- FILE: $f ---"
    head -5 "$f"
    echo "..."
    wc -l "$f"
  fi
done
# Also check the v2 reanalysis script for population reference
echo "=== SCRIPT_REFS ==="
grep -rn 'pop\|population\|1000G\|1000g\|ref_sample\|superpop' ~/v2/global_admixture/*.sh ~/v2/*.sh 2>/dev/null | head -30
# Check if there's a .pop file alongside the bed/bim/fam
echo "=== POP_FILE ==="
ls -la ~/v2/global_admixture/*.pop 2>/dev/null
ls -la ~/v2/global_admixture/pop* 2>/dev/null
# Check if the reanalysis script is there
echo "=== SCRIPTS ==="
ls -la ~/v2/*.sh ~/v2/global_admixture/*.sh 2>/dev/null
# Also try to find 1000G population info 
echo "=== 1000G_POP ==="
for f in ~/v2/global_admixture/integrated_call_samples* ~/v2/global_admixture/igsr_samples* ~/v2/global_admixture/sample_info*; do
  if [ -f "$f" ]; then
    echo "--- FILE: $f ---"
    head -3 "$f"
    echo "..."
    wc -l "$f"
  fi
done
# Last resort: check the reanalysis script for how populations were assigned
echo "=== V2_REANALYSIS ==="
cat ~/v2/v2_reanalysis.sh 2>/dev/null | head -100
echo "=== END ==="
