#!/bin/bash
# Find the v2 reanalysis script and population-related files
echo "=== FIND_V2_SCRIPTS ==="
find ~/v2/ -name "*.sh" -o -name "*.py" -o -name "*.R" 2>/dev/null
echo "=== FIND_POP ==="
find ~/v2/ -name "*pop*" -o -name "*sample*" -o -name "*1000*" 2>/dev/null | head -20
echo "=== V2_LS ==="
ls -la ~/v2/
echo "=== GLOBAL_LS ==="
ls -la ~/v2/global_admixture/
echo "=== REANALYSIS ==="
cat ~/v2/v2_reanalysis.sh 2>/dev/null | grep -A5 -B2 'global\|1000\|pop\|merge\|ref'
echo "=== GLOBAL_LOG ==="
# Check an ADMIXTURE log for dataset info
head -20 ~/v2/global_admixture/global_v2_admix_K2.log 2>/dev/null
echo "=== END ==="
