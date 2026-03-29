#!/usr/bin/env bash
set -euo pipefail

echo '=== HOST ==='
hostname

for root in /data /staging; do
  if [ -d "$root" ]; then
    echo "=== ROOT: $root ==="
    echo '-- evanno/admixture/global_v2_admix candidate files (first 500) --'
    find "$root" -maxdepth 8 -type f \( -iname '*evanno*' -o -iname '*admixture*' -o -iname '*global_v2_admix*' \) 2>/dev/null | sort | sed -n '1,500p'
    echo '-- replicate-hint files (first 500) --'
    find "$root" -maxdepth 8 -type f \( -iname '*evanno*' -o -iname '*admixture*' -o -iname '*global_v2_admix*' \) 2>/dev/null | grep -Ei 'rep|run|seed|K[2-8]_[0-9]+|K[2-8]\.[0-9]+' | sort | sed -n '1,500p' || true
  fi
done

echo '=== DONE ==='
