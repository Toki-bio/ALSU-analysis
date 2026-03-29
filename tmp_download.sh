#!/bin/bash
echo "=== UZB K6-8 Q BASE64 ==="
for k in 6 7 8; do
  echo "--- UZB_K=$k ---"
  base64 ~/v2/admixture/UZB_v2_admix.${k}.Q 2>/dev/null
done
echo "=== GLOBAL Q BASE64 ==="
for k in 2 3 4 5 6 7 8; do
  echo "--- GLOBAL_K=$k ---"
  base64 ~/v2/global_admixture/global_v2_admix.${k}.Q 2>/dev/null
done
echo "=== GLOBAL FAM ==="
base64 ~/v2/global_admixture/global_v2_admix.fam 2>/dev/null
echo "=== sNMF CE ==="
cat ~/v2/snmf_uzbek/cross_entropy_summary.csv
echo "=== GLOBAL LOG-LIKELIHOODS ==="
for k in 2 3 4 5 6 7 8; do
  f=~/v2/global_admixture/global_v2_admix_K${k}.log
  if [ -f "$f" ]; then
    echo -n "K=$k: "
    grep 'Loglikelihood:' "$f" | tail -1
  fi
done
echo "=== END ==="
