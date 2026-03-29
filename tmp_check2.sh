#!/bin/bash
echo "=== GLOBAL CV ERRORS ==="
for k in 2 3 4 5 6 7 8; do
  grep -i 'CV error' ~/v2/global_admixture/global_v2_admix_K${k}.log 2>/dev/null
done
echo "=== GLOBAL LOG-LIKELIHOODS ==="
for k in 2 3 4 5 6 7 8; do
  f=~/v2/global_admixture/global_v2_admix_K${k}.log
  if [ -f "$f" ]; then
    echo -n "K=$k: "
    grep 'Loglikelihood:' "$f" | tail -1
  fi
done
echo "=== UZB CV ERRORS ==="
for k in 2 3 4 5 6 7 8; do
  grep -i 'CV error' ~/v2/admixture/admixture_K${k}.log 2>/dev/null
done
echo "=== UZB LOG-LIKELIHOODS ==="
for k in 2 3 4 5 6 7 8; do
  f=~/v2/admixture/admixture_K${k}.log
  if [ -f "$f" ]; then
    echo -n "K=$k: "
    grep 'Loglikelihood:' "$f" | tail -1
  fi
done
echo "=== GLOBAL FAM COUNT ==="
wc -l ~/v2/global_admixture/global_v2_admix.fam
echo "=== GLOBAL BIM COUNT ==="
wc -l ~/v2/global_admixture/global_v2_admix.bim
echo "DONE"
