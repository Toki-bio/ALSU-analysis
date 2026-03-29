#!/bin/bash
echo "=== ADMIXTURE ==="
tail -5 ~/v2/reanalysis.log 2>/dev/null
echo "=== PROCESSES ==="
ps aux | grep -E 'admixture|Rscript' | grep -v grep
echo "=== Q FILES ==="
ls -lrt ~/v2/admixture/*.Q 2>/dev/null
echo "=== CV ERRORS ==="
for k in 2 3 4 5 6 7 8; do
  grep -i 'CV error' ~/v2/admixture/admixture_K${k}.log 2>/dev/null
done
echo "=== LOG-LIKELIHOODS ==="
for k in 2 3 4 5 6 7 8; do
  f=~/v2/admixture/admixture_K${k}.log
  if [ -f "$f" ]; then
    echo -n "K=$k: "
    grep 'Loglikelihood:' "$f" | tail -1
  fi
done
echo "=== sNMF CE ==="
cat ~/v2/snmf_uzbek/cross_entropy_summary.csv
echo "=== PHASE C ==="
ls ~/v2/global_admixture/ 2>/dev/null
echo "DONE"
