#!/bin/bash
cd /staging/ALSU-analysis/spring2026/admixture
echo "===GLOBAL_LOGLIK==="
for k in 2 3 4 5 6 7 8; do
  echo -n "K=$k: "
  grep 'Loglikelihood:' "admixture_K${k}.log" | tail -1
done

echo "===UZB_LOGLIK==="
for k in 2 3 4 5 6 7 8; do
  echo -n "K=$k: "
  grep 'Loglikelihood:' "uzb_admixture_K${k}.log" 2>/dev/null | tail -1 || echo "N/A"
done

echo "===UZB_COMPLETE==="
for k in 4 5 6 7 8; do
  echo -n "K=$k: "
  if [ -f "uzb_for_admixture.${k}.Q" ]; then
    echo "DONE ($(wc -l < uzb_for_admixture.${k}.Q) lines)"
  else
    echo "NOT DONE"
  fi
done

echo "===UZB_CV==="
for k in 2 3 4 5 6 7 8; do
  echo -n "K=$k: "
  grep 'CV error' "uzb_admixture_K${k}.log" 2>/dev/null || echo "N/A"
done

echo "===RUNNING==="
ps aux | grep admixture | grep -v grep
