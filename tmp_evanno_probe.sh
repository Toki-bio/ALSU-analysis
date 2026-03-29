#!/usr/bin/env bash
set -euo pipefail

echo '=== HOST ==='
hostname

echo '=== DIRS UNDER /home/copilot/v2 ==='
find /home/copilot/v2 -maxdepth 3 -type d | sort

echo '=== GLOBAL ADMIXTURE PRIMARY FILES ==='
if [ -d /home/copilot/v2/global_admixture ]; then
  cd /home/copilot/v2/global_admixture
  pwd
  ls -1 | sort
  echo '--- K log files ---'
  ls -1 global_v2_admix_K*.log 2>/dev/null | sort || true
  echo 'K_LOG_COUNT='$(ls -1 global_v2_admix_K*.log 2>/dev/null | wc -l)
  echo '--- CV lines from K logs ---'
  grep -H 'CV error' global_v2_admix_K*.log 2>/dev/null || true
  echo '--- Log-likelihood lines from K logs ---'
  grep -HiH 'loglikelihood\|log likelihood\|loglik' global_v2_admix_K*.log 2>/dev/null || true
  echo '--- replicate-like filenames in this dir ---'
  ls -1 | grep -Ei 'rep|run|seed|K[2-8]_[0-9]+|K[2-8]\.[0-9]+' || true
fi

echo '=== SEARCH REPLICATE ARTIFACTS UNDER /home/copilot ==='
find /home/copilot -maxdepth 7 -type f \( -iname '*admixture*' -o -iname '*evanno*' -o -iname '*global_v2_admix*' \) | sort

echo '=== FILTER: FILES SUGGESTING REPLICATES ==='
find /home/copilot -maxdepth 7 -type f \( -iname '*admixture*' -o -iname '*evanno*' -o -iname '*global_v2_admix*' \) | grep -Ei 'rep|run|seed|K[2-8]_[0-9]+|K[2-8]\.[0-9]+' | sort || true

echo '=== UZB ADMIXTURE DIR ==='
if [ -d /home/copilot/v2/admixture ]; then
  cd /home/copilot/v2/admixture
  pwd
  ls -1 UZB_v2_admix* 2>/dev/null | sort || true
  echo '--- UZB CV lines ---'
  grep -H 'CV error' UZB_v2_admix*.log 2>/dev/null || true
fi

echo '=== DONE ==='
