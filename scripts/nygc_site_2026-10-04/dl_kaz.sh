#!/bin/bash
# Download the Kazakh GSA genotype VCFs (EVA PRJEB89820, CC BY-NC-ND) into /staging; one file per sample.
set -euo pipefail
cd /staging/ALSU-analysis/admixture_analysis/reference_panels/kazakh_gsa_prjeb89820 || exit 1
mkdir -p vcf_all
cd vcf_all || exit 1
LOG=/staging/tmp/dl_kaz.log
exec >>"$LOG" 2>&1
echo "start $(date) in $(pwd)"
BASE=https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB89820
curl -sS -m 60 "$BASE/" | grep -o 'href="sample_[0-9]*\.snv\.vcf\.gz"' | sed 's/href="//; s/"//' | sort -u > files.txt
wc -l files.txt
n=0
while read -r f; do
  if [ -s "$f" ]; then n=$((n+1)); continue; fi
  curl -sS -fL --retry 4 --retry-delay 5 -m 300 -o "$f.part" "$BASE/$f"
  mv "$f.part" "$f"
  n=$((n+1))
  [ $((n % 20)) -eq 0 ] && echo "$n files $(date)"
done < files.txt
echo "done $n files $(date)"
ls | wc -l
du -sh .
