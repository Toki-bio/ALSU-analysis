#!/bin/bash
# parallel download of all Kazakh sample VCFs (resumable)
set -euo pipefail
cd /staging/ALSU-analysis/admixture_analysis/reference_panels/kazakh_gsa_prjeb89820/vcf_all || exit 1
LOG=/staging/tmp/dl_kaz2.log
exec >>"$LOG" 2>&1
echo "start $(date) in $(pwd)"
BASE=https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB89820
curl -sS -m 60 "$BASE/" | grep -o 'href="sample_[A-Za-z0-9_]*\.snv\.vcf\.gz"' | sed 's/href="//; s/"//' | sort -u > files_all.txt
wc -l files_all.txt
dl() { f=$1; [ -s "$f" ] && return 0; curl -sS -fL --retry 4 --retry-delay 5 -m 600 -o "$f.p2" "https://ftp.ebi.ac.uk/pub/databases/eva/PRJEB89820/$f" && mv "$f.p2" "$f"; }
export -f dl
cat files_all.txt | xargs -P 4 -I{} bash -c 'dl {}'
echo "done $(date)"; ls *.snv.vcf.gz | wc -l; du -sh .
