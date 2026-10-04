#!/bin/bash
set -euo pipefail
cd /staging/tmp || exit 1
SRC=/home/copilot/v2
DST=/staging/ALSU-analysis/spring2026/full_expanded_cohort/alsu_v2
LOG=/staging/tmp/v2_move_2026-10-04.log
exec >>"$LOG" 2>&1
echo "== start $(date)"
[ -d "$SRC" ] && [ ! -L "$SRC" ] || { echo "SRC missing or already a symlink"; exit 1; }
[ ! -e "$DST" ] || { echo "DST exists, abort"; exit 1; }
mkdir "$DST"
rsync -a "$SRC"/ "$DST"/
echo "copy done $(date)"
n1=$(find "$SRC" -type f | wc -l); n2=$(find "$DST" -type f | wc -l)
s1=$(du -sb "$SRC" | cut -f1); s2=$(du -sb "$DST" | cut -f1)
echo "files $n1 $n2  bytes $s1 $s2"
[ "$n1" = "$n2" ] && [ "$s1" = "$s2" ] || { echo "MISMATCH, source kept"; exit 1; }
diffs=$(rsync -anc --itemize-changes "$SRC"/ "$DST"/ | wc -l)
echo "checksum diffs: $diffs"
[ "$diffs" = 0 ] || { echo "CHECKSUM MISMATCH, source kept"; exit 1; }
rm -rf "$SRC"
ln -s "$DST" "$SRC"
echo "== moved, symlink in place $(date)"
df -h /home
