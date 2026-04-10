#!/bin/bash
# Re-export chr10-22 from PLINK binaries and run fixref
# CRITICAL: fixref stderr must NOT go into the pipe
SRCDIR="/tmp/step3_geno10"
READYDIR="/tmp/michigan_ready"
RENAME="/tmp/chr_rename.tsv"
REF="/staging/Genomes/Human/chr/GRCh38.fa"
PLINK_PREFIX="${SRCDIR}/ConvSK_mind20_dedup_snpqc"
LOG="/tmp/fixref_v2.log"
STATS="/tmp/fixref_stats.log"

exec > "$LOG" 2>&1
echo "START $(date)"

for i in $(seq 10 22); do
  OUT="${READYDIR}/impute_chr${i}.vcf.gz"
  if [ -f "$OUT" ] && [ -f "${OUT}.tbi" ]; then
    n=$(bcftools index -n "$OUT" 2>/dev/null || echo 0)
    if [ "$n" -gt 0 ] 2>/dev/null; then
      echo "chr${i}: SKIP_DONE (${n})"
      continue
    fi
  fi
  rm -f "$OUT" "${OUT}.tbi"

  echo "chr${i}: RE-EXPORT $(date +%H:%M:%S)"
  plink --bfile "$PLINK_PREFIX" --chr "$i" \
    --recode vcf-iid bgz \
    --out "/tmp/reexport_chr${i}"

  if [ ! -f "/tmp/reexport_chr${i}.vcf.gz" ]; then
    echo "chr${i}: PLINK_FAIL"
    continue
  fi

  REEXPORT="/tmp/reexport_chr${i}.vcf.gz"

  echo "chr${i}: RENAME $(date +%H:%M:%S)"
  TMPVCF="/tmp/mw_chr${i}.vcf.gz"
  bcftools annotate --rename-chrs "$RENAME" "$REEXPORT" -Oz -o "$TMPVCF"
  if [ $? -ne 0 ]; then
    echo "chr${i}: RENAME_FAIL"
    rm -f "$TMPVCF" "$REEXPORT"
    continue
  fi
  bcftools index "$TMPVCF"

  echo "chr${i}: FIXREF $(date +%H:%M:%S)"
  # CRITICAL: fixref stderr goes to stats log, NOT into the pipe
  bcftools +fixref "$TMPVCF" -Ou -- -f "$REF" -m top 2>>"$STATS" | \
    bcftools view -e 'TYPE="snp" && ((REF="A" && ALT="T") || (REF="T" && ALT="A") || (REF="C" && ALT="G") || (REF="G" && ALT="C"))' \
    -Oz -o "$OUT"

  if [ ! -s "$OUT" ]; then
    echo "chr${i}: FIXREF_FAIL (empty output)"
    rm -f "$OUT"
    rm -f "$TMPVCF" "${TMPVCF}.csi" "$REEXPORT"
    continue
  fi

  bcftools index -t "$OUT"
  n=$(bcftools index -n "$OUT")
  echo "chr${i}: DONE ${n} $(date +%H:%M:%S)"
  rm -f "$TMPVCF" "${TMPVCF}.csi" "$REEXPORT" "/tmp/reexport_chr${i}.log" "/tmp/reexport_chr${i}.nosex"
done

echo ""
echo "=== FINAL SUMMARY ==="
total=0
ok=0
for i in $(seq 1 22); do
  f="${READYDIR}/impute_chr${i}.vcf.gz"
  if [ -f "${f}.tbi" ]; then
    n=$(bcftools index -n "$f" 2>/dev/null || echo 0)
    echo "chr${i}: ${n}"
    total=$((total + n))
    ok=$((ok + 1))
  else
    echo "chr${i}: MISSING"
  fi
done
echo ""
echo "RESULT: ${ok}/22 chromosomes, ${total} total variants"

if [ "$ok" -eq 22 ]; then
  echo "ALL COMPLETE"
  cd /tmp
  tar cf michigan_fixref.tar -C michigan_ready .
  ls -lh michigan_fixref.tar
fi
echo "END $(date)"
