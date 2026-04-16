#!/bin/bash
set -uo pipefail
cd /staging/ALSU-analysis/spring2026/imputation

echo "=== DISK BEFORE ==="
df -h /staging/

echo "=== EXTRACTING ==="
PW='d6rFCYsg&q;RZ0'

for chr in $(seq 1 22); do
  f="chr_${chr}.zip"
  if [ ! -f "$f" ]; then echo "MISSING: $f"; continue; fi
  echo -n "chr${chr}... "
  unzip -P "$PW" -o "$f" "chr${chr}.dose.vcf.gz" "chr${chr}.info.gz" 2>&1 | grep -c 'inflating' || echo "WARN"
done

echo "=== VERIFY FILES ==="
echo "dose_vcf: $(ls chr*.dose.vcf.gz 2>/dev/null | wc -l)"
echo "info_gz: $(ls chr*.info.gz 2>/dev/null | wc -l)"
ls -lhS chr*.dose.vcf.gz 2>/dev/null | head -25

echo "=== SAMPLES ==="
bcftools query -l chr1.dose.vcf.gz | wc -l
echo "First 5:"
bcftools query -l chr1.dose.vcf.gz | head -5
echo "Numeric prefix check:"
bcftools query -l chr1.dose.vcf.gz | grep -cP '^\d+_' || echo "0"
echo "Non-ASCII check:"
bcftools query -l chr1.dose.vcf.gz | grep -cP '[^\x00-\x7F]' || echo "0"

echo "=== INDEX + PER-CHR COUNTS ==="
total=0
for chr in $(seq 1 22); do
  f="chr${chr}.dose.vcf.gz"
  if [ -f "$f" ]; then
    bcftools index -f "$f" 2>/dev/null
    n=$(bcftools index -n "$f" 2>/dev/null)
    total=$((total + n))
    echo "chr${chr}: ${n}"
  else
    echo "chr${chr}: MISSING"
  fi
done
echo "TOTAL: ${total}"

echo "=== R2 PER-CHR ==="
for chr in $(seq 1 22); do
  f="chr${chr}.info.gz"
  if [ -f "$f" ]; then
    zcat "$f" | grep -v '^#' | awk -v c=$chr -F'\t' '
    {
      split($8,kv,";")
      for(i in kv) if(kv[i]~/^R2=/) r2=substr(kv[i],4)+0
      n++; s+=r2
      b=int(r2*10); if(b>=10)b=9; bins[b]++
      if(r2>=0.3)ge3++; if(r2>=0.8)ge8++; if(r2>=0.9)ge9++
    }
    END {
      printf "CHR%d N=%d MEAN=%.4f GE03=%d GE08=%d GE09=%d", c,n,s/n,ge3+0,ge8+0,ge9+0
      for(i=0;i<10;i++) printf " B%d=%d",i,bins[i]+0
      printf "\n"
    }'
  fi
done

echo "=== R2 AGGREGATE ==="
for chr in $(seq 1 22); do
  f="chr${chr}.info.gz"
  [ -f "$f" ] && zcat "$f" | grep -v '^#'
done | awk -F'\t' '
{
  split($8,kv,";")
  for(i in kv) if(kv[i]~/^R2=/) r2=substr(kv[i],4)+0
  n++; s+=r2
  b=int(r2*10); if(b>=10)b=9; bins[b]++
  if(r2>=0.3) ge3++
  if(r2>=0.8) ge8++
  if(r2>=0.9) ge9++
}
END {
  printf "TOTAL N=%d MEAN_R2=%.4f\n", n, s/n
  printf "R2>=0.30: %d (%.1f%%)\n", ge3+0, (ge3+0)/n*100
  printf "R2>=0.80: %d (%.1f%%)\n", ge8+0, (ge8+0)/n*100
  printf "R2>=0.90: %d (%.1f%%)\n", ge9+0, (ge9+0)/n*100
  for(i=0;i<10;i++) printf "Bin[0.%d-0.%d]: %d\n", i, i+1, bins[i]+0
}'

echo "=== DISK AFTER ==="
df -h /staging/
du -sh /staging/ALSU-analysis/spring2026/imputation/

echo "=== DONE ==="
