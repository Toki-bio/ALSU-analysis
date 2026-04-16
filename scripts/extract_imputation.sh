#!/bin/bash
# Extract Michigan imputation results and prepare for post-imputation pipeline
# Run on DRAGEN: /staging/ALSU-analysis/spring2026/imputation/
set -euo pipefail

cd /staging/ALSU-analysis/spring2026/imputation

echo "=== DISK SPACE ==="
df -h /staging/

echo "=== EXTRACTING DOSE + INFO FILES ==="
# Password from Michigan download email
IMPUTE_PASSWORD="d6rFCYsg&q;RZ0"

for chr in $(seq 1 22); do
  f="chr_${chr}.zip"
  if [ ! -f "$f" ]; then
    echo "MISSING: $f"
    continue
  fi
  # Extract only needed files (dose VCF + info + empiricalDose)
  # Skip phased VCFs to save disk (can re-extract later if needed)
  echo "Extracting chr${chr}..."
  unzip -P "$IMPUTE_PASSWORD" -o "$f" \
    "chr${chr}.dose.vcf.gz" \
    "chr${chr}.info.gz" \
    2>/dev/null || echo "WARN: partial extract chr${chr}"
done

echo "=== VERIFY EXTRACTION ==="
echo "dose_vcf: $(ls chr*.dose.vcf.gz 2>/dev/null | wc -l)"
echo "info_gz: $(ls chr*.info.gz 2>/dev/null | wc -l)"

echo "=== SAMPLE COUNT ==="
bcftools query -l chr1.dose.vcf.gz | wc -l
echo "First 5 sample IDs:"
bcftools query -l chr1.dose.vcf.gz | head -5

echo "=== PER-CHR VARIANT COUNTS ==="
total=0
for chr in $(seq 1 22); do
  f="chr${chr}.dose.vcf.gz"
  if [ -f "$f" ]; then
    # Index first for fast counting
    bcftools index -f "$f"
    n=$(bcftools index -n "$f")
    total=$((total + n))
    echo "chr${chr}: ${n}"
  fi
done
echo "TOTAL: ${total}"

echo "=== SAMPLE ID CHECK ==="
# Check for Michigan numeric prefix (e.g., "0_sampleID")
bcftools query -l chr1.dose.vcf.gz | grep -cP '^\d+_' || echo "0 (no numeric prefix)"
# Check for non-ASCII
bcftools query -l chr1.dose.vcf.gz | grep -cP '[^\x00-\x7F]' || echo "0 (all ASCII)"

echo "=== R2 DISTRIBUTION (from .info.gz) ==="
# Compute aggregate R² stats from all chromosomes
for chr in $(seq 1 22); do
  f="chr${chr}.info.gz"
  [ -f "$f" ] && zcat "$f" | grep -v '^#'
done | awk -F'\t' '
{
  split($8,kv,";")
  for(i in kv){
    if(kv[i]~/^R2=/){r2=substr(kv[i],4)+0}
  }
  n++; s+=r2
  b=int(r2*10); if(b>=10) b=9; bins[b]++
  if(r2>=0.3) ge3++
  if(r2>=0.8) ge8++
  if(r2>=0.9) ge9++
}
END {
  printf "Total=%d Mean_R2=%.4f\n", n, s/n
  printf "R2>=0.30: %d (%.1f%%)\n", ge3, ge3/n*100
  printf "R2>=0.80: %d (%.1f%%)\n", ge8, ge8/n*100
  printf "R2>=0.90: %d (%.1f%%)\n", ge9, ge9/n*100
  for(i=0;i<10;i++) printf "Bin[%.1f-%.1f]: %d\n", i/10, (i+1)/10, bins[i]+0
}'

echo "=== DISK AFTER EXTRACTION ==="
du -sh /staging/ALSU-analysis/spring2026/imputation/
df -h /staging/

echo "=== DONE ==="
