#!/bin/bash
# Diagnostic: what's in spring2026 and where did chr{1..22}.vcf.gz come from?
cd /staging/ALSU-analysis/spring2026/

echo "=== Directory listing ==="
ls -lh *.bed *.bim *.fam 2>/dev/null | head -20

echo ""
echo "=== VCF files ==="
ls -lh chr*.vcf.gz 2>/dev/null | head -25

echo ""
echo "=== Variant counts per chr VCF (first 3 + total) ==="
for chr in 1 2 3; do
  n=$(bcftools view -H chr${chr}.vcf.gz 2>/dev/null | wc -l)
  echo "chr${chr}: ${n} variants"
done

echo ""
echo "=== Total across all chr VCFs ==="
total=0
for chr in $(seq 1 22); do
  n=$(bcftools view -H chr${chr}.vcf.gz 2>/dev/null | wc -l)
  total=$((total + n))
done
echo "Total: ${total}"

echo ""
echo "=== Sample count in chr1 VCF ==="
bcftools query -l chr1.vcf.gz 2>/dev/null | wc -l

echo ""
echo "=== VCF header check (chr1, first 30 lines) ==="
bcftools view -h chr1.vcf.gz 2>/dev/null | head -30

echo ""
echo "=== PLINK files present ==="
for prefix in ConvSK_raw ConvSK_mind20 ConvSK_mind20_dedup ConvSK_mind20_dedup_snpqc; do
  if [ -f "${prefix}.bed" ]; then
    fam_n=$(wc -l < ${prefix}.fam)
    bim_n=$(wc -l < ${prefix}.bim)
    echo "${prefix}: ${bim_n} variants x ${fam_n} samples"
  else
    echo "${prefix}: NOT FOUND"
  fi
done

echo ""
echo "=== impute_in files ==="
for f in impute_in.autosomes.vcf.gz impute_in.autosomes.clean.vcf.gz impute_in.autosomes.clean.nodup.vcf.gz impute_in.final.vcf.gz; do
  if [ -f "$f" ]; then
    ls -lh "$f"
  fi
done
