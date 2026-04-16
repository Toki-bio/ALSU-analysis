#!/bin/bash
KG=/staging/ALSU-analysis/winter2025/PLINK_301125_0312/step9_1000g.pca/1000g_grch38/KG_reference_final

echo "========== PRE-ADMIXTURE AUDIT =========="

echo "=== 1. KG_reference BIM variant ID format ==="
head -5 ${KG}.bim
echo "..."
tail -5 ${KG}.bim
echo ""

echo "=== 2. KG_reference FAM: FID distribution ==="
cut -d' ' -f1 ${KG}.fam | sort | uniq -c | sort -rn | head -10
echo ""

echo "=== 3. KG_reference FAM: total and sample name pattern ==="
wc -l ${KG}.fam
head -3 ${KG}.fam
tail -3 ${KG}.fam
echo ""

echo "=== 4. Pop mapping: v2/pca ==="
if [ -f ~/v2/pca/pop_mapping.txt ]; then
  head -5 ~/v2/pca/pop_mapping.txt
  echo "field count per line:"
  head -1 ~/v2/pca/pop_mapping.txt | tr ' ' '\n' | wc -l
  echo "unique col3 values:"
  awk '{print $3}' ~/v2/pca/pop_mapping.txt | sort -u
  echo "col3 counts:"
  awk '{print $3}' ~/v2/pca/pop_mapping.txt | sort | uniq -c | sort -rn
else
  echo "NOT FOUND"
fi
echo ""

echo "=== 5. Pop mapping: v2/global_pca ==="
if [ -f ~/v2/global_pca/pop_mapping.txt ]; then
  head -5 ~/v2/global_pca/pop_mapping.txt
  awk '{print $3}' ~/v2/global_pca/pop_mapping.txt | sort | uniq -c | sort -rn
else
  echo "NOT FOUND"
fi
echo ""

echo "=== 6. Clusters that actually worked (v2/pbs_v2) ==="
if [ -f ~/v2/pbs_v2/clusters.txt ]; then
  head -5 ~/v2/pbs_v2/clusters.txt
  echo ""
  awk '{print $3}' ~/v2/pbs_v2/clusters.txt | sort | uniq -c | sort -rn
else
  echo "NOT FOUND"
fi
echo ""

echo "=== 7. 1000G panel file ==="
find /staging/ALSU-analysis -name '*panel*' -type f 2>/dev/null | head -5
PANEL=$(find /staging/ALSU-analysis -name '*panel*' -type f 2>/dev/null | head -1)
if [ -n "$PANEL" ]; then
  echo "Content:"
  head -3 "$PANEL"
  echo "Superpop counts:"
  awk 'NR>1 {print $3}' "$PANEL" | sort | uniq -c | sort -rn
fi
echo ""

echo "=== 8. KG_reference_final build log ==="
if [ -f ${KG}.log ]; then
  head -40 ${KG}.log
else
  echo "NO LOG"
fi
echo ""

echo "=== 9. What step6 will produce: expected ID format ==="
echo "Step5 VCF sample names (first 5):"
zcat /staging/ALSU-analysis/spring2026/post_imputation/chr22.dose.vcf.gz 2>/dev/null | head -1000 | grep '^#CHROM' | tr '\t' '\n' | tail -n +10 | head -5
echo ""
echo "Step5 VCF variant format (first 5 data lines):"
zcat /staging/ALSU-analysis/spring2026/post_imputation/chr22.dose.vcf.gz 2>/dev/null | grep -v '^#' | head -5 | cut -f1-5
echo ""

echo "=== 10. Step6 progress ==="
tail -5 /tmp/step6_out.txt 2>/dev/null
echo ""
echo "Step6 process:"
ps aux | grep step6 | grep -v grep
echo ""
echo "bcftools:"
ps aux | grep bcftools | grep -v grep
echo ""

echo "=== 11. Disk space ==="
df -h /staging/
echo ""

echo "=== 12. admixture binary ==="
which admixture 2>/dev/null || echo "NOT IN PATH"
admixture 2>&1 | head -3
echo ""

echo "=== 13. plink2 version ==="
plink2 --version 2>/dev/null | head -1
echo ""

echo "========== AUDIT COMPLETE =========="
