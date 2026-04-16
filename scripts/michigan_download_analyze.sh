#!/bin/bash
set -eo pipefail

WORKDIR=/tmp/michigan_results
mkdir -p "$WORKDIR/statistics"
cd "$WORKDIR"

PW='d6rFCYsg&q;RZ0'

echo "=== Michigan Download+Decrypt+Analyze: $(date) ==="
echo "Working directory: $WORKDIR"

# ── 1. Download chromosome zips ──
declare -a U
U[1]="https://imputationserver.sph.umich.edu/share/results/2a87b30430bbd44939595486534fcb9df55426fd81789165d9077e178daccdff/chr_1.zip"
U[2]="https://imputationserver.sph.umich.edu/share/results/ec0f1ecdd1ae72b052894113a2e77f09f8143ac2d6d6e08b6d0e78c02b62bbbf/chr_2.zip"
U[3]="https://imputationserver.sph.umich.edu/share/results/b33f936ad78fa65a3d1d6d7f977c5b5d2836eba61d08cee0b445ae97e47ee887/chr_3.zip"
U[4]="https://imputationserver.sph.umich.edu/share/results/9c756df94ae17767799cdb4e422f16bba0e141c44eb788c33ba3705780d71082/chr_4.zip"
U[5]="https://imputationserver.sph.umich.edu/share/results/00add947fbf01535557af01de599813357bb4bb102ffcce528175d5c03051be8/chr_5.zip"
U[6]="https://imputationserver.sph.umich.edu/share/results/56ae398622942bb46c32f45363771a820b31990d5cff703d43b360ef814f990e/chr_6.zip"
U[7]="https://imputationserver.sph.umich.edu/share/results/d567576ebc4ff599b584002f33d7e7274d41a08e97923ef2a0ef5fbccb883e0b/chr_7.zip"
U[8]="https://imputationserver.sph.umich.edu/share/results/786395b7feab0acfb8a5f1982c5646fadde2b0264acbf3298e44147a67c6d9d3/chr_8.zip"
U[9]="https://imputationserver.sph.umich.edu/share/results/b9cb5a89a9ffc65d3dadc30732f1cefb8b6b0f4cca50f38876c0636964fe5401/chr_9.zip"
U[10]="https://imputationserver.sph.umich.edu/share/results/e50af01426196de2556bdcce7a4412b1d94d62d8a5fcf06e9d81cefa6f1ccabb/chr_10.zip"
U[11]="https://imputationserver.sph.umich.edu/share/results/d7ede4954d11b2e57880f8570b5939dfc163e4109d32897f4f6e40d22bcc75f7/chr_11.zip"
U[12]="https://imputationserver.sph.umich.edu/share/results/8c20885eef6e66c5656d19663e30772714b4bba31509a7d66dc2dc81bbf4f5f7/chr_12.zip"
U[13]="https://imputationserver.sph.umich.edu/share/results/c469a0d38d7e292f25af6c69a295365fe75910c3bd7f82065f98e0ef3bf61194/chr_13.zip"
U[14]="https://imputationserver.sph.umich.edu/share/results/5b107fa2f3eb5eb1a460e65ff1bc9e7fb74b6a9d8e1e976b986bb1aba11c3757/chr_14.zip"
U[15]="https://imputationserver.sph.umich.edu/share/results/82469309f2752b1dce8bbdff97de4c373be06e3eeef5d45c2e48c9febab1b947/chr_15.zip"
U[16]="https://imputationserver.sph.umich.edu/share/results/9171ef87052810bb33ee58a1798c8fdbf9e0606a28dfef2450fc5870394b14d2/chr_16.zip"
U[17]="https://imputationserver.sph.umich.edu/share/results/8c497ee4bea8569f9b87d334c97022bab7ddaf97f86a23932b5d28f5cb20c85d/chr_17.zip"
U[18]="https://imputationserver.sph.umich.edu/share/results/4d78b2e898814c3135a75bccd005e2f33365406a9ac0a449a7e6845dc8fe8043/chr_18.zip"
U[19]="https://imputationserver.sph.umich.edu/share/results/47566eca76332d15b02cdf8eab040ee5bc968e9fe99a088a8a2090cd288366f1/chr_19.zip"
U[20]="https://imputationserver.sph.umich.edu/share/results/e7e3974d1fb72c7d46f700582645db69da15b4471002887a7a3f941398dd7160/chr_20.zip"
U[21]="https://imputationserver.sph.umich.edu/share/results/674be72e77b0690630a9a25f53b2e57a083d93a7f2ad5874c6e42be8c8ef6727/chr_21.zip"
U[22]="https://imputationserver.sph.umich.edu/share/results/9bb42010b77064e9e0c8b2d96d6d1fbce359a91d3636addebd815fb6f344da57/chr_22.zip"

echo "=== Downloading 22 chromosome zips ==="
for chr in $(seq 1 22); do
  out="chr_${chr}.zip"
  if [ -f "$out" ] && [ "$(stat -c%s "$out" 2>/dev/null)" -gt 1000000 ]; then
    echo "SKIP $out ($(du -h "$out" | cut -f1) exists)"
  else
    echo -n "Downloading $out ... "
    curl -sL -o "$out" "${U[$chr]}"
    echo "$(du -h "$out" | cut -f1)"
  fi
done

# QC files
curl -sL -o qc_report.txt "https://imputationserver.sph.umich.edu/share/results/b6f88356af6ca105e8170aa9baf694739683793ef64633db611fe7d02f25442f/qc_report.txt"
curl -sL -o "statistics/chunks-excluded.txt" "https://imputationserver.sph.umich.edu/share/results/5d697668a3cd0d5bd62b728e5ab61bec415d8d630f632b570beee788df15efe8/statistics/chunks-excluded.txt"
curl -sL -o "statistics/snps-excluded.txt" "https://imputationserver.sph.umich.edu/share/results/1194b94cb7852fd808ff6258f5ef4e747c030b3cc1ae6c4a02eb799f6919efc2/statistics/snps-excluded.txt"
curl -sL -o "statistics/snps-typed-only.txt" "https://imputationserver.sph.umich.edu/share/results/35e268282089532828d645e5f47b3f2f582c065d6d1bb5fd99011538214db8fb/statistics/snps-typed-only.txt"
echo "QC files done"

echo ""
echo "=== Downloads complete: $(date) ==="
echo ""

# ── 2. Decrypt all zips ──
echo "=== Decrypting ==="
for chr in $(seq 1 22); do
  if [ -f "chr${chr}.dose.vcf.gz" ]; then
    echo "chr${chr} already decrypted"
    continue
  fi
  if command -v 7z &>/dev/null; then
    7z x -p"$PW" -aoa "chr_${chr}.zip" >/dev/null 2>&1 && echo "chr_${chr} decrypted (7z)" || echo "chr_${chr} DECRYPT FAILED (7z)"
  else
    unzip -P "$PW" -o "chr_${chr}.zip" >/dev/null 2>&1 && echo "chr_${chr} decrypted (unzip)" || echo "chr_${chr} DECRYPT FAILED (unzip)"
  fi
done
echo "=== Decryption complete: $(date) ==="
echo ""

# ── 3. Verify files exist ──
echo "=== Extracted files ==="
ls -lh chr*.dose.vcf.gz 2>/dev/null | awk '{print $NF, $5}'
ls -lh chr*.info.gz 2>/dev/null | awk '{print $NF, $5}'
echo ""

# ── 4. Info header ──
echo "=== INFO file header ==="
zcat chr1.info.gz 2>/dev/null | head -2
echo ""

# ── 5. Sample count ──
echo "=== Sample count ==="
bcftools query -l chr1.dose.vcf.gz 2>/dev/null | wc -l
echo ""

# ── 6. Per-chromosome variant counts ──
echo "=== Per-chromosome variant counts ==="
printf "%-6s %10s %10s\n" "CHR" "TOTAL" "TYPED"
grand_total=0
grand_typed=0
for chr in $(seq 1 22); do
  if [ -f "chr${chr}.info.gz" ]; then
    read total typed < <(zcat "chr${chr}.info.gz" | tail -n +2 | awk '
      {n++; if ($8=="Typed" || $8=="Genotyped") t++}
      END {printf "%d %d\n", n, t}
    ')
    grand_total=$((grand_total + total))
    grand_typed=$((grand_typed + typed))
    printf "chr%-4s %10d %10d\n" "$chr" "$total" "$typed"
  else
    printf "chr%-4s %10s %10s\n" "$chr" "MISSING" "MISSING"
  fi
done
printf "%-6s %10d %10d\n" "TOTAL" "$grand_total" "$grand_typed"
echo ""

# ── 7. INFO/R² distribution ──
echo "=== INFO/R² distribution ==="
for chr in $(seq 1 22); do
  [ -f "chr${chr}.info.gz" ] && zcat "chr${chr}.info.gz" | tail -n +2 | awk '{print $7}'
done | awk '
BEGIN {
  for (i=0; i<10; i++) bins[i] = 0
  n = 0; sum = 0
  ge30 = 0; ge50 = 0; ge80 = 0; ge90 = 0
}
{
  v = $1 + 0
  sum += v; n++
  bin = int(v * 10)
  if (bin >= 10) bin = 9
  bins[bin]++
  if (v >= 0.3) ge30++
  if (v >= 0.5) ge50++
  if (v >= 0.8) ge80++
  if (v >= 0.9) ge90++
}
END {
  printf "Total variants: %d\n", n
  printf "Mean R²: %.6f\n", (n>0 ? sum/n : 0)
  printf "\nHistogram (10 bins):\n"
  for (i=0; i<10; i++) {
    lo = i/10; hi = (i+1)/10
    pct = (n>0 ? bins[i]/n*100 : 0)
    printf "  %.1f-%.1f: %12d  (%5.1f%%)\n", lo, hi, bins[i], pct
  }
  printf "\nThresholds:\n"
  printf "  R² >= 0.30: %12d  (%5.1f%%)\n", ge30, (n>0 ? ge30/n*100 : 0)
  printf "  R² >= 0.50: %12d  (%5.1f%%)\n", ge50, (n>0 ? ge50/n*100 : 0)
  printf "  R² >= 0.80: %12d  (%5.1f%%)\n", ge80, (n>0 ? ge80/n*100 : 0)
  printf "  R² >= 0.90: %12d  (%5.1f%%)\n", ge90, (n>0 ? ge90/n*100 : 0)
}
' | tee info_distribution.txt
echo ""

# ── 8. Typed vs Imputed counts ──
echo "=== Typed vs Imputed ==="
for chr in $(seq 1 22); do
  [ -f "chr${chr}.info.gz" ] && zcat "chr${chr}.info.gz" | tail -n +2 | awk '{print $8}'
done | sort | uniq -c | sort -rn
echo ""

# ── 9. MAF distribution of imputed variants ──
echo "=== MAF distribution (imputed only) ==="
for chr in $(seq 1 22); do
  [ -f "chr${chr}.info.gz" ] && zcat "chr${chr}.info.gz" | tail -n +2 | awk '$8!="Typed" && $8!="Genotyped" {print $5}'
done | awk '
BEGIN { for (i=0; i<10; i++) bins[i] = 0 }
{
  v = $1 + 0
  bin = int(v * 20)  # 5% bins up to 50%
  if (bin >= 10) bin = 9
  bins[bin]++; n++
}
END {
  printf "MAF distribution (imputed, N=%d):\n", n
  labels[0]="0-5%"; labels[1]="5-10%"; labels[2]="10-15%"; labels[3]="15-20%"; labels[4]="20-25%"
  labels[5]="25-30%"; labels[6]="30-35%"; labels[7]="35-40%"; labels[8]="40-45%"; labels[9]="45-50%"
  for (i=0; i<10; i++) printf "  %s: %d (%.1f%%)\n", labels[i], bins[i], bins[i]/n*100
}
'
echo ""

# ── 10. Winter 2025 comparison ──
echo "=== Winter 2025 comparison ==="
winter_bim="/staging/ALSU-analysis/winter2025/3_post-imputation/ConvSK_final_clean.bim"
winter_fam="/staging/ALSU-analysis/winter2025/3_post-imputation/ConvSK_final_clean.fam"
if [ -f "$winter_bim" ]; then
  w_vars=$(wc -l < "$winter_bim")
  w_samp=$(wc -l < "$winter_fam" 2>/dev/null || echo "N/A")
  echo "Winter 2025 post-imputation PLINK (after INFO filtering):"
  echo "  Variants: $w_vars"
  echo "  Samples: $w_samp"
else
  echo "Winter 2025 PLINK files not found at $winter_bim"
fi
echo ""

# Check if winter raw imputation files exist
echo "Winter raw imputed VCFs:"
ls /staging/ALSU-analysis/winter2025/PLINK_301125_0312/imputed/*.dose.vcf.gz 2>/dev/null | wc -l | xargs -I{} echo "{} files"
ls /staging/ALSU-analysis/winter2025/PLINK_301125_0312/imputed/*.info.gz 2>/dev/null | wc -l | xargs -I{} echo "{} info files"

# If winter info files exist, get quick summary
winter_info="/staging/ALSU-analysis/winter2025/PLINK_301125_0312/imputed/chr1.info.gz"
if [ -f "$winter_info" ]; then
  echo ""
  echo "Winter INFO distribution (quick summary):"
  for chr in $(seq 1 22); do
    f="/staging/ALSU-analysis/winter2025/PLINK_301125_0312/imputed/chr${chr}.info.gz"
    [ -f "$f" ] && zcat "$f" | tail -n +2 | awk '{print $7}'
  done | awk '
  { n++; sum += $1; if ($1>=0.3) a++; if ($1>=0.8) b++; if ($1>=0.9) c++ }
  END {
    printf "  Total: %d, Mean R²: %.4f\n", n, sum/n
    printf "  >=0.3: %d (%.1f%%), >=0.8: %d (%.1f%%), >=0.9: %d (%.1f%%)\n", a,a/n*100, b,b/n*100, c,c/n*100
  }'
fi

echo ""
echo "=== COMPLETE: $(date) ==="
