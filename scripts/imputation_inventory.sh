#!/bin/bash
cd /staging/ALSU-analysis/spring2026/imputation 2>/dev/null || { echo "NODIR"; exit 1; }

echo "=== LISTING ==="
ls -lhS | head -80

echo "=== FILE_COUNTS ==="
echo "zip: $(ls *.zip 2>/dev/null | wc -l)"
echo "dose_vcf: $(ls chr*.dose.vcf.gz 2>/dev/null | wc -l)"
echo "info_gz: $(ls chr*.info.gz 2>/dev/null | wc -l)"
echo "tbi: $(ls *.tbi 2>/dev/null | wc -l)"
echo "txt: $(ls *.txt 2>/dev/null | wc -l)"
echo "html: $(ls *.html 2>/dev/null | wc -l)"

echo "=== QC_REPORT ==="
cat qc_report.txt 2>/dev/null || echo "NO_QC_REPORT"

echo "=== STATISTICS ==="
cat statistics.txt 2>/dev/null || echo "NO_STATISTICS"

echo "=== SAMPLE_COUNT ==="
f=$(ls chr*.dose.vcf.gz 2>/dev/null | head -1)
if [ -n "$f" ]; then
  bcftools query -l "$f" | wc -l
else
  echo "NO_DOSE_VCF"
fi

echo "=== PER_CHR_VARIANTS ==="
for chr in $(seq 1 22); do
  f="chr${chr}.dose.vcf.gz"
  if [ -f "$f" ]; then
    n=$(bcftools index -n "$f" 2>/dev/null || bcftools view -H "$f" | wc -l)
    echo "chr${chr} ${n}"
  else
    echo "chr${chr} MISSING"
  fi
done

echo "=== INFO_R2_PER_CHR ==="
for chr in $(seq 1 22); do
  f="chr${chr}.info.gz"
  if [ -f "$f" ]; then
    zcat "$f" | grep -v '^#' | awk -v c=$chr -F'\t' '
    {
      split($8,kv,";")
      for(i in kv){
        if(kv[i]~/^R2=/){r2=substr(kv[i],4)+0}
      }
      n++; s+=r2
      b=int(r2*10); if(b>=10)b=9; bins[b]++
      if(r2>=0.3) ge3++
      if(r2>=0.8) ge8++
    }
    END {
      printf "CHR%d N=%d MEAN_R2=%.4f GE03=%d GE08=%d", c, n, s/n, ge3, ge8
      for(i=0;i<10;i++) printf " B%d=%d", i, bins[i]+0
      printf "\n"
    }'
  else
    echo "chr${chr}.info.gz MISSING"
  fi
done

echo "=== TOTAL_DISK ==="
du -sh /staging/ALSU-analysis/spring2026/imputation/ 2>/dev/null

echo "=== DONE ==="
