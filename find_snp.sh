#!/bin/bash
# Search for rs5822325 / chr17:80107730 across all .bim files
echo "=== Searching for rs5822325 / chr17:80107730 ==="
echo "--- Finding .bim files ---"
find /home /data /staging /work /tmp -name "*.bim" -type f 2>/dev/null | head -50

echo ""
echo "--- Searching by rsID (rs5822325) ---"
for f in $(find /home /data /staging /work /tmp -name "*.bim" -type f 2>/dev/null | head -50); do
  hit=$(grep -i "5822325" "$f" 2>/dev/null)
  if [ -n "$hit" ]; then
    echo "FOUND in $f:"
    echo "$hit"
  fi
done

echo ""
echo "--- Searching by position chr17:80107730 ---"
for f in $(find /home /data /staging /work /tmp -name "*.bim" -type f 2>/dev/null | head -50); do
  hit=$(awk '$1==17 && $4==80107730' "$f" 2>/dev/null)
  if [ -n "$hit" ]; then
    echo "FOUND in $f:"
    echo "$hit"
  fi
done

echo ""
echo "--- Checking for VCF files with this position ---"
for f in $(find /home /data /staging /work /tmp -name "*.vcf.gz" -type f 2>/dev/null | head -20); do
  echo "VCF: $f"
done

echo "=== DONE ==="
