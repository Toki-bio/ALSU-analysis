#!/bin/bash
set -eo pipefail
PW='d6rFCYsg&q;RZ0'
WORKDIR=/tmp/michigan_results
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Check tools
echo "=== Tool check ==="
which 7z 2>/dev/null && echo "7z: yes" || echo "7z: no"
which unzip 2>/dev/null && echo "unzip: yes" || echo "unzip: no"
df -h /tmp | tail -1

# Delete incomplete chr_2
rm -f chr_2.zip
echo "Deleted incomplete chr_2.zip"
df -h /tmp | tail -1

# Test: extract just .info.gz from the complete chr_1.zip
echo ""
echo "=== Testing extraction of chr1.info.gz from chr_1.zip ==="
if command -v 7z &>/dev/null; then
  7z l -p"$PW" chr_1.zip 2>/dev/null | head -30
  echo "---"
  7z e -p"$PW" -aoa chr_1.zip "chr1.info.gz" 2>/dev/null && echo "Extracted chr1.info.gz" || echo "7z extract failed"
elif command -v unzip &>/dev/null; then
  unzip -l chr_1.zip 2>/dev/null | head -20
  echo "---"
  unzip -P "$PW" -o chr_1.zip "chr1.info.gz" 2>/dev/null && echo "Extracted chr1.info.gz" || echo "unzip extract failed"
fi

echo ""
echo "=== Extracted files ==="
ls -lh chr1.info.gz 2>/dev/null || echo "no chr1.info.gz"

echo ""
echo "=== Check format ==="
zcat chr1.info.gz 2>/dev/null | head -3 || file chr1.info.gz 2>/dev/null

echo "=== DONE ==="
