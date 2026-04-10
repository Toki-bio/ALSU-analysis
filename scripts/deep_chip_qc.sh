#!/bin/bash
# Deep dive into bad vs good chip QC data, control probes, scan quality
set -e

CONV="/staging/ALSU-analysis/Conversion/OUTPUT"
BAD="208993030034 208993030080 208993030112 208993030109"
GOOD="207591080076 207585740001 207591070002"

echo "=== 1. QC FILES — BAD CHIPS ==="
for chip in $BAD; do
    echo "--- ${chip}_qc.txt ---"
    cat "$CONV/$chip/${chip}_qc.txt" 2>/dev/null || echo "(not found)"
    echo ""
done

echo "=== 2. QC FILES — GOOD CHIPS (for comparison) ==="
for chip in $GOOD; do
    echo "--- ${chip}_qc.txt ---"
    cat "$CONV/$chip/${chip}_qc.txt" 2>/dev/null || echo "(not found)"
    echo ""
done

echo "=== 3. SDF FILES — BAD vs GOOD ==="
for chip in $BAD $GOOD; do
    echo "--- $chip.sdf (first 30 lines) ---"
    head -30 "$CONV/$chip/${chip}.sdf" 2>/dev/null || echo "(not found)"
    echo ""
done

echo "=== 4. FILE LISTING — BAD CHIP (full) ==="
for chip in $BAD; do
    echo "--- $chip ---"
    ls -la "$CONV/$chip/" 2>/dev/null | head -80
    echo ""
done

echo "=== 5. FILE LISTING — GOOD CHIP (full) ==="
ls -la "$CONV/207591080076/" 2>/dev/null | head -80

echo ""
echo "=== 6. IDAT FILE COUNTS AND SIZES ==="
for chip in $BAD $GOOD; do
    n=$(find "$CONV/$chip/" -name "*.idat" 2>/dev/null | wc -l)
    sz=$(find "$CONV/$chip/" -name "*.idat" -exec du -ch {} + 2>/dev/null | tail -1 || echo "0")
    echo "$chip: $n idat files, total $sz"
done

echo ""
echo "=== 7. CONTROL PROBES FROM MANIFEST ==="
grep -A 200 "^\[Controls\]" /staging/ALSU-analysis/Conversion/GSA-24v3-0_A2.csv 2>/dev/null | head -60

echo ""
echo "=== 8. XML SCAN METADATA — BAD CHIP SAMPLE ==="
# Grab one XML from a bad chip to see scan parameters
head -50 "$CONV/208993030034/208993030034_R01C01_1_Green.xml" 2>/dev/null || echo "(not found)"

echo ""
echo "=== 9. XML SCAN METADATA — GOOD CHIP SAMPLE ==="
head -50 "$CONV/207591080076/207591080076_R01C01_1_Green.xml" 2>/dev/null || echo "(not found)"

echo ""
echo "=== 10. SCAN JPEG SIZES — BAD vs GOOD ==="
echo "BAD chip 208993030034:"
ls -la "$CONV/208993030034/"*.jpg 2>/dev/null | awk '{print $5, $9}' | head -20
echo ""
echo "GOOD chip 207591080076:"
ls -la "$CONV/207591080076/"*.jpg 2>/dev/null | awk '{print $5, $9}' | head -20

echo ""
echo "=== 11. CHECK diag DIRECTORY ==="
ls -la "$CONV/diag/" 2>/dev/null | head -20
ls -la "$CONV/diag/207591080076/" 2>/dev/null | head -20

echo ""
echo "=== 12. BSC FILE (ConvSK.bsc) — first 20 lines ==="
head -20 "$CONV/ConvSK/ConvSK.bsc" 2>/dev/null || echo "(not found)"

echo ""
echo "DONE"
