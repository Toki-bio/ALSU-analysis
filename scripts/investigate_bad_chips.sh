#!/bin/bash
# Investigate the 4 catastrophic chips: idats, metadata, controls, scan images
set -e

CONV="/staging/ALSU-analysis/Conversion/OUTPUT"
BAD_CHIPS="208993030034 208993030080 208993030112 208993030109"
# Also check a known-good chip for comparison
GOOD_CHIP="207591080076"

echo "=== 1. DIRECTORY STRUCTURE OF CONVERSION OUTPUT ==="
ls -la "$CONV/" 2>/dev/null | head -20

echo ""
echo "=== 2. LOOK FOR IDAT FILES FOR BAD CHIPS ==="
for chip in $BAD_CHIPS; do
    echo "--- $chip ---"
    find "$CONV" -maxdepth 3 -name "${chip}*" -type f 2>/dev/null | head -20
done

echo ""
echo "=== 3. LOOK FOR IDAT FILES FOR GOOD CHIP (compare) ==="
find "$CONV" -maxdepth 3 -name "${GOOD_CHIP}*" -type f 2>/dev/null | head -10

echo ""
echo "=== 4. CHECK FOR JPEG/PNG SCAN IMAGES ==="
find "$CONV" -maxdepth 3 \( -name "*.jpg" -o -name "*.jpeg" -o -name "*.png" -o -name "*.tif" -o -name "*.tiff" -o -name "*.bmp" \) 2>/dev/null | head -20

echo ""
echo "=== 5. CHECK FOR SAMPLE SHEETS / METADATA ==="
find "$CONV" -maxdepth 3 \( -name "*.csv" -o -name "*.xml" -o -name "SampleSheet*" -o -name "*manifest*" -o -name "*Sample*" \) 2>/dev/null | head -30

echo ""
echo "=== 6. LOOK FOR CONTROL PROBE FILES ==="
find "$CONV" -maxdepth 3 \( -name "*control*" -o -name "*Control*" -o -name "*.ctl" -o -name "*staining*" -o -name "*QC*" \) -type f 2>/dev/null | head -20

echo ""
echo "=== 7. CHECK ConvSK SUBDIRECTORY ==="
ls -la "$CONV/ConvSK/" 2>/dev/null | head -30

echo ""
echo "=== 8. BROADER SEARCH: ALL CHIP DIRECTORIES ==="
# Check if chips have their own directories
for chip in $BAD_CHIPS $GOOD_CHIP; do
    echo "--- $chip ---"
    find "$CONV" -maxdepth 2 -type d -name "*${chip}*" 2>/dev/null
    # Also check just the directory listing for the chip
    ls -la "$CONV/ConvSK/${chip}/" 2>/dev/null | head -5 || echo "(no dir)"
done

echo ""
echo "=== 9. CHECK /staging/ALSU-analysis/ FOR IDAT STORAGE ==="
find /staging/ALSU-analysis/ -maxdepth 2 -type d \( -name "*idat*" -o -name "*IDAT*" -o -name "*raw*" -o -name "*idats*" \) 2>/dev/null | head -10

echo ""
echo "=== 10. FULL LISTING OF CONVERSION OUTPUT (first 3 levels) ==="
find "$CONV" -maxdepth 3 -type d 2>/dev/null | head -50

echo ""
echo "=== 11. SAMPLE SHEET: CHECK FOR CONTROL SAMPLES ==="
# The sample sheet should list control positions
find /staging/ALSU-analysis/ -maxdepth 3 -name "*.csv" -type f 2>/dev/null | while read f; do
    if grep -liq "control\|staining\|negative\|positive" "$f" 2>/dev/null; then
        echo "MATCH: $f"
        grep -i "control\|staining\|negative\|positive" "$f" | head -10
    fi
done

echo ""
echo "DONE"
