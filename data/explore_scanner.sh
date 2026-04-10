#!/bin/bash
cd /staging/ALSU-analysis/Conversion/OUTPUT

echo "=== TOP LEVEL ==="
ls -la | head -40

echo "=== SUBDIRS ==="
find . -maxdepth 1 -type d | sort | head -30

echo "=== FILE TYPES ==="
find . -maxdepth 2 -type f | sed 's/.*\.//' | sort | uniq -c | sort -rn | head -20

echo "=== SAMPLE SHEETS ==="
find . -maxdepth 3 -name "*.csv" -o -name "*ample*heet*" -o -name "*manifest*" -o -name "*sample*" 2>/dev/null | head -20

echo "=== IDAT FILES ==="
find . -maxdepth 3 -name "*.idat" 2>/dev/null | head -10
echo "IDAT count:"
find . -maxdepth 3 -name "*.idat" 2>/dev/null | wc -l

echo "=== GTC FILES ==="
find . -maxdepth 3 -name "*.gtc" 2>/dev/null | head -10
echo "GTC count:"
find . -maxdepth 3 -name "*.gtc" 2>/dev/null | wc -l

echo "=== LOOK FOR SENTRIX PATTERNS ==="
find . -maxdepth 3 -type f -name "*Sentrix*" -o -name "*sentrix*" -o -name "*barcode*" 2>/dev/null | head -10

echo "=== CHECK FOR GENOMESTUDIO PROJECT ==="
find . -maxdepth 3 -name "*.bpm" -o -name "*.egt" -o -name "*.bsc" -o -name "*.xml" 2>/dev/null | head -20

echo "=== PARENT DIR ==="
ls -la /staging/ALSU-analysis/Conversion/ | head -20

echo "=== DONE ==="
