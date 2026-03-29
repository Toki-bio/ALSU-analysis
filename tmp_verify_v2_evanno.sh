#!/usr/bin/env bash
set -euo pipefail

echo "=============================================="
echo "V2 EVANNO LAUNCH VERIFICATION — $(date)"
echo "=============================================="

echo ""
echo "=== 1. IS ANYTHING RUNNING? ==="
pgrep -af 'admixture|run_evanno' || echo "NOTHING RUNNING!"

echo ""
echo "=== 2. WHAT PID IS ADMIXTURE, WHAT IS ITS CWD AND CMDLINE? ==="
APID=$(pgrep -f 'admixture --cv' | head -1)
if [[ -n "$APID" ]]; then
  echo "PID: $APID"
  echo "CMD: $(cat /proc/$APID/cmdline | tr '\0' ' ')"
  echo "CWD: $(readlink /proc/$APID/cwd)"
  echo "CPU: $(ps -p $APID -o %cpu= 2>/dev/null)%"
  echo "ELAPSED: $(ps -p $APID -o etime= 2>/dev/null)"
  echo "STATE: $(ps -p $APID -o stat= 2>/dev/null)"
  
  # What directory is it actually in?
  CWD=$(readlink /proc/$APID/cwd)
  
  echo ""
  echo "=== 3. FILES IN ACTIVE WORKING DIR ==="
  ls -la "$CWD/"
  
  echo ""
  echo "=== 4. SYMLINK TARGETS — WHERE DO .bed/.bim/.fam POINT? ==="
  for ext in bed bim fam; do
    f=$(ls "$CWD/"*.$ext 2>/dev/null | head -1)
    if [[ -L "$f" ]]; then
      TARGET=$(readlink -f "$f")
      echo "$ext -> $TARGET ($(stat -c%s "$TARGET") bytes)"
    elif [[ -f "$f" ]]; then
      echo "$ext: $f (regular file, $(stat -c%s "$f") bytes)"
    else
      echo "$ext: NOT FOUND!"
    fi
  done
  
  echo ""
  echo "=== 5. SAMPLE & SNP COUNTS IN ACTIVE INPUT ==="
  BIM=$(ls "$CWD/"*.bim 2>/dev/null | head -1)
  FAM=$(ls "$CWD/"*.fam 2>/dev/null | head -1)
  BED=$(ls "$CWD/"*.bed 2>/dev/null | head -1)
  
  BIMTARGET=$(readlink -f "$BIM" 2>/dev/null || echo "$BIM")
  FAMTARGET=$(readlink -f "$FAM" 2>/dev/null || echo "$FAM")
  
  NSNP=$(wc -l < "$BIMTARGET")
  NSAMP=$(wc -l < "$FAMTARGET")
  echo "SNPs (bim lines): $NSNP"
  echo "Samples (fam lines): $NSAMP"
  echo "First 3 fam lines:"
  head -3 "$FAMTARGET"
  echo "First 3 bim lines:"
  head -3 "$BIMTARGET"
  
  echo ""
  echo "=== 6. COMPARE WITH CANONICAL V2 FILES ==="
  echo "--- Global V2 canonical ---"
  GFAM=/home/copilot/v2/global_admixture/global_v2_admix.fam
  GBIM=/home/copilot/v2/global_admixture/global_v2_admix.bim
  GBED=/home/copilot/v2/global_admixture/global_v2_admix.bed
  echo "Samples: $(wc -l < $GFAM), SNPs: $(wc -l < $GBIM), BED: $(stat -c%s $GBED) bytes"
  
  echo "--- UZB V2 canonical ---"
  UFAM=/home/copilot/v2/admixture/UZB_v2_admix.fam
  UBIM=/home/copilot/v2/admixture/UZB_v2_admix.bim
  UBED=/home/copilot/v2/admixture/UZB_v2_admix.bed
  echo "Samples: $(wc -l < $UFAM), SNPs: $(wc -l < $UBIM), BED: $(stat -c%s $UBED) bytes"
  
  echo ""
  echo "--- Byte-identical check (symlink target vs canonical) ---"
  if cmp -s "$BIMTARGET" "$GBIM"; then
    echo "BIM: MATCHES global_v2_admix.bim ✓"
  elif cmp -s "$BIMTARGET" "$UBIM"; then
    echo "BIM: MATCHES UZB_v2_admix.bim ✓"
  else
    echo "BIM: MATCHES NEITHER! ✗✗✗"
  fi
  if cmp -s "$FAMTARGET" "$GFAM"; then
    echo "FAM: MATCHES global_v2_admix.fam ✓"
  elif cmp -s "$FAMTARGET" "$UFAM"; then
    echo "FAM: MATCHES UZB_v2_admix.fam ✓"
  else
    echo "FAM: MATCHES NEITHER! ✗✗✗"
  fi
  BEDTARGET=$(readlink -f "$BED" 2>/dev/null || echo "$BED")
  if cmp -s "$BEDTARGET" "$GBED"; then
    echo "BED: MATCHES global_v2_admix.bed ✓"
  elif cmp -s "$BEDTARGET" "$UBED"; then
    echo "BED: MATCHES UZB_v2_admix.bed ✓"
  else
    echo "BED: MATCHES NEITHER! ✗✗✗"
  fi
  
  # Also check it's NOT the old V1 file
  echo ""
  echo "--- NOT-V1 check ---"
  V1BIM=/staging/ALSU-analysis/admixture_analysis/UZB_for_admixture.bim
  V1FAM=/staging/ALSU-analysis/admixture_analysis/UZB_for_admixture.fam
  if [[ -f "$V1BIM" ]]; then
    if cmp -s "$BIMTARGET" "$V1BIM"; then
      echo "BIM: MATCHES OLD V1! ✗✗✗ ABORT!"
    else
      echo "BIM: Does NOT match V1 ✓ (good)"
    fi
  fi
  if [[ -f "$V1FAM" ]]; then
    if cmp -s "$FAMTARGET" "$V1FAM"; then
      echo "FAM: MATCHES OLD V1! ✗✗✗ ABORT!"
    else
      echo "FAM: Does NOT match V1 ✓ (good)"
    fi
  fi
  
else
  echo "No admixture process found!"
fi

echo ""
echo "=== 7. WHAT HAS THE RUNNER PRODUCED SO FAR? ==="
echo "--- Global evanno dir ---"
GDIR=/home/copilot/v2/global_admixture/evanno
if [[ -d "$GDIR" ]]; then
  for K in 2 3 4 5 6 7 8; do
    q=$(find "$GDIR" -maxdepth 2 -name "global_v2_admix.${K}.Q" -size +0c 2>/dev/null | wc -l)
    p=$(find "$GDIR" -maxdepth 2 -name "global_v2_admix.${K}.P" -size +0c 2>/dev/null | wc -l)
    echo "Global K${K}: Q=${q}/10 P=${p}/10"
  done
else
  echo "Global evanno dir does not exist!"
fi

echo "--- UZB evanno dir ---"
UDIR=/home/copilot/v2/admixture/evanno
if [[ -d "$UDIR" ]]; then
  for K in 2 3 4 5 6 7 8; do
    q=$(find "$UDIR" -maxdepth 2 -name "UZB_v2_admix.${K}.Q" -size +0c 2>/dev/null | wc -l)
    p=$(find "$UDIR" -maxdepth 2 -name "UZB_v2_admix.${K}.P" -size +0c 2>/dev/null | wc -l)
    echo "UZB K${K}: Q=${q}/10 P=${p}/10"
  done
else
  echo "UZB evanno dir not yet created (queued after global)"
fi

echo ""
echo "=== 8. LOG TAIL ==="
tail -20 /home/copilot/v2/global_admixture/evanno/evanno_global.log 2>/dev/null || echo "No global log yet"
tail -5 /home/copilot/v2/evanno_v2_nohup.out 2>/dev/null || echo "No nohup log"

echo ""
echo "=== 9. FIRST COMPLETED Q FILE SANITY ==="
FIRST_Q=$(find "$GDIR" -maxdepth 2 -name "*.Q" -size +0c 2>/dev/null | head -1)
if [[ -n "$FIRST_Q" ]]; then
  echo "File: $FIRST_Q"
  echo "Lines (should = samples): $(wc -l < "$FIRST_Q")"
  echo "Columns: $(head -1 "$FIRST_Q" | awk '{print NF}')"
  echo "First 3 lines:"
  head -3 "$FIRST_Q"
  echo "Sum of first row (should ≈ 1.0):"
  head -1 "$FIRST_Q" | awk '{s=0; for(i=1;i<=NF;i++) s+=$i; printf "%.6f\n", s}'
else
  echo "No Q files produced yet"
fi

echo ""
echo "=== 10. FIRST COMPLETED LOG - CV ERROR LINE ==="
FIRST_LOG=$(find "$GDIR" -maxdepth 2 -name "admixture_K*_rep*.log" -size +0c 2>/dev/null | sort | head -1)
if [[ -n "$FIRST_LOG" ]]; then
  echo "File: $FIRST_LOG"
  grep -i 'CV error\|loglikelihood' "$FIRST_LOG" | tail -3
else
  echo "No admixture logs yet"
fi

echo ""
echo "=== VERIFICATION COMPLETE ==="
