#!/usr/bin/env bash
set -euo pipefail

echo "========================================================"
echo "  FULL PIPELINE AUDIT — $(date)"
echo "========================================================"

# ── 1. RAW / EARLY QC FILES ──────────────────────────────────
echo ""
echo "=== 1. EARLY QC FILES (winter2025) ==="
BASE=/staging/ALSU-analysis/winter2025
if [ -d "$BASE" ]; then
  for f in ConvSK_mind20 ConvSK_mind20_dedup ConvSK_mind20_dedup_snpqc; do
    bed="$BASE/$f.bed"
    bim="$BASE/$f.bim"
    fam="$BASE/$f.fam"
    if [ -f "$bim" ] && [ -f "$fam" ]; then
      snps=$(wc -l < "$bim")
      samps=$(wc -l < "$fam")
      echo "  $f: ${samps} samples, ${snps} SNPs"
    else
      echo "  $f: MISSING"
    fi
  done
else
  echo "  /staging/ALSU-analysis/winter2025/ NOT FOUND"
fi

# ── 2. V2 PLINK QC CHAIN ─────────────────────────────────────
echo ""
echo "=== 2. V2 PLINK QC CHAIN ==="
V2=/home/copilot/v2/plink
if [ -d "$V2" ]; then
  for f in UZB_v2_clean UZB_v2_qc UZB_v2_ldpruned; do
    bim="$V2/$f.bim"
    fam="$V2/$f.fam"
    if [ -f "$bim" ] && [ -f "$fam" ]; then
      snps=$(wc -l < "$bim")
      samps=$(wc -l < "$fam")
      echo "  $f: ${samps} samples, ${snps} SNPs"
    else
      echo "  $f: MISSING"
    fi
  done
  # Check for dot/period in sample IDs (the dot glitch)
  echo "  --- Dot glitch check (periods in IID) ---"
  dots=$(awk '$2 ~ /\./' "$V2/UZB_v2_qc.fam" | wc -l)
  echo "  Samples with '.' in IID: $dots"
  # Check for Cyrillic homoglyphs
  echo "  --- Cyrillic homoglyph check ---"
  cyrillic=$(grep -cP '[\x{0400}-\x{04FF}]' "$V2/UZB_v2_qc.fam" 2>/dev/null || echo 0)
  echo "  Samples with Cyrillic chars: $cyrillic"
else
  echo "  /home/copilot/v2/plink/ NOT FOUND"
fi

# ── 3. 1000G MERGE ───────────────────────────────────────────
echo ""
echo "=== 3. GLOBAL MERGE (1000G + UZB) ==="
GP=/home/copilot/v2/global_pca
if [ -d "$GP" ]; then
  for f in v2_common kg_common UZB_1kG_v2_merged; do
    bim="$GP/$f.bim"
    fam="$GP/$f.fam"
    if [ -f "$bim" ] && [ -f "$fam" ]; then
      snps=$(wc -l < "$bim")
      samps=$(wc -l < "$fam")
      echo "  $f: ${samps} samples, ${snps} SNPs"
    else
      echo "  $f: MISSING"
    fi
  done
  # Verify merged sample composition
  echo "  --- Merged FAM population breakdown ---"
  total=$(wc -l < "$GP/UZB_1kG_v2_merged.fam")
  uzb_like=$(awk '$2 !~ /^(HG|NA)/' "$GP/UZB_1kG_v2_merged.fam" | wc -l)
  kg_like=$(awk '$2 ~ /^(HG|NA)/' "$GP/UZB_1kG_v2_merged.fam" | wc -l)
  echo "  Total: $total  UZB-like: $uzb_like  1000G-like: $kg_like"
  # Check SNP overlap is exact
  echo "  --- SNP ID consistency ---"
  head -3 "$GP/UZB_1kG_v2_merged.bim"
  echo "  ..."
  tail -3 "$GP/UZB_1kG_v2_merged.bim"
else
  echo "  /home/copilot/v2/global_pca/ NOT FOUND"
fi

# ── 4. GLOBAL ADMIXTURE (single-run) ─────────────────────────
echo ""
echo "=== 4. GLOBAL ADMIXTURE (single-run K2-K8) ==="
GA=/home/copilot/v2/global_admixture
if [ -d "$GA" ]; then
  echo "  Input BED/BIM/FAM:"
  for ext in bed bim fam; do
    f="$GA/global_v2_admix.$ext"
    if [ -f "$f" ]; then
      sz=$(stat -c %s "$f")
      echo "    global_v2_admix.$ext: ${sz} bytes"
    else
      echo "    global_v2_admix.$ext: MISSING"
    fi
  done
  ga_snps=$(wc -l < "$GA/global_v2_admix.bim")
  ga_samps=$(wc -l < "$GA/global_v2_admix.fam")
  echo "  global_v2_admix: ${ga_samps} samples, ${ga_snps} SNPs"
  echo "  --- Are global_v2_admix and UZB_1kG_v2_merged the SAME data? ---"
  if diff <(sort "$GA/global_v2_admix.bim") <(sort "$GP/UZB_1kG_v2_merged.bim") > /dev/null 2>&1; then
    echo "  BIM files: IDENTICAL ✓"
  else
    echo "  BIM files: DIFFER ✗"
    diff <(wc -l < "$GA/global_v2_admix.bim") <(wc -l < "$GP/UZB_1kG_v2_merged.bim") || true
  fi
  if diff <(sort "$GA/global_v2_admix.fam") <(sort "$GP/UZB_1kG_v2_merged.fam") > /dev/null 2>&1; then
    echo "  FAM files: IDENTICAL ✓"
  else
    echo "  FAM files: DIFFER ✗"
  fi
  # CV errors from logs
  echo "  --- CV errors ---"
  for K in 2 3 4 5 6 7 8; do
    cv=$(grep -oP 'CV error \(K=\d+\): \K[\d.]+' "$GA/global_v2_admix_K${K}.log" 2>/dev/null || echo "N/A")
    ll=$(grep -oP 'Loglikelihood: \K[-\d.e+]+' "$GA/global_v2_admix_K${K}.log" 2>/dev/null | tail -1 || echo "N/A")
    echo "  K=$K: CV=$cv LogL=$ll"
  done
  # Q/P file sizes
  echo "  --- Q/P files ---"
  for K in 2 3 4 5 6 7 8; do
    q="$GA/global_v2_admix.${K}.Q"
    p="$GA/global_v2_admix.${K}.P"
    if [ -f "$q" ] && [ -f "$p" ]; then
      ql=$(wc -l < "$q")
      pl=$(wc -l < "$p")
      echo "  K=$K: Q=${ql} rows, P=${pl} rows"
    else
      echo "  K=$K: MISSING Q/P"
    fi
  done
else
  echo "  /home/copilot/v2/global_admixture/ NOT FOUND"
fi

# ── 5. EVANNO REPLICATE FILES ────────────────────────────────
echo ""
echo "=== 5. EVANNO REPLICATE INPUTS ==="
ER=/staging/ALSU-analysis/admixture_analysis
if [ -d "$ER" ]; then
  echo "  --- Master input files ---"
  for ext in bed bim fam; do
    f="$ER/UZB_for_admixture.$ext"
    if [ -f "$f" ]; then
      sz=$(stat -c %s "$f")
      echo "    UZB_for_admixture.$ext: ${sz} bytes"
    else
      echo "    UZB_for_admixture.$ext: MISSING"
    fi
  done
  ev_snps=$(wc -l < "$ER/UZB_for_admixture.bim" 2>/dev/null || echo "N/A")
  ev_samps=$(wc -l < "$ER/UZB_for_admixture.fam" 2>/dev/null || echo "N/A")
  echo "  UZB_for_admixture: ${ev_samps} samples, ${ev_snps} SNPs"

  echo ""
  echo "  --- CRITICAL: Does UZB_for_admixture match global_v2_admix? ---"
  if [ -f "$GA/global_v2_admix.bim" ] && [ -f "$ER/UZB_for_admixture.bim" ]; then
    if diff <(sort "$ER/UZB_for_admixture.bim") <(sort "$GA/global_v2_admix.bim") > /dev/null 2>&1; then
      echo "  BIM: IDENTICAL ✓ (same SNPs)"
    else
      echo "  BIM: DIFFER ✗ — THIS IS A PROBLEM"
      echo "    Evanno BIM: $(wc -l < "$ER/UZB_for_admixture.bim") SNPs"
      echo "    Global BIM: $(wc -l < "$GA/global_v2_admix.bim") SNPs"
    fi
    if diff <(sort "$ER/UZB_for_admixture.fam") <(sort "$GA/global_v2_admix.fam") > /dev/null 2>&1; then
      echo "  FAM: IDENTICAL ✓ (same samples)"
    else
      echo "  FAM: DIFFER ✗ — THIS IS A PROBLEM"
      echo "    Evanno FAM: $(wc -l < "$ER/UZB_for_admixture.fam") samples"
      echo "    Global FAM: $(wc -l < "$GA/global_v2_admix.fam") samples"
    fi
    # BED byte-level check (file size must match)
    ev_bed_sz=$(stat -c %s "$ER/UZB_for_admixture.bed")
    ga_bed_sz=$(stat -c %s "$GA/global_v2_admix.bed")
    if [ "$ev_bed_sz" = "$ga_bed_sz" ]; then
      echo "  BED size: IDENTICAL ✓ ($ev_bed_sz bytes)"
    else
      echo "  BED size: DIFFER ✗ (Evanno=$ev_bed_sz, Global=$ga_bed_sz) — THIS IS A PROBLEM"
    fi
  fi

  echo ""
  echo "  --- Symlink check in K7_rep4 (currently running) ---"
  rep4="$ER/evanno_runs/K7_rep4"
  if [ -d "$rep4" ]; then
    for ext in bed bim fam; do
      target=$(readlink -f "$rep4/UZB_for_admixture.$ext" 2>/dev/null || echo "NOT_A_LINK")
      expected="$ER/UZB_for_admixture.$ext"
      expected_resolved=$(readlink -f "$expected" 2>/dev/null || echo "N/A")
      if [ "$target" = "$expected_resolved" ]; then
        echo "    $ext: symlink → $target ✓"
      else
        echo "    $ext: MISMATCH → target=$target expected=$expected_resolved ✗"
      fi
    done
  else
    echo "    K7_rep4 dir NOT FOUND"
  fi

  echo ""
  echo "  --- Per-K replicate counts ---"
  for K in 2 3 4 5 6 7 8; do
    dirs=$(find "$ER/evanno_runs" -maxdepth 1 -type d -name "K${K}_rep*" | wc -l)
    q=$(find "$ER/evanno_runs" -maxdepth 2 -type f -path "*/K${K}_rep*/UZB_for_admixture.${K}.Q" | wc -l)
    p=$(find "$ER/evanno_runs" -maxdepth 2 -type f -path "*/K${K}_rep*/UZB_for_admixture.${K}.P" | wc -l)
    logs=$(find "$ER/evanno_runs" -maxdepth 2 -type f -path "*/K${K}_rep*/admixture_K${K}_rep*.log" | wc -l)
    echo "  K=$K: dirs=$dirs logs=$logs Q=$q P=$p"
  done

  echo ""
  echo "  --- Check Q file row counts match sample count ---"
  for K in 2 3 4 5 6 7 8; do
    for qf in "$ER"/evanno_runs/K${K}_rep*/UZB_for_admixture.${K}.Q; do
      [ -f "$qf" ] || continue
      ql=$(wc -l < "$qf")
      if [ "$ql" != "$ev_samps" ]; then
        echo "  MISMATCH: $qf has $ql rows, expected $ev_samps ✗"
      fi
    done
  done
  echo "  Q row count check done (no output = all match ✓)"

  echo ""
  echo "  --- Check for zero-size Q/P files ---"
  bad=0
  for f in "$ER"/evanno_runs/K*_rep*/UZB_for_admixture.*.Q "$ER"/evanno_runs/K*_rep*/UZB_for_admixture.*.P; do
    [ -f "$f" ] || continue
    if [ ! -s "$f" ]; then
      echo "  EMPTY FILE: $f ✗"
      bad=$((bad+1))
    fi
  done
  [ "$bad" -eq 0 ] && echo "  No empty Q/P files ✓"

  echo ""
  echo "  --- Completed replicate CV errors (spot-check K2, K5, K7) ---"
  for K in 2 5 7; do
    for logf in "$ER"/evanno_runs/K${K}_rep*/admixture_K${K}_rep*.log; do
      [ -f "$logf" ] || continue
      cv=$(grep -oP 'CV error \(K=\d+\): \K[\d.]+' "$logf" 2>/dev/null || echo "NO_CV")
      rep=$(basename "$(dirname "$logf")")
      echo "  $rep: CV=$cv"
    done
  done
else
  echo "  /staging/ALSU-analysis/admixture_analysis/ NOT FOUND"
fi

# ── 6. RUNNING PROCESS ───────────────────────────────────────
echo ""
echo "=== 6. ACTIVE ADMIXTURE PROCESS ==="
pgrep -af 'admixture --cv UZB_for_admixture' || echo "NO ACTIVE PROCESS"
echo ""
echo "=== 7. RUNNER STATUS ==="
pgrep -af 'tmp_evanno_complete_qp' || echo "No runner process"
tail -n 10 "$ER/evanno_runs/evanno_qp_nohup.out" 2>/dev/null || true

echo ""
echo "========================================================"
echo "  AUDIT COMPLETE — $(date)"
echo "========================================================"
