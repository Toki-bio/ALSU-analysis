#!/bin/bash
set -e

# =======================================================================
# V2 PBS + LD Recomputation
# Uses UZB_1kG_v2_merged (77,111 SNPs x 3,595 samples)
# DRAGEN server, March 2026
# =======================================================================

WORKDIR=~/v2/pbs_v2
MERGED=~/v2/global_pca/UZB_1kG_v2_merged
POPMAP=~/v2/pca/pop_mapping.txt
UZB_QC=~/v2/plink/UZB_v2_qc

mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "===== V2 PBS + LD RECOMPUTATION ====="
echo "Start: $(date)"
echo "Merged data: $MERGED"
echo ""

# -----------------------------------------------------------------------
# 1. Build PLINK cluster file from pop_mapping + merged FAM
# -----------------------------------------------------------------------
echo "=== Step 1: Population cluster file ==="

# pop_mapping.txt: SampleID<tab>POP (EUR/EAS/SAS/AFR/AMR/UZB)
# merged FAM: FID IID ...
# PLINK --within needs: FID IID CLUSTER

awk 'NR==FNR {pop[$1]=$2; next} {
  fid=$1; iid=$2;
  if (iid in pop) p=pop[iid];
  else if (fid in pop) p=pop[fid];
  else p="UZB";
  print fid, iid, p
}' "$POPMAP" "${MERGED}.fam" > clusters.txt

# Report population counts
echo "Population counts:"
awk '{print $3}' clusters.txt | sort | uniq -c | sort -rn
echo ""

N_TOTAL=$(wc -l < clusters.txt)
N_UZB=$(grep -c ' UZB$' clusters.txt || echo 0)
N_EUR=$(grep -c ' EUR$' clusters.txt || echo 0)
N_EAS=$(grep -c ' EAS$' clusters.txt || echo 0)
N_SAS=$(grep -c ' SAS$' clusters.txt || echo 0)
N_AFR=$(grep -c ' AFR$' clusters.txt || echo 0)
N_AMR=$(grep -c ' AMR$' clusters.txt || echo 0)
N_SNPS=$(wc -l < "${MERGED}.bim")

echo "Total: $N_TOTAL samples, $N_SNPS SNPs"
echo "UZB=$N_UZB EUR=$N_EUR EAS=$N_EAS SAS=$N_SAS AFR=$N_AFR AMR=$N_AMR"
echo ""

# -----------------------------------------------------------------------
# 2. Extract per-population sample lists
# -----------------------------------------------------------------------
echo "=== Step 2: Per-population sample lists ==="
for POP in UZB EUR EAS SAS AFR; do
    awk -v p="$POP" '$3==p {print $1, $2}' clusters.txt > keep_${POP}.txt
    echo "  $POP: $(wc -l < keep_${POP}.txt) samples"
done
echo ""

# -----------------------------------------------------------------------
# 3. Compute per-population allele frequencies
# -----------------------------------------------------------------------
echo "=== Step 3: Allele frequencies ==="
for POP in UZB EUR EAS SAS AFR; do
    plink --bfile "$MERGED" \
          --keep keep_${POP}.txt \
          --freq \
          --out freq_${POP} \
          --allow-no-sex --silent
    echo "  freq_${POP}.frq written ($(wc -l < freq_${POP}.frq) lines)"
done
echo ""

# -----------------------------------------------------------------------
# 4. Compute pairwise per-SNP FST (3 pairs for PBS + 2 extra)
# -----------------------------------------------------------------------
echo "=== Step 4: Pairwise per-SNP FST ==="

compute_fst() {
    POP1=$1; POP2=$2
    echo "  Computing FST: ${POP1} vs ${POP2}..."

    # Create merged keep list
    cat keep_${POP1}.txt keep_${POP2}.txt > keep_${POP1}_${POP2}.txt

    # Create within-cluster file
    awk -v p="$POP1" '{print $1, $2, p}' keep_${POP1}.txt > within_${POP1}_${POP2}.txt
    awk -v p="$POP2" '{print $1, $2, p}' keep_${POP2}.txt >> within_${POP1}_${POP2}.txt

    plink --bfile "$MERGED" \
          --keep keep_${POP1}_${POP2}.txt \
          --fst \
          --within within_${POP1}_${POP2}.txt \
          --out fst_${POP1}_${POP2} \
          --allow-no-sex --silent

    # Extract genome-wide FST from log
    WFST=$(grep "Weighted" fst_${POP1}_${POP2}.log | awk '{print $NF}')
    MFST=$(grep "Mean" fst_${POP1}_${POP2}.log | awk '{print $NF}')
    echo "    Weighted FST = $WFST, Mean FST = $MFST"
}

# PBS triangle
compute_fst UZB EUR
compute_fst UZB EAS
compute_fst EUR EAS

# Extra pairs for ΔAF context
compute_fst UZB SAS
compute_fst UZB AFR

echo ""

# -----------------------------------------------------------------------
# 5. PBS + ΔAF + Near-private (Python)
# -----------------------------------------------------------------------
echo "=== Step 5: PBS calculation ==="

cat > pbs_calc.py << 'PYEOF'
import sys, math, json, os
from collections import defaultdict

def read_frq(path):
    """Read PLINK .frq file -> dict[snp] = maf"""
    freqs = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.split()
            if len(parts) >= 5:
                snp = parts[1]
                maf = float(parts[4]) if parts[4] != 'NA' else None
                freqs[snp] = maf
    return freqs

def read_fst(path):
    """Read PLINK .fst file -> dict[snp] = fst"""
    fst = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.split()
            if len(parts) >= 5:
                snp = parts[1]
                val = float(parts[4]) if parts[4] != 'nan' else 0.0
                fst[snp] = val
    return fst

def fst_to_T(fst_val):
    """Convert FST to divergence time, capping FST at 0.999"""
    f = max(0.0, min(fst_val, 0.999))
    return -math.log(1.0 - f)

# Read frequencies
print("  Loading frequencies...")
freq_UZB = read_frq("freq_UZB.frq")
freq_EUR = read_frq("freq_EUR.frq")
freq_EAS = read_frq("freq_EAS.frq")
freq_SAS = read_frq("freq_SAS.frq")
freq_AFR = read_frq("freq_AFR.frq")

# Read per-SNP FST
print("  Loading per-SNP FST...")
fst_UZB_EUR = read_fst("fst_UZB_EUR.fst")
fst_UZB_EAS = read_fst("fst_UZB_EAS.fst")
fst_EUR_EAS = read_fst("fst_EUR_EAS.fst")
fst_UZB_SAS = read_fst("fst_UZB_SAS.fst")
fst_UZB_AFR = read_fst("fst_UZB_AFR.fst")

# Read BIM for chromosome/position info
print("  Loading BIM...")
bim = {}
bim_order = []
with open(os.path.expanduser("~/v2/global_pca/UZB_1kG_v2_merged.bim")) as f:
    for line in f:
        parts = line.split()
        chrom, snp, cm, bp = parts[0], parts[1], parts[2], int(parts[3])
        bim[snp] = (chrom, bp)
        bim_order.append(snp)

# Compute PBS for all SNPs present in all three FST files
print("  Computing PBS...")
common_snps = set(fst_UZB_EUR.keys()) & set(fst_UZB_EAS.keys()) & set(fst_EUR_EAS.keys())
common_snps = common_snps & set(freq_UZB.keys()) & set(freq_EUR.keys()) & set(freq_EAS.keys())
common_snps = common_snps & set(freq_SAS.keys()) & set(freq_AFR.keys())

results = []
for snp in common_snps:
    T_ue = fst_to_T(fst_UZB_EUR[snp])
    T_ua = fst_to_T(fst_UZB_EAS[snp])
    T_ea = fst_to_T(fst_EUR_EAS[snp])

    pbs_uzb = (T_ue + T_ua - T_ea) / 2.0
    pbs_eur = (T_ue + T_ea - T_ua) / 2.0
    pbs_eas = (T_ua + T_ea - T_ue) / 2.0

    maf_uzb = freq_UZB[snp] if freq_UZB[snp] is not None else 0
    maf_eur = freq_EUR[snp] if freq_EUR[snp] is not None else 0
    maf_eas = freq_EAS[snp] if freq_EAS[snp] is not None else 0
    maf_sas = freq_SAS[snp] if freq_SAS[snp] is not None else 0
    maf_afr = freq_AFR[snp] if freq_AFR[snp] is not None else 0

    # Delta AF: min across 4 comparisons
    delta_af = min(
        abs(maf_uzb - maf_eur),
        abs(maf_uzb - maf_eas),
        abs(maf_uzb - maf_sas),
        abs(maf_uzb - maf_afr)
    )

    # Near-private: UZB MAF >= 5%, all others <= 1%
    near_private = (maf_uzb >= 0.05 and
                    maf_eur <= 0.01 and maf_eas <= 0.01 and
                    maf_sas <= 0.01 and maf_afr <= 0.01)

    chrom, bp = bim.get(snp, ('?', 0))

    results.append({
        'snp': snp, 'chr': chrom, 'bp': bp,
        'pbs_uzb': pbs_uzb, 'pbs_eur': pbs_eur, 'pbs_eas': pbs_eas,
        'maf_uzb': maf_uzb, 'maf_eur': maf_eur, 'maf_eas': maf_eas,
        'maf_sas': maf_sas, 'maf_afr': maf_afr,
        'delta_af': delta_af, 'near_private': near_private,
        'tier1': pbs_uzb >= 0.3,
        'tier2': delta_af >= 0.3,
        'tier3': near_private
    })

print(f"  Total SNPs analyzed: {len(results)}")

# ---- Summary statistics ----
pbs_vals = [r['pbs_uzb'] for r in results]
pbs_vals.sort()
n = len(pbs_vals)

tier1 = [r for r in results if r['tier1']]
tier2 = [r for r in results if r['tier2']]
tier3 = [r for r in results if r['tier3']]
candidates = [r for r in results if r['tier1'] or r['tier2'] or r['tier3']]

import statistics
stats = {
    'n_snps': n,
    'mean': statistics.mean(pbs_vals),
    'median': statistics.median(pbs_vals),
    'stdev': statistics.stdev(pbs_vals),
    'min': min(pbs_vals),
    'max': max(pbs_vals),
    'p95': pbs_vals[int(0.95*n)],
    'p99': pbs_vals[int(0.99*n)],
    'p999': pbs_vals[int(0.999*n)],
    'n_pbs_ge_03': len(tier1),
    'n_pbs_ge_01': len([v for v in pbs_vals if v >= 0.1]),
    'n_negative': len([v for v in pbs_vals if v < 0]),
    'pct_negative': len([v for v in pbs_vals if v < 0]) / n * 100,
    'n_tier1': len(tier1),
    'n_tier2': len(tier2),
    'n_tier3': len(tier3),
    'n_candidates': len(candidates)
}

print(f"\n=== PBS SUMMARY ===")
for k, v in stats.items():
    if isinstance(v, float):
        print(f"  {k}: {v:.6f}")
    else:
        print(f"  {k}: {v}")

# ---- Save stats as JSON ----
with open("pbs_stats.json", "w") as f:
    json.dump(stats, f, indent=2)

# ---- Save all candidates ----
candidates.sort(key=lambda r: -r['pbs_uzb'])
with open("pbs_candidates.json", "w") as f:
    json.dump(candidates, f, indent=2)

# ---- Save ALL results for histogram ----
with open("pbs_all.tsv", "w") as f:
    f.write("SNP\tCHR\tBP\tPBS_UZB\tPBS_EUR\tPBS_EAS\tMAF_UZB\tMAF_EUR\tMAF_EAS\tMAF_SAS\tMAF_AFR\tDELTA_AF\tNEAR_PRIVATE\tTIER1\tTIER2\tTIER3\n")
    for r in sorted(results, key=lambda x: -x['pbs_uzb']):
        f.write(f"{r['snp']}\t{r['chr']}\t{r['bp']}\t{r['pbs_uzb']:.6f}\t{r['pbs_eur']:.6f}\t{r['pbs_eas']:.6f}\t{r['maf_uzb']:.4f}\t{r['maf_eur']:.4f}\t{r['maf_eas']:.4f}\t{r['maf_sas']:.4f}\t{r['maf_afr']:.4f}\t{r['delta_af']:.4f}\t{1 if r['near_private'] else 0}\t{1 if r['tier1'] else 0}\t{1 if r['tier2'] else 0}\t{1 if r['tier3'] else 0}\n")

# ---- Top 20 table ----
print(f"\n=== TOP 20 PBS HITS ===")
print(f"{'#':>3} {'Chr':>3} {'Position':>12} {'PBS_UZB':>10} {'MAF_UZB':>8} {'MAF_EUR':>8} {'MAF_EAS':>8} {'MAF_SAS':>8} {'MAF_AFR':>8} {'PBS_EUR':>10} {'PBS_EAS':>10}")
for i, r in enumerate(candidates[:20]):
    print(f"{i+1:3d} {r['chr']:>3} {r['bp']:>12d} {r['pbs_uzb']:10.4f} {r['maf_uzb']:8.4f} {r['maf_eur']:8.4f} {r['maf_eas']:8.4f} {r['maf_sas']:8.4f} {r['maf_afr']:8.4f} {r['pbs_eur']:10.4f} {r['pbs_eas']:10.4f}")

# ---- PBS histogram data (for JS) ----
import numpy as np
bins = np.arange(-1.5, 5.05, 0.05)
hist, edges = np.histogram(pbs_vals, bins=bins)
hist_data = [{"x": round(float(edges[i]), 3), "y": int(hist[i])} for i in range(len(hist)) if hist[i] > 0]
with open("pbs_histogram.json", "w") as f:
    json.dump(hist_data, f)

# ---- Pairwise FST summary ----
print(f"\n=== PAIRWISE FST (genome-wide weighted) ===")
for pair in ["UZB_EUR", "UZB_EAS", "UZB_SAS", "UZB_AFR", "EUR_EAS"]:
    logf = f"fst_{pair}.log"
    if os.path.exists(logf):
        with open(logf) as lf:
            for line in lf:
                if "Weighted" in line:
                    print(f"  {pair}: Weighted {line.split()[-1]}")
                elif "Mean" in line and "Fst" in line:
                    print(f"  {pair}: Mean {line.split()[-1]}")

# ---- Create PLINK-compatible clump input ----
# Map PBS to proxy p-values: higher PBS -> smaller p
max_pbs = max(r['pbs_uzb'] for r in candidates) if candidates else 1
with open("pbs_clump_pval.txt", "w") as f:
    f.write("SNP\tP\n")
    for r in candidates:
        p_proxy = 1.0 - (r['pbs_uzb'] / (max_pbs + 0.001))
        p_proxy = max(p_proxy, 1e-300)
        f.write(f"{r['snp']}\t{p_proxy:.15e}\n")

# Also write candidate SNP list
with open("pbs_candidate_snps.txt", "w") as f:
    for r in candidates:
        f.write(f"{r['snp']}\n")

print(f"\nFiles written: pbs_stats.json, pbs_candidates.json, pbs_all.tsv, pbs_histogram.json, pbs_clump_pval.txt")
print(f"Candidate SNPs: {len(candidates)}")
PYEOF

python3 pbs_calc.py
echo ""

# -----------------------------------------------------------------------
# 6. LD clumping of PBS candidates
# -----------------------------------------------------------------------
echo "=== Step 6: LD Clumping ==="

N_CAND=$(wc -l < pbs_candidate_snps.txt)
echo "  Input candidates: $N_CAND"

if [ "$N_CAND" -gt 0 ]; then
    # Use the LD-pruned UZB data for LD reference
    # (same SNP set as the PBS analysis)
    plink --bfile "$MERGED" \
          --keep keep_UZB.txt \
          --clump pbs_clump_pval.txt \
          --clump-snp-field SNP \
          --clump-field P \
          --clump-p1 1 \
          --clump-p2 1 \
          --clump-r2 0.5 \
          --clump-kb 1000 \
          --allow-no-sex --silent \
          --out pbs_clumped 2>&1 || echo "  Clumping completed with warnings"

    if [ -f pbs_clumped.clumped ]; then
        N_LOCI=$(awk 'NR>1 && NF>0 {count++} END {print count+0}' pbs_clumped.clumped)
        echo "  Independent loci after clumping: $N_LOCI"
    else
        echo "  No clumped output (check logs)"
        N_LOCI=0
    fi

    # Compute pairwise LD among all PBS candidates
    plink --bfile "$MERGED" \
          --keep keep_UZB.txt \
          --extract pbs_candidate_snps.txt \
          --r2 \
          --ld-window-kb 5000 \
          --ld-window 99999 \
          --ld-window-r2 0.1 \
          --allow-no-sex --silent \
          --out pbs_ld 2>&1 || echo "  LD calc completed with warnings"
else
    echo "  No candidates to clump"
fi
echo ""

# -----------------------------------------------------------------------
# 7. LD Decay (genome-wide on UZB_v2_qc)
# -----------------------------------------------------------------------
echo "=== Step 7: LD Decay ==="

# Randomly sample 3000 SNPs from UZB_v2_qc
echo "  Sampling 3000 SNPs from $(wc -l < ${UZB_QC}.bim) total..."
awk '{print $2}' ${UZB_QC}.bim | shuf -n 3000 --random-source=<(yes 42) > decay_snps.txt 2>/dev/null || \
awk 'BEGIN{srand(42)} {if(rand()<3000/NR) print $2}' ${UZB_QC}.bim | head -3000 > decay_snps.txt

echo "  Selected $(wc -l < decay_snps.txt) SNPs for LD decay"

# Compute pairwise r² for decay SNPs
plink --bfile "$UZB_QC" \
      --extract decay_snps.txt \
      --r2 \
      --ld-window-kb 2000 \
      --ld-window 99999 \
      --ld-window-r2 0.0 \
      --allow-no-sex --silent \
      --out ld_decay 2>&1 || echo "  LD decay calc completed with warnings"

# Bin into 50kb windows using awk
if [ -f ld_decay.ld ]; then
    echo "  Binning LD decay results..."
    awk 'NR>1 {
        dist = ($5 > $2) ? ($5 - $2) : ($2 - $5);
        if ($1 == $4) {  # same chromosome only
            bin = int(dist / 50000) * 50;
            sum[bin] += $7;
            count[bin]++;
        }
    } END {
        for (b in sum) {
            printf "%d\t%.6f\t%d\n", b, sum[b]/count[b], count[b]
        }
    }' ld_decay.ld | sort -n > ld_decay_binned.tsv
    echo "  LD decay bins: $(wc -l < ld_decay_binned.tsv)"
else
    echo "  No LD decay output file"
fi
echo ""

# -----------------------------------------------------------------------
# 8. Parse clumping results to JSON
# -----------------------------------------------------------------------
echo "=== Step 8: Compile JSON results ==="

cat > compile_results.py << 'PYEOF2'
import json, os

# ---- Clumping results ----
clumps = []
if os.path.exists("pbs_clumped.clumped"):
    # Load PBS candidates for cross-reference
    candidates = json.load(open("pbs_candidates.json"))
    cand_map = {c['snp']: c for c in candidates}

    with open("pbs_clumped.clumped") as f:
        header = f.readline()
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            chrom = parts[0]
            snp = parts[2]
            bp = int(parts[3]) if len(parts) > 3 else 0
            n_total = int(parts[5]) if len(parts) > 5 else 0

            # Parse proxy SNPs (SP2 column, last field)
            proxies = []
            if len(parts) > 11 and parts[11] != "NONE":
                sp2 = parts[11]
                # SP2 format: "snp1(1),snp2(1),...""
                for token in sp2.split(","):
                    proxy_snp = token.split("(")[0].strip()
                    if proxy_snp and proxy_snp != "NONE":
                        proxies.append(proxy_snp)

            c = cand_map.get(snp, {})
            clumps.append({
                "chr": chrom,
                "bp": bp,
                "snp": snp,
                "n_total": len(proxies),
                "pbs": c.get('pbs_uzb', 0),
                "maf": c.get('maf_uzb', 0),
                "proxies": proxies
            })

    clumps.sort(key=lambda x: -x['pbs'])

# ---- LD pairs ----
ld_pairs = []
if os.path.exists("pbs_ld.ld"):
    with open("pbs_ld.ld") as f:
        header = f.readline()
        for line in f:
            parts = line.split()
            if len(parts) >= 7:
                snp_a = parts[2]
                snp_b = parts[5]
                r2 = float(parts[6])
                bp_a = int(parts[1])
                bp_b = int(parts[4])
                kb = abs(bp_b - bp_a) / 1000.0
                ld_pairs.append({
                    "snp_a": snp_a, "snp_b": snp_b,
                    "kb": round(kb, 1), "r2": round(r2, 3)
                })
    ld_pairs.sort(key=lambda x: -x['r2'])

# ---- LD decay ----
decay_curve = []
if os.path.exists("ld_decay_binned.tsv"):
    with open("ld_decay_binned.tsv") as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                kb = int(parts[0])
                mean_r2 = float(parts[1])
                n = int(parts[2])
                if kb <= 2000:
                    decay_curve.append({"kb": kb, "mean_r2": round(mean_r2, 4), "n": n})
    decay_curve.sort(key=lambda x: x['kb'])

# ---- Chr distribution ----
from collections import Counter
chr_counts = Counter(c['chr'] for c in clumps)

# ---- Compile LD_DATA object ----
stats = json.load(open("pbs_stats.json"))

ld_data = {
    "n_input": len(json.load(open("pbs_candidates.json"))) if os.path.exists("pbs_candidates.json") else 0,
    "n_loci": len(clumps),
    "clumps": clumps,
    "decay_curve": decay_curve,
    "pbs_ld_pairs": ld_pairs[:300],  # top 300 pairs
    "pbs_ld_pairs_total": len(ld_pairs),
    "clumps_with_proxies": len([c for c in clumps if c['n_total'] > 0]),
    "chr_loci_counts": {str(i): chr_counts.get(str(i), 0) for i in range(1, 23)},
    "total_snps_in_clumps": len(clumps),
    "mean_r2_0kb": decay_curve[0]['mean_r2'] if decay_curve else 0,
    "mean_r2_500kb": next((d['mean_r2'] for d in decay_curve if d['kb'] == 500), 0)
}

with open("ld_data.json", "w") as f:
    json.dump(ld_data, f, indent=2)

print(f"  ld_data.json: {len(clumps)} loci, {len(ld_pairs)} LD pairs, {len(decay_curve)} decay bins")

# ---- Also save FST summary ----
fst_summary = {}
for pair in ["UZB_EUR", "UZB_EAS", "UZB_SAS", "UZB_AFR", "EUR_EAS"]:
    logf = f"fst_{pair}.log"
    if os.path.exists(logf):
        with open(logf) as lf:
            for line in lf:
                if "Weighted" in line:
                    fst_summary[f"{pair}_weighted"] = float(line.split()[-1])
                elif "Mean" in line and "Fst" in line:
                    fst_summary[f"{pair}_mean"] = float(line.split()[-1])

with open("fst_summary.json", "w") as f:
    json.dump(fst_summary, f, indent=2)

print(f"  fst_summary.json: {len(fst_summary)} entries")
print(f"\nAll results compiled successfully.")
PYEOF2

python3 compile_results.py
echo ""

# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
echo "===== V2 PBS + LD RECOMPUTATION COMPLETE ====="
echo "End: $(date)"
echo ""
echo "Output files in $WORKDIR:"
ls -lh *.json *.tsv 2>/dev/null
echo ""
echo "Key files to download:"
echo "  pbs_stats.json       - PBS summary statistics"
echo "  pbs_candidates.json  - All candidate SNPs with details"
echo "  pbs_histogram.json   - PBS distribution for charts"
echo "  ld_data.json         - LD clumping + decay + pairs (for step13)"
echo "  fst_summary.json     - Pairwise FST values"
echo "  pbs_all.tsv          - Complete per-SNP results"
