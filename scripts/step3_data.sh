#!/bin/bash
set -euo pipefail
cd /staging/ALSU-analysis/spring2026/
echo "=== STEP 3 DATA EXTRACTION ==="
echo "Date: $(date)"

# Check if Step 2 output exists
if [ ! -f ConvSK_mind20_dedup.fam ]; then
    echo "ConvSK_mind20_dedup not found. Running Step 2 dedup first..."
    # Generate removal list from IBD clusters
    python3 << 'PYEOF'
from collections import defaultdict
pairs = []
fid_of = {}
with open('ConvSK_mind20.genome') as f:
    header = f.readline()
    for line in f:
        parts = line.split()
        fid1, iid1 = int(parts[0]), parts[1]
        fid2, iid2 = int(parts[2]), parts[3]
        pairs.append((iid1, iid2))
        fid_of[iid1] = fid1
        fid_of[iid2] = fid2
adj = defaultdict(set)
for s1, s2 in pairs:
    adj[s1].add(s2)
    adj[s2].add(s1)
visited = set()
clusters = []
for node in sorted(adj.keys()):
    if node in visited:
        continue
    comp = set()
    q = [node]
    while q:
        n = q.pop(0)
        if n in visited:
            continue
        visited.add(n)
        comp.add(n)
        q.extend(nb for nb in adj[n] if nb not in visited)
    clusters.append(comp)
to_remove = []
for cluster in clusters:
    members = sorted(cluster, key=lambda s: fid_of[s])
    to_remove.extend(members[1:])
with open('duplicates_pihat098.txt', 'w') as f:
    for iid in sorted(to_remove, key=lambda s: fid_of[s]):
        f.write(f'{fid_of[iid]}\t{iid}\n')
print(f'Pairs: {len(pairs)}, Clusters: {len(clusters)}, Remove: {len(to_remove)}')
PYEOF
    plink --bfile ConvSK_mind20 \
      --remove duplicates_pihat098.txt \
      --make-bed --out ConvSK_mind20_dedup
    echo "Step 2 dedup complete."
fi

echo ""
echo "=== INPUT ==="
echo "Samples: $(wc -l < ConvSK_mind20_dedup.fam)"
echo "Variants: $(wc -l < ConvSK_mind20_dedup.bim)"

# Compute pre-QC statistics
echo ""
echo "=== PRE-QC STATS ==="
plink --bfile ConvSK_mind20_dedup --freq --hardy --missing --out step3_pre 2>&1 | tail -5

# MAF distribution
echo ""
echo "=== MAF_DIST ==="
awk 'NR>1 {
    maf = $5+0
    if (maf < 0.01) bin = 0
    else if (maf < 0.05) bin = 1
    else if (maf < 0.10) bin = 2
    else if (maf < 0.15) bin = 3
    else if (maf < 0.20) bin = 4
    else if (maf < 0.25) bin = 5
    else if (maf < 0.30) bin = 6
    else if (maf < 0.35) bin = 7
    else if (maf < 0.40) bin = 8
    else if (maf < 0.45) bin = 9
    else bin = 10
    counts[bin]++
}
END {
    split("0.00,0.01,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45", lo, ",")
    split("0.01,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50", hi, ",")
    for (i=0; i<=10; i++) printf "%s %s %d\n", lo[i+1], hi[i+1], counts[i]+0
}' step3_pre.frq

# HWE -log10(p) distribution
echo ""
echo "=== HWE_DIST ==="
awk 'NR>1 && $3=="ALL" {
    p = $9+0
    if (p <= 0) logp = 20
    else logp = -log(p)/log(10)
    if (logp < 1) bin = 0
    else if (logp < 2) bin = 1
    else if (logp < 3) bin = 2
    else if (logp < 4) bin = 3
    else if (logp < 5) bin = 4
    else if (logp < 6) bin = 5
    else if (logp < 10) bin = 6
    else bin = 7
    counts[bin]++
}
END {
    split("0,1,2,3,4,5,6,10", lo, ",")
    split("1,2,3,4,5,6,10,20", hi, ",")
    for (i=0; i<=7; i++) printf "%s %s %d\n", lo[i+1], hi[i+1], counts[i]+0
}' step3_pre.hwe

# Per-SNP missingness distribution
echo ""
echo "=== LMISS_DIST ==="
awk 'NR>1 {
    miss = $5+0
    if (miss < 0.005) bin = 0
    else if (miss < 0.01) bin = 1
    else if (miss < 0.02) bin = 2
    else if (miss < 0.05) bin = 3
    else if (miss < 0.10) bin = 4
    else if (miss < 0.15) bin = 5
    else if (miss < 0.20) bin = 6
    else bin = 7
    counts[bin]++
}
END {
    split("0.000,0.005,0.010,0.020,0.050,0.100,0.150,0.200", lo, ",")
    split("0.005,0.010,0.020,0.050,0.100,0.150,0.200,1.000", hi, ",")
    for (i=0; i<=7; i++) printf "%s %s %d\n", lo[i+1], hi[i+1], counts[i]+0
}' step3_pre.lmiss

# Run combined SNP QC
echo ""
echo "=== RUNNING SNP QC ==="
plink --bfile ConvSK_mind20_dedup \
  --maf 0.01 --hwe 1e-6 --geno 0.01 --chr 1-22 \
  --make-bed --out ConvSK_mind20_dedup_snpqc 2>&1 | tail -10

echo ""
echo "=== OUTPUT ==="
echo "Samples: $(wc -l < ConvSK_mind20_dedup_snpqc.fam)"
echo "Variants: $(wc -l < ConvSK_mind20_dedup_snpqc.bim)"

# Per-filter breakdown (individual runs)
echo ""
echo "=== PER_FILTER ==="
plink --bfile ConvSK_mind20_dedup --chr 1-22 --make-bed --out _tf_chr 2>/dev/null
echo "chr_only: $(wc -l < _tf_chr.bim)"
plink --bfile ConvSK_mind20_dedup --maf 0.01 --chr 1-22 --make-bed --out _tf_maf 2>/dev/null
echo "maf_only: $(wc -l < _tf_maf.bim)"
plink --bfile ConvSK_mind20_dedup --hwe 1e-6 --chr 1-22 --make-bed --out _tf_hwe 2>/dev/null
echo "hwe_only: $(wc -l < _tf_hwe.bim)"
plink --bfile ConvSK_mind20_dedup --geno 0.01 --chr 1-22 --make-bed --out _tf_geno 2>/dev/null
echo "geno_only: $(wc -l < _tf_geno.bim)"
rm -f _tf_*.bed _tf_*.bim _tf_*.fam _tf_*.log _tf_*.nosex _tf_*.hh

# Verification
echo ""
echo "=== VERIFY ==="
plink --bfile ConvSK_mind20_dedup_snpqc --freq --hardy --missing --out step3_verify 2>/dev/null
echo "Min_MAF: $(awk 'NR>1{print $5}' step3_verify.frq | sort -g | head -1)"
echo "Max_lmiss: $(awk 'NR>1{print $5}' step3_verify.lmiss | sort -g | tail -1)"
echo "Min_HWE_p: $(awk 'NR>1 && $3=="ALL"{print $9}' step3_verify.hwe | sort -g | head -1)"
echo "MAF_violations: $(awk 'NR>1 && $5+0 < 0.01' step3_verify.frq | wc -l)"
echo "GENO_violations: $(awk 'NR>1 && $5+0 > 0.01' step3_verify.lmiss | wc -l)"
echo "HWE_violations: $(awk 'NR>1 && $3=="ALL" && $9+0 < 1e-6' step3_verify.hwe | wc -l)"

# I/D alleles and duplicate positions in QC'd output
echo ""
echo "=== VCF_CLEANUP ==="
echo "ID_alleles: $(awk '($5=="I" || $5=="D" || $6=="I" || $6=="D")' ConvSK_mind20_dedup_snpqc.bim | wc -l)"
echo "Dup_positions: $(awk '{print $1,$4}' ConvSK_mind20_dedup_snpqc.bim | sort | uniq -d | wc -l)"

# Total autosomal in input
echo ""
echo "=== AUTOSOMAL ==="
echo "Input_autosomal: $(awk '$1>=1 && $1<=22' ConvSK_mind20_dedup.bim | wc -l)"
echo "Non_autosomal: $(awk '$1<1 || $1>22' ConvSK_mind20_dedup.bim | wc -l)"

# Per-chromosome variant counts in QC output
echo ""
echo "=== CHR_COUNTS ==="
awk '{counts[$1]++} END {for (c=1; c<=22; c++) printf "chr%d: %d\n", c, counts[c]+0}' ConvSK_mind20_dedup_snpqc.bim

# Clean up temp files
rm -f step3_pre.frq step3_pre.hwe step3_pre.lmiss step3_pre.log step3_pre.nosex step3_pre.hh
rm -f step3_verify.frq step3_verify.hwe step3_verify.lmiss step3_verify.log step3_verify.nosex step3_verify.hh

echo ""
echo "=== STEP 3 EXTRACTION COMPLETE ==="
