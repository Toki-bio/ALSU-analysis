#!/bin/bash
set -euo pipefail
cd /staging/ALSU-analysis/spring2026/

echo "=== Generating removal list from connected components ==="

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

# Strategy: keep the sample with the lowest FID (earliest registered) per cluster
to_remove = []
with open('dup_clusters_summary.tsv', 'w') as cf:
    cf.write('cluster_id\tsize\tkept_fid\tkept_iid\tremoved\n')
    for i, cluster in enumerate(sorted(clusters, key=lambda c: min(fid_of[s] for s in c)), 1):
        members = sorted(cluster, key=lambda s: fid_of[s])
        keep = members[0]
        remove = members[1:]
        to_remove.extend(remove)
        removed_str = ','.join(f'{fid_of[s]}:{s}' for s in remove)
        cf.write(f'{i}\t{len(cluster)}\t{fid_of[keep]}\t{keep}\t{removed_str}\n')

with open('duplicates_pihat098.txt', 'w') as f:
    for iid in sorted(to_remove, key=lambda s: fid_of[s]):
        f.write(f'{fid_of[iid]}\t{iid}\n')

print(f'Pairs: {len(pairs)}')
print(f'Clusters: {len(clusters)}')
print(f'Samples involved: {len(adj)}')
print(f'To remove: {len(to_remove)}')
print(f'To keep: {len(clusters)}')
print(f'Expected final: 1148 - {len(to_remove)} = {1148 - len(to_remove)}')
PYEOF

echo ""
echo "=== Removal list ==="
wc -l duplicates_pihat098.txt
head -5 duplicates_pihat098.txt
echo "..."

echo ""
echo "=== Cluster summary ==="
wc -l dup_clusters_summary.tsv
head -5 dup_clusters_summary.tsv

echo ""
echo "=== Applying plink --remove ==="
plink --bfile ConvSK_mind20 \
  --remove duplicates_pihat098.txt \
  --make-bed \
  --out ConvSK_mind20_dedup

echo ""
echo "=== Verification ==="
wc -l ConvSK_mind20_dedup.fam
wc -l ConvSK_mind20_dedup.bim
ls -lh ConvSK_mind20_dedup.bed ConvSK_mind20_dedup.bim ConvSK_mind20_dedup.fam
