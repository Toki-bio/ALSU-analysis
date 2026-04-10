#!/usr/bin/env python3
"""
Analyze IBD clusters from ConvSK_mind20.genome and build deduplication strategy.
Uses graph-based connected components. Within each cluster, keep the sample with
the lowest FID (simplest, reproducible). Prints full analysis for step2.html.
"""
import re
from collections import defaultdict

# All 63 pairs from server data (IID1, IID2, PI_HAT)
pairs_raw = """
02-29 01-29 0.9999
02-29 01-29t 0.9998
08-131 08-131d 0.9998
01-29 01-29t 0.9999
02-39 06-30 0.9999
01-17 01-17d 0.9999
09-37 02-104 0.9999
03-154 03-155d 0.9999
03-155 03-156 0.9998
03-155 03-156d 0.9999
03-156 03-156d 0.9998
01-50 01-50d 0.9999
02-36 07-19 0.9999
02-45 08-107 0.9999
02-45 06-43d 1.0000
02-52 01-53 0.9998
02-59 02-45d 0.9999
02-59 02-45t 0.9988
04-07 06-34d 0.9998
04-07 06-15d 0.9999
04-08 07-01d 0.9997
08-493 08-744 0.9993
04-13 01-59d 0.9999
04-14 02-49d 1.0000
04-22 06-23d 0.9999
04-25 07-02d 0.9991
04-45 06-06d 0.9982
06-04 02-52d 1.0000
06-28 08-107d 0.9997
06-29 02-104d 0.9999
06-29 02-104t 0.9999
06-38 06-42d 0.9999
06-39 07-10d 0.9999
07-01 02-64d 1.0000
07-04 04-20d 1.0000
07-10 04-54d 1.0000
08-799 06-41d 0.9846
08-816 02-36d 0.9995
08-129 08-493d 1.0000
08-541 04-22d 0.9999
08-509 04-36 1.0000
08-509 04-40 0.9999
08-774 04-55d 1.0000
08-107 06-43d 1.0000
07-15 08-160 1.0000
07-16 02-90 0.9998
08-436 07-17 0.9998
08-179 08-194 0.9992
08-179 03-37 0.9992
08-194 03-37 0.9991
08-498 08-795d 1.0000
08-128 08-124 0.9994
06-34d 06-15d 0.9999
08-817 08-81d 0.9999
02-45d 02-45t 0.9989
02-104d 02-104t 0.9999
04-23 04-23d 0.9999
04-36 04-40 0.9999
08-770 08-770d 1.0000
01-18 09-76 1.0000
08-265 08-267 1.0000
08-181 08-45 1.0000
12-04 12-05 0.9998
""".strip()

# Also store FIDs from server data
fid_map_raw = """
11 02-29
52 01-29
338 01-29t
42 08-131
381 08-131d
71 02-39
206 06-30
90 01-17
291 01-17d
91 09-37
150 02-104
129 03-154
352 03-155d
130 03-155
131 03-156
353 03-156d
145 01-50
327 01-50d
152 02-36
270 07-19
154 02-45
237 08-107
325 06-43d
157 02-52
290 01-53
158 02-59
287 02-45d
347 02-45t
165 04-07
273 06-34d
276 06-15d
166 04-08
280 07-01d
168 08-493
422 08-744
169 04-13
284 01-59d
170 04-14
282 02-49d
176 04-22
266 06-23d
178 04-25
263 07-02d
187 04-45
262 06-06d
194 06-04
293 02-52d
204 06-28
380 08-107d
205 06-29
292 02-104d
342 02-104t
213 06-38
268 06-42d
214 06-39
278 07-10d
218 07-01
283 02-64d
220 07-04
297 04-20d
222 07-10
279 04-54d
225 08-799
300 06-41d
229 08-816
326 02-36d
231 08-129
415 08-493d
232 08-541
294 04-22d
233 08-509
296 04-36
369 04-40
234 08-774
281 04-55d
243 07-15
244 08-160
245 07-16
247 02-90
246 08-436
248 07-17
249 08-179
261 08-194
265 03-37
269 08-498
272 08-795d
271 08-128
275 08-124
277 08-817
306 08-81d
295 04-23
329 04-23d
335 08-770
519 08-770d
580 01-18
581 09-76
757 08-265
758 08-267
832 08-181
849 08-45
908 12-04
909 12-05
"""

# Parse FID map
fid_map = {}
for line in fid_map_raw.strip().split('\n'):
    parts = line.split()
    fid_map[parts[1]] = int(parts[0])

# Parse pairs
pairs = []
for line in pairs_raw.split('\n'):
    parts = line.split()
    pairs.append((parts[0], parts[1], float(parts[2])))

print(f"=== PAIR SUMMARY ===")
print(f"Total pairs: {len(pairs)}")
print(f"PI_HAT range: {min(p[2] for p in pairs):.4f} - {max(p[2] for p in pairs):.4f}")

# Build adjacency graph
adj = defaultdict(set)
for s1, s2, _ in pairs:
    adj[s1].add(s2)
    adj[s2].add(s1)

all_samples = set(adj.keys())
print(f"Unique samples in IBD pairs: {len(all_samples)}")

# Find connected components (clusters)
visited = set()
clusters = []

def bfs(start):
    queue = [start]
    component = set()
    while queue:
        node = queue.pop(0)
        if node in visited:
            continue
        visited.add(node)
        component.add(node)
        for neighbor in adj[node]:
            if neighbor not in visited:
                queue.append(neighbor)
    return component

for sample in sorted(all_samples):
    if sample not in visited:
        component = bfs(sample)
        clusters.append(sorted(component, key=lambda s: fid_map.get(s, 9999)))

clusters.sort(key=lambda c: fid_map.get(c[0], 9999))

print(f"\n=== CLUSTERS ({len(clusters)}) ===")
print(f"{'Cluster':>8}  {'Size':>5}  Members")
print("-" * 80)

to_remove = []
to_keep = []

for i, cluster in enumerate(clusters, 1):
    # Strategy: keep the sample with the lowest FID (earliest registered)
    members_with_fid = [(fid_map.get(s, 9999), s) for s in cluster]
    members_with_fid.sort()
    
    keep = members_with_fid[0][1]
    remove = [s for _, s in members_with_fid[1:]]
    
    to_keep.append(keep)
    to_remove.extend(remove)
    
    members_str = ', '.join(f"{s}(FID:{fid_map.get(s,'?')})" for _, s in members_with_fid)
    print(f"  C{i:>3}     {len(cluster):>3}    keep={keep}  remove=[{', '.join(remove)}]")

print(f"\n=== DEDUP SUMMARY ===")
print(f"Clusters: {len(clusters)}")
print(f"Samples involved: {len(all_samples)}")
print(f"To keep: {len(to_keep)}")
print(f"To remove: {len(to_remove)}")
print(f"Expected final: 1148 - {len(to_remove)} = {1148 - len(to_remove)}")

# Cluster size distribution
sizes = [len(c) for c in clusters]
size_dist = defaultdict(int)
for s in sizes:
    size_dist[s] += 1
print(f"\nCluster size distribution:")
for sz in sorted(size_dist):
    print(f"  Size {sz}: {size_dist[sz]} clusters")

# Print removal list in FID\tIID format
print(f"\n=== REMOVAL_LIST ({len(to_remove)} samples) ===")
for s in sorted(to_remove, key=lambda x: fid_map.get(x, 9999)):
    print(f"{fid_map.get(s, '?')}\t{s}")

# Check for d/t suffix patterns
d_suffix = [s for s in all_samples if s.endswith('d')]
t_suffix = [s for s in all_samples if s.endswith('t')]
no_suffix = [s for s in all_samples if not s.endswith('d') and not s.endswith('t')]
print(f"\n=== SUFFIX ANALYSIS ===")
print(f"Samples with 'd' suffix (explicit duplicates): {len(d_suffix)}")
print(f"Samples with 't' suffix (explicit triplicates): {len(t_suffix)}")
print(f"Samples without suffix: {len(no_suffix)}")

# Check: are all d/t suffix samples in the removal list?
d_in_remove = [s for s in d_suffix if s in to_remove]
d_kept = [s for s in d_suffix if s not in to_remove]
print(f"'d' suffix samples removed: {len(d_in_remove)}/{len(d_suffix)}")
if d_kept:
    print(f"  'd' suffix samples KEPT (lower FID than partner): {d_kept}")

t_in_remove = [s for s in t_suffix if s in to_remove]
t_kept = [s for s in t_suffix if s not in to_remove]
print(f"'t' suffix samples removed: {len(t_in_remove)}/{len(t_suffix)}")
if t_kept:
    print(f"  't' suffix samples KEPT: {t_kept}")

# Also identify non-obvious pairs (different IID roots)
print(f"\n=== NON-OBVIOUS PAIRS (different ID roots) ===")
for s1, s2, pihat in pairs:
    root1 = re.sub(r'[dt]$', '', s1)
    root2 = re.sub(r'[dt]$', '', s2)
    if root1 != root2:
        print(f"  {s1} (FID:{fid_map.get(s1,'?')}) -- {s2} (FID:{fid_map.get(s2,'?')})  PI_HAT={pihat}")
