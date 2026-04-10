#!/usr/bin/env python3
"""
ALSU Sample Investigation — Comprehensive Analysis
Investigates: duplicate IDs, high-relatedness pairs, missingness patterns
Uses ONLY real data from server extracts.
"""
import re
import json
from collections import defaultdict

DUMP = r"C:\work\ALSU-analysis\data\investigation_dump.txt"
IBD_FILE = r"C:\work\ALSU-analysis\data\ibd_results.txt"

def parse_sections(path):
    """Parse the multi-section dump file into dictionaries."""
    sections = {}
    current = None
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            m = re.match(r'^===(\w+.*?)===$', line)
            if m:
                current = m.group(1)
                sections[current] = []
                continue
            if current:
                sections[current].append(line)
    return sections

def parse_fam(lines):
    """Parse .fam lines -> list of (fid, iid, index_in_fam)"""
    samples = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 2:
            samples.append({'fid': parts[0], 'iid': parts[1], 'idx': int(parts[0])})
    return samples

def parse_imiss(lines):
    """Parse .imiss lines -> dict of iid -> F_MISS"""
    result = {}
    for line in lines:
        parts = line.split()
        if len(parts) >= 6 and parts[0] != 'FID':
            fid, iid = parts[0], parts[1]
            try:
                fmiss = float(parts[5])
                result[f"{fid}:{iid}"] = {'fid': fid, 'iid': iid, 'fmiss': fmiss}
            except ValueError:
                pass
    return result

def parse_high_pairs(path):
    """Parse the IBD results file for high PI_HAT pairs."""
    pairs = []
    in_section = False
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line == '=== HIGH_PIHAT_PAIRS ===':
                in_section = True
                continue
            if line.startswith('===') and in_section:
                break
            if in_section and line and not line.startswith('FID1'):
                parts = line.split()
                if len(parts) >= 8:
                    pairs.append({
                        'fid1': parts[0], 'iid1': parts[1],
                        'fid2': parts[2], 'iid2': parts[3],
                        'z0': float(parts[4]), 'z1': float(parts[5]),
                        'z2': float(parts[6]), 'pihat': float(parts[7])
                    })
    return pairs

def strip_suffix(iid):
    """Remove d/t suffix to get base sample ID."""
    m = re.match(r'^(.+?)[dt]$', iid)
    return m.group(1) if m else iid

def get_suffix(iid):
    """Get the suffix (d/t or empty)."""
    if iid.endswith('d'):
        return 'd'
    elif iid.endswith('t'):
        return 't'
    return ''

# ============================================================
# PARSE ALL DATA
# ============================================================
print("=" * 70)
print("ALSU SAMPLE INVESTIGATION REPORT")
print("=" * 70)

sections = parse_sections(DUMP)
fam_data = parse_fam(sections.get('FAM', []))
imiss_raw = parse_imiss(sections.get('IMISS_RAW', []))
imiss_mind20 = parse_imiss(sections.get('IMISS_MIND20', []))

# Parse removed samples
removed_miss20 = []
for line in sections.get('REMOVE_MISS20', []):
    parts = line.split()
    if len(parts) >= 2:
        removed_miss20.append({'fid': parts[0], 'iid': parts[1]})

# Parse dup clusters
dup_clusters = []
for line in sections.get('DUP_CLUSTERS', []):
    if line.strip() and not line.startswith('cluster'):
        dup_clusters.append(line)

# Parse all IBD pairs with PI_HAT >= 0.15
all_pairs = parse_high_pairs(IBD_FILE)

# Only 0.98+ pairs
pairs_098 = [p for p in all_pairs if p['pihat'] >= 0.98]

print(f"\nTotal samples in .fam: {len(fam_data)}")
print(f"Unique base IDs: {len(set(strip_suffix(s['iid']) for s in fam_data))}")
print(f"Samples with d suffix: {sum(1 for s in fam_data if s['iid'].endswith('d'))}")
print(f"Samples with t suffix: {sum(1 for s in fam_data if s['iid'].endswith('t'))}")
print(f"Samples removed for >20% missingness: {len(removed_miss20)}")
print(f"Samples in mind20 dataset: {len(imiss_mind20)}")
print(f"Total pairs with PI_HAT >= 0.15: {len(all_pairs)}")
print(f"Total pairs with PI_HAT >= 0.98: {len(pairs_098)}")

# ============================================================
# ANALYSIS 1: d/t suffix samples — which pair with their original?
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 1: DO d/t SUFFIX SAMPLES MATCH THEIR ORIGINALS?")
print("=" * 70)

# Build lookup: iid -> fid
iid_to_fid = {}
for s in fam_data:
    iid_to_fid[s['iid']] = s['fid']

# For each d/t sample, find if there's a 0.98+ pair with its base
dt_samples = [s for s in fam_data if get_suffix(s['iid']) in ('d', 't')]
print(f"\nTotal d/t suffix samples: {len(dt_samples)}")

# Build a lookup for 0.98+ pairs
pair_lookup_098 = defaultdict(dict)
for p in pairs_098:
    pair_lookup_098[p['iid1']][p['iid2']] = p['pihat']
    pair_lookup_098[p['iid2']][p['iid1']] = p['pihat']

# Also build lookup for ALL pairs (>= 0.15)
pair_lookup_all = defaultdict(dict)
for p in all_pairs:
    pair_lookup_all[p['iid1']][p['iid2']] = p['pihat']
    pair_lookup_all[p['iid2']][p['iid1']] = p['pihat']

# Track: which d/t samples were removed in step 1 (missingness)?
removed_iids = set(r['iid'] for r in removed_miss20)

dt_correct = []   # d/t matches its base at PI_HAT >= 0.98
dt_wrong = []     # d/t matches a DIFFERENT sample at PI_HAT >= 0.98
dt_no_match = []  # d/t doesn't match anything at PI_HAT >= 0.98
dt_removed = []   # d/t or base was removed before IBD

for dt in dt_samples:
    iid = dt['iid']
    base = strip_suffix(iid)
    suffix = get_suffix(iid)
    
    # Check if either was removed in missingness filter
    if iid in removed_iids or base in removed_iids:
        dt_removed.append({
            'iid': iid, 'base': base,
            'dt_removed': iid in removed_iids,
            'base_removed': base in removed_iids
        })
        continue
    
    # Check if base exists in fam
    if base not in iid_to_fid:
        dt_no_match.append({'iid': iid, 'base': base, 'reason': 'base not in fam'})
        continue
    
    # Check PI_HAT between dt and base
    pihat_with_base = pair_lookup_098.get(iid, {}).get(base, None)
    
    # What DOES this dt match at >= 0.98?
    matches = {k: v for k, v in pair_lookup_098.get(iid, {}).items()}
    
    if pihat_with_base is not None and pihat_with_base >= 0.98:
        dt_correct.append({'iid': iid, 'base': base, 'pihat': pihat_with_base, 'all_matches': matches})
    elif matches:
        dt_wrong.append({'iid': iid, 'base': base, 'pihat_base': pihat_with_base, 'matches': matches})
    else:
        # Check pair_lookup_all for lower PI_HAT with base
        pihat_low = pair_lookup_all.get(iid, {}).get(base, None)
        dt_no_match.append({'iid': iid, 'base': base, 'pihat_base_low': pihat_low})

print(f"\n--- d/t samples that CORRECTLY match their base (PI_HAT >= 0.98): {len(dt_correct)} ---")
for item in sorted(dt_correct, key=lambda x: x['iid']):
    extra = [f"{k}={v:.4f}" for k, v in item['all_matches'].items() if k != item['base']]
    extra_str = f"  (also matches: {', '.join(extra)})" if extra else ""
    print(f"  {item['iid']} ↔ {item['base']}: PI_HAT = {item['pihat']:.4f}{extra_str}")

print(f"\n--- d/t samples that match a DIFFERENT sample (WRONG identity): {len(dt_wrong)} ---")
for item in sorted(dt_wrong, key=lambda x: x['iid']):
    base_info = f"PI_HAT with {item['base']} = {item['pihat_base']:.4f}" if item['pihat_base'] else f"PI_HAT with {item['base']} < 0.15"
    matches_str = ", ".join(f"{k}={v:.4f}" for k, v in item['matches'].items())
    print(f"  {item['iid']} — Expected: {item['base']} ({base_info})")
    print(f"    Actual match(es): {matches_str}")

print(f"\n--- d/t samples with NO match at PI_HAT >= 0.98: {len(dt_no_match)} ---")
for item in sorted(dt_no_match, key=lambda x: x['iid']):
    if item.get('pihat_base_low'):
        print(f"  {item['iid']} — PI_HAT with {item['base']} = {item['pihat_base_low']:.4f} (low!)")
    elif item.get('reason'):
        print(f"  {item['iid']} — {item['reason']}")
    else:
        print(f"  {item['iid']} — PI_HAT with {item['base']} < 0.15 (completely unrelated!)")

print(f"\n--- d/t samples REMOVED in missingness filter (couldn't compare): {len(dt_removed)} ---")
for item in sorted(dt_removed, key=lambda x: x['iid']):
    who = []
    if item['dt_removed']:
        who.append(f"{item['iid']} removed")
    if item['base_removed']:
        who.append(f"{item['base']} removed")
    print(f"  {item['iid']} / {item['base']}: {', '.join(who)}")

# ============================================================
# ANALYSIS 2: Unexpected identical pairs (different base IDs)
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 2: UNEXPECTED IDENTICAL PAIRS (DIFFERENT BASE IDs, PI_HAT >= 0.98)")
print("=" * 70)

unexpected = []
expected_dt = []
for p in pairs_098:
    base1 = strip_suffix(p['iid1'])
    base2 = strip_suffix(p['iid2'])
    suf1 = get_suffix(p['iid1'])
    suf2 = get_suffix(p['iid2'])
    
    if base1 == base2:
        expected_dt.append(p)
    else:
        unexpected.append(p)

print(f"\nExpected (same base ID): {len(expected_dt)} pairs")
print(f"Unexpected (different base IDs): {len(unexpected)} pairs")

print(f"\n--- Unexpected identical pairs: ---")
for p in sorted(unexpected, key=lambda x: (x['iid1'], x['iid2'])):
    s1, s2 = get_suffix(p['iid1']), get_suffix(p['iid2'])
    tag1 = f" [{s1}]" if s1 else ""
    tag2 = f" [{s2}]" if s2 else ""
    print(f"  {p['fid1']:>4} {p['iid1']}{tag1} ↔ {p['fid2']:>4} {p['iid2']}{tag2}  PI_HAT={p['pihat']:.4f}")

# ============================================================
# ANALYSIS 3: Build connected components of identity clusters
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 3: IDENTITY CLUSTERS (connected components at PI_HAT >= 0.98)")
print("=" * 70)

from collections import deque

# Build adjacency from 0.98+ pairs
adj = defaultdict(set)
for p in pairs_098:
    adj[p['iid1']].add(p['iid2'])
    adj[p['iid2']].add(p['iid1'])

visited = set()
clusters = []
for node in adj:
    if node not in visited:
        cluster = []
        queue = deque([node])
        while queue:
            n = queue.popleft()
            if n in visited:
                continue
            visited.add(n)
            cluster.append(n)
            for neighbor in adj[n]:
                if neighbor not in visited:
                    queue.append(neighbor)
        clusters.append(sorted(cluster))

clusters.sort(key=lambda c: -len(c))

print(f"\nTotal identity clusters: {len(clusters)}")
print(f"Sizes: {[len(c) for c in clusters]}")

for i, cluster in enumerate(clusters):
    bases = set(strip_suffix(s) for s in cluster)
    originals = [s for s in cluster if not get_suffix(s)]
    dups = [s for s in cluster if get_suffix(s)]
    
    # Check missingness
    miss_info = []
    for s in cluster:
        key = f"{iid_to_fid.get(s, '?')}:{s}"
        fmiss = imiss_raw.get(key, {}).get('fmiss', None)
        miss_info.append(f"{s} (F_MISS={fmiss:.4f})" if fmiss is not None else f"{s} (filtered/unknown)")
    
    print(f"\n  Cluster {i+1} [{len(cluster)} members, {len(bases)} unique base IDs]:")
    print(f"    Members: {', '.join(miss_info)}")
    print(f"    Originals: {originals}")
    print(f"    d/t copies: {dups}")
    
    # Identify if any d/t maps to wrong original
    for d in dups:
        base = strip_suffix(d)
        if base not in [strip_suffix(o) for o in originals]:
            print(f"    ⚠ {d} has base '{base}' but '{base}' is NOT in this cluster!")
        elif base in [strip_suffix(o) for o in originals]:
            print(f"    ✓ {d} correctly clusters with its base {base}")

# ============================================================
# ANALYSIS 4: Missingness vs sample position / FID (index)
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 4: MISSINGNESS vs SAMPLE INDEX (FID)")
print("=" * 70)

# All 1247 samples in raw, sorted by FID (position)
fam_by_fid = sorted(fam_data, key=lambda x: x['idx'])

high_miss = []
low_miss = []
for s in fam_by_fid:
    key = f"{s['fid']}:{s['iid']}"
    info = imiss_raw.get(key, None)
    if info:
        entry = {**s, 'fmiss': info['fmiss']}
        if info['fmiss'] > 0.20:
            high_miss.append(entry)
        else:
            low_miss.append(entry)

print(f"\nSamples with F_MISS > 0.20: {len(high_miss)}")
print(f"Samples with F_MISS <= 0.20: {len(low_miss)}")

# Check if high-miss samples cluster at high FID numbers
high_fids = [s['idx'] for s in high_miss]
low_fids = [s['idx'] for s in low_miss]

print(f"\nHigh-missingness sample FID range: {min(high_fids)} - {max(high_fids)}")
print(f"Low-missingness sample FID range: {min(low_fids)} - {max(low_fids)}")

# Distribution across FID ranges
for threshold in [500, 700, 900, 1000, 1100, 1150, 1200]:
    above = [s for s in high_miss if s['idx'] >= threshold]
    total_above = sum(1 for s in fam_by_fid if s['idx'] >= threshold)
    pct = 100 * len(above) / len(high_miss) if high_miss else 0
    print(f"  FID >= {threshold}: {len(above)}/{len(high_miss)} high-miss samples ({pct:.1f}%), total samples >= {threshold}: {total_above}")

# Are high-miss samples d/t or original?
high_miss_dt = [s for s in high_miss if get_suffix(s['iid']) in ('d', 't')]
high_miss_orig = [s for s in high_miss if not get_suffix(s['iid'])]
print(f"\nOf high-missingness samples:")
print(f"  d/t suffix: {len(high_miss_dt)}")
print(f"  Original: {len(high_miss_orig)}")

# Show the high-miss samples sorted by FID
print(f"\n--- High-missingness samples (sorted by FID): ---")
for s in sorted(high_miss, key=lambda x: x['idx']):
    suf = get_suffix(s['iid'])
    tag = f" [{suf}]" if suf else ""
    print(f"  FID={s['idx']:>4}  IID={s['iid']}{tag}  F_MISS={s['fmiss']:.4f}")

# ============================================================
# ANALYSIS 5: First-degree relatives (PI_HAT 0.40 - 0.60) 
# that are NOT expected duplicates
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 5: POTENTIAL FIRST-DEGREE RELATIVES (0.40 <= PI_HAT < 0.60)")
print("=" * 70)

first_deg = [p for p in all_pairs if 0.40 <= p['pihat'] < 0.60 
             and strip_suffix(p['iid1']) != strip_suffix(p['iid2'])]
print(f"\nPairs in first-degree range (different base IDs): {len(first_deg)}")
for p in sorted(first_deg, key=lambda x: -x['pihat'])[:20]:
    print(f"  {p['iid1']} ↔ {p['iid2']}: PI_HAT={p['pihat']:.4f} Z0={p['z0']:.4f} Z1={p['z1']:.4f} Z2={p['z2']:.4f}")

# ============================================================
# ANALYSIS 6: Second-degree relatives (0.20 <= PI_HAT < 0.40)
# ============================================================
print("\n" + "=" * 70)
print("ANALYSIS 6: POTENTIAL 2ND-DEGREE RELATIVES (0.20 <= PI_HAT < 0.40)")
print("=" * 70)

second_deg = [p for p in all_pairs if 0.20 <= p['pihat'] < 0.40
              and strip_suffix(p['iid1']) != strip_suffix(p['iid2'])]
print(f"\n2nd-degree pairs (different base IDs): {len(second_deg)}")

# Group by how many connections each sample has
connections = defaultdict(int)
for p in second_deg:
    connections[p['iid1']] += 1
    connections[p['iid2']] += 1

if connections:
    top_connected = sorted(connections.items(), key=lambda x: -x[1])[:10]
    print(f"\nMost connected (2nd degree range):")
    for iid, cnt in top_connected:
        print(f"  {iid}: {cnt} connections")

# ============================================================
# SUMMARY OUTPUT FOR VISUALIZATION 
# ============================================================
# Save structured data for diagram generation
output = {
    'fam': [{'fid': s['fid'], 'iid': s['iid'], 'idx': s['idx']} for s in fam_data],
    'imiss_raw': {k: v['fmiss'] for k, v in imiss_raw.items()},
    'pairs_098': [{'iid1': p['iid1'], 'iid2': p['iid2'], 'pihat': p['pihat']} for p in pairs_098],
    'clusters': clusters,
    'high_miss': [{'iid': s['iid'], 'idx': s['idx'], 'fmiss': s['fmiss']} for s in high_miss],
    'unexpected_pairs': [{'iid1': p['iid1'], 'iid2': p['iid2'], 'pihat': p['pihat']} for p in unexpected]
}

with open(r"C:\work\ALSU-analysis\data\investigation_data.json", 'w') as f:
    json.dump(output, f, indent=2)

print("\n\nStructured data saved to data/investigation_data.json")
print("Ready for visualization.")
