import json

# Load ADMIXTURE data (1,047 samples with IDs and K2 proportions)
d = json.load(open('data/admix_v2_data.json', 'r', encoding='utf-8-sig'))
ids = d['ids']
k2 = d['k2']  # each entry is [q1, q2]

# Load ROH data
roh = {}
with open('data/UZB_v2_ROH.hom.indiv', 'r') as f:
    header = f.readline()  # skip header
    for line in f:
        parts = line.split()
        iid = parts[1]  # IID
        kb = float(parts[4])  # KB column
        roh[iid] = kb

# Autosomal genome length ~2,881,033 KB (hg38)
genome_kb = 2881033.0

# Build scatter data: Q1 (European-like) vs F_ROH
points = []
matched = 0
unmatched = 0
for i, sid in enumerate(ids):
    q1 = k2[i][0]  # first component
    if sid in roh:
        froh = roh[sid] / genome_kb
        points.append({'q': round(q1, 6), 'f': round(froh, 6), 'id': sid})
        matched += 1
    else:
        unmatched += 1

print(f'Matched: {matched}, Unmatched (no ROH): {unmatched}')

q_vals = [p['q'] for p in points]
f_vals = [p['f'] for p in points]
print(f'Q1 range: {min(q_vals):.4f} - {max(q_vals):.4f}')
print(f'F_ROH range: {min(f_vals):.6f} - {max(f_vals):.6f}')

# Count consanguineous (FROH > 0.0625)
consang = sum(1 for p in points if p['f'] > 0.0625)
print(f'Consanguineous (FROH > 0.0625): {consang}')

# Save as compact JSON for embedding
with open('data/admix_roh_scatter.json', 'w') as f:
    json.dump(points, f, separators=(',', ':'))
print(f'Saved {len(points)} points to data/admix_roh_scatter.json')

# Print size
import os
sz = os.path.getsize('data/admix_roh_scatter.json')
print(f'File size: {sz:,} bytes')
