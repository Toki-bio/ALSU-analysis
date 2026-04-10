import json, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

d = json.load(open('data/investigation_data_v2.json'))
es = d['enriched_samples']

# Check if idx is sequential
idxes = sorted([s['idx'] for s in es])
print(f"idx values: min={min(idxes)}, max={max(idxes)}, count={len(idxes)}")
print(f"All unique? {len(idxes) == len(set(idxes))}")
print(f"Sequential? {idxes == list(range(1, len(idxes)+1))}")
print()

# Find 20-06m
target_iids = ['20-06m', '08-701', '08-825', '08-719']
print("Sample positions:")
for iid in target_iids:
    found = [s for s in es if s['iid'] == iid]
    if found:
        s = found[0]
        print(f"  {iid:14s}: idx={s.get('idx'):4d}, barcode={s.get('barcode', '?')[-4:]}, fmiss={s.get('fmiss', 0):.4f}")
    else:
        print(f"  {iid}: NOT FOUND")

print()
print("idx might be: a) sequential FID assignment, b) plate/well order, c) sequential ID number")
print("If samples processed out of order in GenomeStudio, idx wouldn't match processing sequence")
print()

# Check the first few samples to understand idx pattern
print("First 10 samples (by idx):")
sorted_by_idx = sorted(es, key=lambda x: x['idx'])
for s in sorted_by_idx[:10]:
    print(f"  idx={s['idx']:4d} iid={s['iid']:14s} fmiss={s.get('fmiss',0):.4f}")
