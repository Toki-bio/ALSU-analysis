import json, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

d = json.load(open('data/investigation_data_v2.json'))
es = d['enriched_samples']

catast_barcodes = {'208993030034', '208993030080', '208993030112', '208993030109'}

# All samples from catastrophic chips
catast_all = [s for s in es if s.get('barcode') in catast_barcodes]

print(f"Total from catastrophic chips: {len(catast_all)}")
print()

# Group by barcode
by_barcode = {}
for s in catast_all:
    bc = s.get('barcode')
    if bc not in by_barcode:
        by_barcode[bc] = []
    by_barcode[bc].append(s)

for bc in sorted(by_barcode.keys()):
    samples = by_barcode[bc]
    idxes = [s['idx'] for s in samples]
    print(f"{bc}: {len(samples)} samples, idx range {min(idxes)}-{max(idxes)}")

print()
print("=" * 60)
print("OUR 3 SAMPLES:")
target_iids = {'08-719': 1096, '08-825': 898, '08-701': 499}
for iid in ['08-701', '08-719', '08-825']:
    s = [s for s in es if s['iid'] == iid][0]
    bc = s.get('barcode')
    idx = s.get('idx')
    fmiss = s.get('fmiss')
    print(f"{iid}: barcode={bc}, idx={idx}, fmiss={fmiss:.4f}")

print()
print("Expected idx range for 208993030109 samples:")
s109 = [s for s in by_barcode.get('208993030109', []) if s['iid'] != '08-701']
if s109:
    idxes = [s['idx'] for s in s109]
    print(f"  {len(s109)+1} samples total, idx range {min(min(idxes), 499)}-{max(max(idxes), 499)}")
else:
    print("  Only 08-701")

print()
print("Expected idx range for 208993030112 samples:")
s112 = [s for s in by_barcode.get('208993030112', []) if s['iid'] != '08-825']
if s112:
    idxes = [s['idx'] for s in s112]
    print(f"  {len(s112)+1} samples total, idx range {min(min(idxes), 898)}-{max(max(idxes), 898)}")
else:
    print("  Only 08-825")

print()
print("Expected idx range for 208993030109 (all):")
s109_all = by_barcode.get('208993030109', [])
idxes = [s['idx'] for s in s109_all]
print(f"  {len(s109_all)} samples, idx range {min(idxes)}-{max(idxes)}")
print(f"  08-701 at idx={499} is OUTSIDE this range? {499 < min([x for x in idxes if x != 499]) or 499 > max([x for x in idxes if x != 499])}")

print()
print("===  HYPOTHESIS CHECK ===")
print("These 3 samples might be plotted but VISUALLY SEPARATE")
print("because their idx (sample order) is different from the main")
print("catastrophic cluster. The 1100-1200 dense block is where MOST")
print("catastrophic samples cluster, but these 3 are outliers.")
