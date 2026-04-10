import json, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

d = json.load(open('data/investigation_data_v2.json'))
es = d['enriched_samples']

# Find any sample from catastrophic chips
catast_chip_ids = ['208993030034', '208993030080', '208993030112', '208993030109']
catast_samples = [s for s in es if s.get('barcode') in catast_chip_ids][:5]

print(f"Total samples from catastrophic chips: {len([s for s in es if s.get('barcode') in catast_chip_ids])}")
print()
print("Sample structure (from catastrophic chip):")
if catast_samples:
    s = catast_samples[0]
    for k, v in list(s.items())[:15]:
        print(f"  {k}: {v}")
print()

# Now check our 3 specific samples
print("=" * 60)
print("THREE SPECIFIC SAMPLES:")
for iid in ['08-719', '08-825', '08-701']:
    found = [s for s in es if s['iid'] == iid]
    if found:
        s = found[0]
        print(f"\n{iid}:")
        for k in ['iid', 'idx', 'fmiss', 'het_rate', 'barcode', 'position']:
            print(f"  {k}: {repr(s.get(k))}")
    else:
        print(f"\n{iid}: NOT FOUND")

# Check the types
print()
print("=" * 60)
print("DATA TYPES:")
if catast_samples:
    s = catast_samples[0]
    print(f"barcode type: {type(s.get('barcode'))} value={repr(s.get('barcode'))}")
    print(f"barcode.startswith check: {s.get('barcode').startswith('20899303') if s.get('barcode') else 'N/A'}")

s = [s for s in es if s['iid'] == '08-719'][0]
print(f"\n08-719 barcode type: {type(s.get('barcode'))} value={repr(s.get('barcode'))}")
print(f"08-719 barcode==208993030109? {s.get('barcode') == '208993030109'}")
print(f"08-719 barcode.startswith('20899303')? {s.get('barcode', '').startswith('20899303')}")
