import json, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

d = json.load(open('data/investigation_data_v2.json'))
v_map = {s['iid']: s for s in d['sample_verdicts']}
e_map = {s['iid']: s for s in d['enriched_samples']}

for iid in ['08-719', '08-825', '08-701']:
    v = v_map.get(iid, {})
    e = e_map.get(iid, {})
    print(f'{iid}:')
    print(f'  action: {v.get("action", "?")}')
    print(f'  reasons: {v.get("reasons", "?")[:100]}')
    print(f'  fmiss: {e.get("fmiss", "?")}')
    print(f'  het_rate: {e.get("het_rate", "?")}')
    print(f'  barcode: {e.get("barcode", "?")}')
    print()

# Also check the catastrophic chips list
print("=" * 60)
print("CATASTROPHIC CHIPS (from dashboard):")
print("208993030034, 208993030080, 208993030112, 208993030109")
print()
print("Samples from these chips marked REMOVE:")
catast_removes = [s for s in d['sample_verdicts'] 
                  if s['action'] == 'REMOVE' and 
                  any(cat in s.get('reasons', '') for cat in 
                      ['208993030034', '208993030080', '208993030112', '208993030109'])]
print(f"Count: {len(catast_removes)}")
for s in catast_removes[:10]:
    e = e_map.get(s['iid'], {})
    print(f"  {s['iid']:14s} fmiss={e.get('fmiss',0):.4f} barcode={e.get('barcode','')}")
if len(catast_removes) > 10:
    print(f"  ... ({len(catast_removes)} total)")

# Now check: which contaminated samples are there?
print()
print("=" * 60)
print("CONTAMINATED samples (08-495, 08-25, 08-701):")
for iid in ['08-495', '08-25', '08-701']:
    if iid in v_map:
        v = v_map[iid]
        e = e_map.get(iid, {})
        print(f"  {iid}: action={v.get('action')} fmiss={e.get('fmiss',0):.4f} barcode={e.get('barcode','')}")
