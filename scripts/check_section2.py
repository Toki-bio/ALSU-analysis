import json, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

d = json.load(open('data/investigation_data_v2.json'))
es = d['enriched_samples']

print(f'Total enriched_samples: {len(es)}')
print()

for iid in ['08-719', '08-825', '08-701']:
    found = [s for s in es if s['iid'] == iid]
    if found:
        s = found[0]
        fmiss = s.get('fmiss', 0)
        barcode = s.get('barcode', '')
        idx = s.get('idx')
        print(f'{iid}: idx={idx}, fmiss={fmiss:.4f}, barcode={barcode}')
        
        # Check if it would be in outBad (top panel, red color)
        in_top_panel = fmiss > 0.06
        is_bad_chip = barcode and barcode.startswith('20899303')
        print(f'  fmiss > 0.06? {in_top_panel}')
        print(f'  barcode starts with 20899303? {is_bad_chip}')
        print(f'  → Should be RED in Section 2 top panel? {in_top_panel and is_bad_chip}')
    else:
        print(f'{iid}: NOT FOUND in enriched_samples')
    print()
