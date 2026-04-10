import json
d = json.load(open('data/investigation_data_v2.json'))
es = {s['iid']: s for s in d['enriched_samples']}

same_count = 0
diff_count = 0

for ds in d['dt_sex']:
    iid = ds['iid']
    base = ds['base']
    sdt = es.get(iid, {})
    sbase = es.get(base, {})
    bc_dt = sdt.get('barcode', '')
    bc_base = sbase.get('barcode', '')
    pos_dt = sdt.get('position', '')
    pos_base = sbase.get('position', '')
    same = bc_dt == bc_base and pos_dt == pos_base and bc_dt != ''
    if same:
        same_count += 1
    else:
        diff_count += 1
    b1 = bc_dt[-4:] if bc_dt else '?'
    b2 = bc_base[-4:] if bc_base else '?'
    flag = 'SAME' if same else 'DIFF'
    print(f'{iid:>12} {b1} {pos_dt:>6} | {base:>12} {b2} {pos_base:>6}  {flag}')

print(f'\nSame position: {same_count}, Different position: {diff_count}')
