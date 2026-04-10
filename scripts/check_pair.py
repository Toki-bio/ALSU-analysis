import json
d = json.load(open('data/investigation_data_v2.json'))
es = {s['iid']: s for s in d['enriched_samples']}

for p in d['pairs_098']:
    if '08-265' in (p['iid1'], p['iid2']) or '08-267' in (p['iid1'], p['iid2']):
        print(f"PAIR: {p}")

for iid in ['08-265', '08-267']:
    s = es.get(iid, {})
    print(f"{iid}: fmiss={s.get('fmiss')}, het={s.get('het_rate')}, bc={s.get('barcode')}, pos={s.get('position')}")

for ci, cl in enumerate(d['clusters']):
    if '08-265' in cl or '08-267' in cl:
        print(f"Cluster #{ci}: {cl}")
