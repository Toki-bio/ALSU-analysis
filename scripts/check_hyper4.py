import json

d = json.load(open('data/investigation_data_v2.json'))
targets = ['08-123', '20-05', '08-86', '07-04d']
es_map = {s['iid']: s for s in d['enriched_samples']}
pairs = d['pairs_098']
clusters = d['clusters']
het_mean = d['global_het_mean']
het_std = d['global_het_std']

for iid in targets:
    s = es_map.get(iid, {})
    het = s.get('het_rate')
    fmiss = s.get('fmiss')
    sex_f = s.get('sex_f')
    z = (het - het_mean) / het_std if het else None
    print(f"=== {iid} ===")
    print(f"  fmiss={fmiss}, het_rate={het}, z_het={z:.2f}" if z else f"  fmiss={fmiss}, het_rate={het}")
    print(f"  sex_f={sex_f}, snpsex={s.get('snpsex')}")
    print(f"  barcode={s.get('barcode')}, position={s.get('position')}")

    # Clusters
    for ci, cl in enumerate(clusters):
        if iid in cl:
            print(f"  Cluster #{ci}: {cl}")

    # All PI_HAT >= 0.98 pairs
    my_pairs = [p for p in pairs if p['iid1'] == iid or p['iid2'] == iid]
    print(f"  Pairs at PI_HAT>=0.98: {len(my_pairs)}")
    for p in my_pairs:
        other = p['iid2'] if p['iid1'] == iid else p['iid1']
        o = es_map.get(other, {})
        print(f"    -> {other} (pihat={p['pihat']:.4f}, other_fmiss={o.get('fmiss')}, other_het={o.get('het_rate')})")
    print()

# Summary verdict
print("=" * 60)
print("COMPARISON with contaminated samples:")
for iid in ['08-495', '08-25', '08-701']:
    s = es_map.get(iid, {})
    het = s.get('het_rate')
    z = (het - het_mean) / het_std if het else None
    n_pairs = len([p for p in pairs if p['iid1'] == iid or p['iid2'] == iid])
    print(f"  {iid}: het={het}, z={z:.1f}, fmiss={s.get('fmiss')}, sex_f={s.get('sex_f')}, pairs={n_pairs}")
