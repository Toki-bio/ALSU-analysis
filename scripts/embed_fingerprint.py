"""Embed fingerprinting plan data into investigation_data_v2.json."""
import json, re

d = json.load(open('data/investigation_data_v2.json'))
es_map = {s['iid']: s for s in d['enriched_samples']}
verdicts = {v['iid']: v for v in d['sample_verdicts']}
ia = d['identity_audit']
unexpected = ia['unexpected_pair_details']

# Build the pair list
def strip_suffix(iid):
    m = re.match(r'^(.+?)[dt]$', iid)
    return m.group(1) if m else iid

# Collect all pairs to verify
pairs_to_verify = []
seen = set()

# From unverified KEEP samples
unverified = [v for v in d['sample_verdicts'] if v['identity_status'] == 'unverified']
for v in unverified:
    iid = v['iid']
    for up in unexpected:
        if up['iid1'] == iid or up['iid2'] == iid:
            partner = up['iid2'] if up['iid1'] == iid else up['iid1']
            key = tuple(sorted([iid, partner]))
            if key not in seen:
                seen.add(key)
                sa = es_map.get(iid, {})
                sb = es_map.get(partner, {})
                pairs_to_verify.append({
                    'iid1': iid, 'iid2': partner,
                    'reason': 'unverified_pair',
                    'bc1': sa.get('barcode', ''), 'bc2': sb.get('barcode', ''),
                    'pos1': sa.get('position', ''), 'pos2': sb.get('position', ''),
                    'fmiss1': sa.get('fmiss'), 'fmiss2': sb.get('fmiss'),
                    'same_chip': up['same_chip'],
                    'action1': verdicts.get(iid, {}).get('action', '?'),
                    'action2': verdicts.get(partner, {}).get('action', '?')
                })

# KING new duplicates
king_new = [('07-04d', '04-50'), ('20-08', '08-123'), ('08-754', '08-749')]
for a, b in king_new:
    key = tuple(sorted([a, b]))
    if key not in seen:
        seen.add(key)
        sa = es_map.get(a, {})
        sb = es_map.get(b, {})
        pairs_to_verify.append({
            'iid1': a, 'iid2': b,
            'reason': 'king_new_duplicate',
            'bc1': sa.get('barcode', ''), 'bc2': sb.get('barcode', ''),
            'pos1': sa.get('position', ''), 'pos2': sb.get('position', ''),
            'fmiss1': sa.get('fmiss'), 'fmiss2': sb.get('fmiss'),
            'same_chip': sa.get('barcode') and sa.get('barcode') == sb.get('barcode'),
            'action1': verdicts.get(a, {}).get('action', '?'),
            'action2': verdicts.get(b, {}).get('action', '?')
        })

# Unique samples
all_samples = set()
for p in pairs_to_verify:
    all_samples.add(p['iid1'])
    all_samples.add(p['iid2'])

d['fingerprint_plan'] = {
    'pairs': pairs_to_verify,
    'total_samples': len(all_samples),
    'total_pairs': len(pairs_to_verify),
    'unverified_pairs': sum(1 for p in pairs_to_verify if p['reason'] == 'unverified_pair'),
    'king_new_pairs': sum(1 for p in pairs_to_verify if p['reason'] == 'king_new_duplicate'),
    'cost_low': len(all_samples) * 8,
    'cost_high': len(all_samples) * 25
}

with open('data/investigation_data_v2.json', 'w') as f:
    json.dump(d, f)

print(f"Fingerprint plan: {len(pairs_to_verify)} pairs, {len(all_samples)} samples")
print(f"Cost: ${len(all_samples)*8}–${len(all_samples)*25}")
