"""Categorize all 1,247 samples with status, action, reason, and identity verification."""
import json, re

d = json.load(open('data/investigation_data_v2.json'))
ES = d['enriched_samples']
CS = d['chip_stats']
pairs = d['pairs_098']
clusters = d['clusters']
het_mean = d['global_het_mean']
het_std = d['global_het_std']

es_map = {s['iid']: s for s in ES}

# Catastrophic chips
catastrophic_barcodes = set(c['barcode'] for c in CS if c['mean_fmiss'] > 0.20)
# Contaminated
contaminated_ids = {'08-495', '08-25', '08-701'}

# Build identity pair graph
pair_partners = {}
for p in pairs:
    pair_partners.setdefault(p['iid1'], []).append((p['iid2'], p['pihat']))
    pair_partners.setdefault(p['iid2'], []).append((p['iid1'], p['pihat']))

# d/t samples
def strip_suffix(iid):
    m = re.match(r'^(.+?)[dt]$', iid)
    return m.group(1) if m else iid

def get_suffix(iid):
    if iid.endswith('d'): return 'd'
    if iid.endswith('t'): return 't'
    return ''

dt_ids = set(s['iid'] for s in ES if get_suffix(s['iid']) in ('d', 't'))

# ============================================================
# Classify all 65 pairs into expected vs unexpected
# ============================================================
expected_pairs = set()
unexpected_pairs = []
for p in pairs:
    base1, base2 = strip_suffix(p['iid1']), strip_suffix(p['iid2'])
    if base1 == base2:
        expected_pairs.add((p['iid1'], p['iid2']))
    else:
        s1 = es_map.get(p['iid1'], {})
        s2 = es_map.get(p['iid2'], {})
        same_chip = s1.get('barcode') and s1.get('barcode') == s2.get('barcode')
        unexpected_pairs.append({
            'iid1': p['iid1'], 'iid2': p['iid2'], 'pihat': p['pihat'],
            'same_chip': same_chip
        })

# Track all samples involved in unexpected pairs
# (both members, regardless of which one we keep/remove)
unexpected_pair_members = set()
# Map: sample -> list of unexpected partners
unexpected_partners = {}
for up in unexpected_pairs:
    unexpected_pair_members.add(up['iid1'])
    unexpected_pair_members.add(up['iid2'])
    unexpected_partners.setdefault(up['iid1'], []).append({
        'partner': up['iid2'], 'pihat': up['pihat'], 'same_chip': up['same_chip']
    })
    unexpected_partners.setdefault(up['iid2'], []).append({
        'partner': up['iid1'], 'pihat': up['pihat'], 'same_chip': up['same_chip']
    })

# For unexpected pairs, mark the higher-fmiss member for removal
unexpected_pair_remove = set()
for up in unexpected_pairs:
    s1 = es_map.get(up['iid1'], {})
    s2 = es_map.get(up['iid2'], {})
    f1 = s1.get('fmiss', 1)
    f2 = s2.get('fmiss', 1)
    # Remove the one with higher fmiss
    if f1 >= f2:
        unexpected_pair_remove.add(up['iid1'])
    else:
        unexpected_pair_remove.add(up['iid2'])

# Categorize every sample
categories = {}
for s in ES:
    iid = s['iid']
    fmiss = s.get('fmiss')
    het = s.get('het_rate')
    bc = s.get('barcode', '')
    reasons = []
    
    on_catastrophic = bc in catastrophic_barcodes
    is_contaminated = iid in contaminated_ids
    is_dt = iid in dt_ids
    is_ibd_remove = iid in unexpected_pair_remove
    high_fmiss = fmiss is not None and fmiss > 0.20
    
    if on_catastrophic:
        reasons.append(f'Catastrophic chip {bc[-4:]} (mean >20% missing)')
    if is_contaminated:
        reasons.append(f'Contaminated: two-person DNA mix (het={het:.3f}, z={((het-het_mean)/het_std):.1f})')
    if is_dt:
        reasons.append(f'Duplicate/triplicate suffix ({get_suffix(iid)}) — same idat as base {strip_suffix(iid)}')
    if is_ibd_remove:
        partners = pair_partners.get(iid, [])
        partner_str = ', '.join(f'{p[0]}' for p in partners[:3])
        # Check if same chip
        partner_s = es_map.get(partners[0][0], {}) if partners else {}
        same_chip = bc and partner_s.get('barcode') == bc
        if same_chip:
            reasons.append(f'PI_HAT=1.0 with {partner_str} (same chip, adjacent positions — same DNA in two wells or sample-sheet mapping error). One copy removed.')
        else:
            reasons.append(f'PI_HAT>=0.98 with {partner_str} (cross-chip — same DNA on different chips, tube-level swap). One copy removed.')
    if high_fmiss and not on_catastrophic:
        reasons.append(f'Individual high missingness (F_MISS={fmiss:.3f})')
    
    remove = on_catastrophic or is_contaminated or is_dt or is_ibd_remove or high_fmiss
    
    if remove:
        action = 'REMOVE'
        status = 'contaminated' if is_contaminated else 'chip_failure' if on_catastrophic else 'dt_artifact' if is_dt else 'ibd_duplicate' if is_ibd_remove else 'high_fmiss'
    else:
        action = 'KEEP'
        status = 'pass'
        reasons.append('Passes all QC thresholds')
    
    categories[iid] = {
        'action': action,
        'status': status,
        'reasons': reasons
    }

# Count summary
from collections import Counter
action_counts = Counter(c['action'] for c in categories.values())
status_counts = Counter(c['status'] for c in categories.values())

print(f"Total samples: {len(categories)}")
print(f"KEEP: {action_counts['KEEP']}, REMOVE: {action_counts['REMOVE']}")
print(f"\nRemoval breakdown:")
for status, count in sorted(status_counts.items(), key=lambda x: -x[1]):
    if status != 'pass':
        print(f"  {status}: {count}")

# Check for overlaps
remove_samples = [iid for iid, c in categories.items() if c['action'] == 'REMOVE']
multi_reason = [iid for iid, c in categories.items() if len(c['reasons']) > 1 and c['action'] == 'REMOVE']
print(f"\nSamples with multiple removal reasons: {len(multi_reason)}")
for iid in sorted(multi_reason)[:10]:
    print(f"  {iid}: {categories[iid]['reasons']}")

# Unique removals (deduplicated)
print(f"\nUnique samples to remove: {len(remove_samples)}")
print(f"Samples to keep: {action_counts['KEEP']}")

# Save as JSON for the dashboard
out = []
for s in ES:
    iid = s['iid']
    c = categories[iid]
    
    # Identity verification status
    if c['action'] == 'REMOVE':
        id_status = 'removed'
        id_note = ''
    elif iid in unexpected_pair_members:
        id_status = 'unverified'
        partners_info = unexpected_partners.get(iid, [])
        partner_names = [p['partner'] for p in partners_info]
        id_note = f"Identity unverifiable: PI_HAT≥0.98 with {', '.join(partner_names)}. Cannot confirm this is person {iid}."
    else:
        id_status = 'verified'
        id_note = ''
    
    out.append({
        'iid': iid,
        'fmiss': s.get('fmiss'),
        'het_rate': s.get('het_rate'),
        'sex_f': s.get('sex_f'),
        'barcode': s.get('barcode', ''),
        'position': s.get('position', ''),
        'action': c['action'],
        'status': c['status'],
        'reasons': ' | '.join(c['reasons']),
        'identity_status': id_status,
        'identity_note': id_note
    })

# Count identity stats
id_counts = {'verified': 0, 'unverified': 0, 'removed': 0}
for v in out:
    id_counts[v['identity_status']] = id_counts.get(v['identity_status'], 0) + 1

# Identity audit data for Section 16
identity_audit = {
    'expected_pairs': len(expected_pairs),
    'unexpected_same_chip': sum(1 for up in unexpected_pairs if up['same_chip']),
    'unexpected_cross_chip': sum(1 for up in unexpected_pairs if not up['same_chip']),
    'keep_unverified': id_counts['unverified'],
    'keep_verified': id_counts['verified'],
    'unexpected_pair_details': [{
        'iid1': up['iid1'], 'iid2': up['iid2'],
        'pihat': up['pihat'], 'same_chip': up['same_chip'],
        'fmiss1': es_map.get(up['iid1'], {}).get('fmiss'),
        'fmiss2': es_map.get(up['iid2'], {}).get('fmiss'),
        'bc1': es_map.get(up['iid1'], {}).get('barcode', ''),
        'bc2': es_map.get(up['iid2'], {}).get('barcode', ''),
        'pos1': es_map.get(up['iid1'], {}).get('position', ''),
        'pos2': es_map.get(up['iid2'], {}).get('position', '')
    } for up in unexpected_pairs]
}

# Embed in investigation_data_v2.json
d['sample_verdicts'] = out
d['verdict_summary'] = {
    'total': len(out),
    'keep': action_counts['KEEP'],
    'remove': action_counts['REMOVE'],
    'breakdown': dict(status_counts),
    'identity': id_counts
}
d['identity_audit'] = identity_audit
with open('data/investigation_data_v2.json', 'w') as f:
    json.dump(d, f)
print(f"\nSaved sample_verdicts ({len(out)} entries) into investigation_data_v2.json")
print(f"Identity: {id_counts['verified']} verified, {id_counts['unverified']} UNVERIFIED, {id_counts['removed']} removed")
