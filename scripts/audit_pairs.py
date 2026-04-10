"""Audit ALL unexpected identity pairs — flag the identity ambiguity problem."""
import json, re

d = json.load(open('data/investigation_data_v2.json'))
es = {s['iid']: s for s in d['enriched_samples']}
pairs = d['pairs_098']
clusters = d['clusters']

def strip_suffix(iid):
    m = re.match(r'^(.+?)[dt]$', iid)
    return m.group(1) if m else iid

def get_suffix(iid):
    if iid.endswith('d'): return 'd'
    if iid.endswith('t'): return 't'
    return ''

print("=" * 100)
print("FULL AUDIT OF ALL 65 PI_HAT >= 0.98 PAIRS")
print("=" * 100)

# Categorize each pair
expected = []
unexpected_same_chip = []
unexpected_cross_chip = []

for p in pairs:
    i1, i2, pi = p['iid1'], p['iid2'], p['pihat']
    s1, s2 = es.get(i1, {}), es.get(i2, {})
    bc1, bc2 = s1.get('barcode', '?'), s2.get('barcode', '?')
    pos1, pos2 = s1.get('position', '?'), s2.get('position', '?')
    f1, f2 = s1.get('fmiss'), s2.get('fmiss')
    
    base1, base2 = strip_suffix(i1), strip_suffix(i2)
    is_expected = base1 == base2  # d/t of same person
    same_chip = bc1 == bc2 and bc1 != '?'
    
    entry = {
        'iid1': i1, 'iid2': i2, 'pihat': pi,
        'bc1': bc1, 'bc2': bc2, 'pos1': pos1, 'pos2': pos2,
        'fmiss1': f1, 'fmiss2': f2,
        'same_chip': same_chip, 'expected': is_expected,
        'suf1': get_suffix(i1), 'suf2': get_suffix(i2)
    }
    
    if is_expected:
        expected.append(entry)
    elif same_chip:
        unexpected_same_chip.append(entry)
    else:
        unexpected_cross_chip.append(entry)

print(f"\nExpected d/t pairs (base matches): {len(expected)}")
print(f"Unexpected SAME-CHIP pairs: {len(unexpected_same_chip)}")
print(f"Unexpected CROSS-CHIP pairs: {len(unexpected_cross_chip)}")

# EXPECTED: these are fine — d/t is supposed to match base
print("\n" + "=" * 100)
print("EXPECTED PAIRS (d/t matching base) — these are OK")
print("=" * 100)
for e in expected:
    print(f"  {e['iid1']:12s} <-> {e['iid2']:12s}  PI_HAT={e['pihat']:.4f}  chip={e['bc1'][-4:]}/{e['bc2'][-4:]}  pos={e['pos1']}/{e['pos2']}")

# UNEXPECTED SAME-CHIP: identity is AMBIGUOUS
print("\n" + "=" * 100)
print("UNEXPECTED SAME-CHIP — IDENTITY AMBIGUOUS (we don't know whose DNA this is)")
print("=" * 100)
for e in unexpected_same_chip:
    # Neither is a d/t of the other
    adjacent = False
    if e['pos1'] and e['pos2'] and e['pos1'] != '?' and e['pos2'] != '?':
        r1 = int(e['pos1'][1:3]) if len(e['pos1']) >= 3 else 0
        r2 = int(e['pos2'][1:3]) if len(e['pos2']) >= 3 else 0
        c1 = e['pos1'][3:] if len(e['pos1']) > 3 else ''
        c2 = e['pos2'][3:] if len(e['pos2']) > 3 else ''
        adjacent = abs(r1 - r2) <= 1
    
    problem = "IDENTITY UNKNOWN"
    if e['suf1'] or e['suf2']:
        # One is a d/t but of a DIFFERENT person
        dt_id = e['iid1'] if e['suf1'] else e['iid2']
        other_id = e['iid2'] if e['suf1'] else e['iid1']
        dt_base = strip_suffix(dt_id)
        problem = f"d/t '{dt_id}' (base={dt_base}) matches DIFFERENT person '{other_id}' — mislabeled d/t"
    
    adj_str = " ADJACENT" if adjacent else ""
    print(f"  {e['iid1']:12s} <-> {e['iid2']:12s}  PI_HAT={e['pihat']:.4f}  chip={e['bc1'][-4:]}  pos={e['pos1']}/{e['pos2']}{adj_str}")
    print(f"    fmiss: {e['fmiss1']:.4f} / {e['fmiss2']:.4f}  |  {problem}")

# UNEXPECTED CROSS-CHIP: even worse — DNA from different processing batches
print("\n" + "=" * 100)
print("UNEXPECTED CROSS-CHIP — WRONG DNA IN TUBE (both samples' identities questionable)")
print("=" * 100)
for e in unexpected_cross_chip:
    problem = "BOTH IDENTITIES QUESTIONABLE — tube swap"
    if e['suf1'] or e['suf2']:
        dt_id = e['iid1'] if e['suf1'] else e['iid2']
        other_id = e['iid2'] if e['suf1'] else e['iid1']
        dt_base = strip_suffix(dt_id)
        problem = f"d/t '{dt_id}' (base={dt_base}) matches DIFFERENT person '{other_id}'"
    
    print(f"  {e['iid1']:12s} <-> {e['iid2']:12s}  PI_HAT={e['pihat']:.4f}  chips={e['bc1'][-4:]}/{e['bc2'][-4:]}  pos={e['pos1']}/{e['pos2']}")
    print(f"    fmiss: {e['fmiss1']:.4f} / {e['fmiss2']:.4f}  |  {problem}")

# Now the critical question: for the "KEEP" samples from unexpected pairs,
# can we trust their identity?
print("\n" + "=" * 100)
print("VERDICT PROBLEM: samples we marked KEEP but whose identity is UNVERIFIED")
print("=" * 100)

# Load current verdicts
verdicts = {v['iid']: v for v in d['sample_verdicts']}

kept_but_ambiguous = []
for e in unexpected_same_chip + unexpected_cross_chip:
    for iid in [e['iid1'], e['iid2']]:
        v = verdicts.get(iid, {})
        if v.get('action') == 'KEEP':
            kept_but_ambiguous.append({
                'iid': iid,
                'partner': e['iid2'] if iid == e['iid1'] else e['iid1'],
                'same_chip': e['same_chip'],
                'pihat': e['pihat'],
                'fmiss': v.get('fmiss')
            })

print(f"\nSamples marked KEEP but involved in unexpected identity pairs: {len(kept_but_ambiguous)}")
for s in sorted(kept_but_ambiguous, key=lambda x: x['iid']):
    chip_str = "same-chip" if s['same_chip'] else "cross-chip"
    print(f"  {s['iid']:12s}  partner={s['partner']:12s}  PI_HAT={s['pihat']:.4f}  {chip_str}  fmiss={s['fmiss']:.4f}")
    print(f"    PROBLEM: We kept this sample but cannot confirm it IS person {s['iid']}. It might be {s['partner']}'s DNA.")
