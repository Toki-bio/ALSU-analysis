"""Extract KING .kin0 data from the raw plink output and save as clean TSV.
Also compare KING pairs against PLINK PI_HAT pairs."""
import json, re

# Read the raw output
lines = open('data/king_output.txt', encoding='utf-16').read().splitlines()

# Find the KING duplicate pairs section
start = None
end = None
for i, line in enumerate(lines):
    if '=== KING DUPLICATE-LEVEL PAIRS ===' in line:
        start = i + 2  # skip the wc -l line, go to header
    if start and i > start and line.startswith('==='):
        end = i
        break
if not end:
    end = len(lines)

print(f"KING data: lines {start} to {end}")
king_lines = [l for l in lines[start:end] if l.strip()]
print(f"Total KING lines (incl header): {len(king_lines)}")

# Parse KING pairs
# Columns: #FID1 IID1 FID2 IID2 NSNP HETHET IBS0 KINSHIP
king_pairs = []
header = king_lines[0].split('\t') if king_lines else []
print(f"Header: {header}")

for line in king_lines[1:]:
    parts = line.split('\t')
    if len(parts) < 8:
        parts = line.split()
    if len(parts) < 8:
        continue
    king_pairs.append({
        'fid1': parts[0], 'iid1': parts[1],
        'fid2': parts[2], 'iid2': parts[3],
        'nsnp': int(parts[4]),
        'hethet': float(parts[5]),
        'ibs0': float(parts[6]),
        'kinship': float(parts[7])
    })

print(f"\nTotal KING pairs at kinship >= 0.354: {len(king_pairs)}")

# Load PLINK PI_HAT pairs
d = json.load(open('data/investigation_data_v2.json'))
plink_pairs = d['pairs_098']  # 65 pairs

# Load enriched samples for chip info
es_map = {s['iid']: s for s in d['enriched_samples']}
catastrophic_barcodes = set(c['barcode'] for c in d['chip_stats'] if c['mean_fmiss'] > 0.20)

# Make sets for comparison
# Normalize pair key (sorted IIDs)
def pair_key(a, b):
    return tuple(sorted([a, b]))

plink_set = {}
for p in plink_pairs:
    pk = pair_key(p['iid1'], p['iid2'])
    plink_set[pk] = p['pihat']

king_set = {}
for p in king_pairs:
    pk = pair_key(p['iid1'], p['iid2'])
    king_set[pk] = p

# Overlap analysis
both = set(plink_set.keys()) & set(king_set.keys())
plink_only = set(plink_set.keys()) - set(king_set.keys())
king_only = set(king_set.keys()) - set(plink_set.keys())

print(f"\n=== OVERLAP ANALYSIS ===")
print(f"In BOTH PLINK and KING: {len(both)}")
print(f"PLINK only (PI_HAT>=0.98 but KING<0.354): {len(plink_only)}")
print(f"KING only (kinship>=0.354 but PI_HAT<0.98): {len(king_only)}")

# Show PLINK-only pairs (most interesting - did KING disagree?)
if plink_only:
    print(f"\n=== PLINK-ONLY PAIRS (KING doesn't confirm) ===")
    for pk in sorted(plink_only):
        pihat = plink_set[pk]
        print(f"  {pk[0]:12s} <-> {pk[1]:12s}  PI_HAT={pihat:.4f}")

# Analyze KING-only pairs: how many involve catastrophic chips?
king_only_catastrophic = 0
king_only_clean = 0
for pk in king_only:
    kp = king_set[pk]
    bc1 = es_map.get(kp['iid1'], {}).get('barcode', '')
    bc2 = es_map.get(kp['iid2'], {}).get('barcode', '')
    on_bad = bc1 in catastrophic_barcodes or bc2 in catastrophic_barcodes
    if on_bad:
        king_only_catastrophic += 1
    else:
        king_only_clean += 1

print(f"\n=== KING-ONLY BREAKDOWN ===")
print(f"  Involve catastrophic chip: {king_only_catastrophic}")
print(f"  Both on clean chips: {king_only_clean}")

# Show sample of KING-only clean pairs (genuinely new findings)
if king_only_clean > 0:
    print(f"\n=== KING-ONLY CLEAN PAIRS (genuine new findings, up to 20) ===")
    count = 0
    for pk in sorted(king_only):
        kp = king_set[pk]
        bc1 = es_map.get(kp['iid1'], {}).get('barcode', '')
        bc2 = es_map.get(kp['iid2'], {}).get('barcode', '')
        if bc1 not in catastrophic_barcodes and bc2 not in catastrophic_barcodes:
            fm1 = es_map.get(kp['iid1'], {}).get('fmiss', -1)
            fm2 = es_map.get(kp['iid2'], {}).get('fmiss', -1)
            print(f"  {kp['iid1']:12s} <-> {kp['iid2']:12s}  KING={kp['kinship']:.4f}  NSNP={kp['nsnp']}  IBS0={kp['ibs0']:.6f}  fmiss={fm1:.4f}/{fm2:.4f}  chips={bc1[-4:]}/{bc2[-4:]}")
            count += 1
            if count >= 20:
                break

# For pairs in BOTH: show the comparison
print(f"\n=== CONCORDANCE: PLINK PI_HAT vs KING KINSHIP (first 20 of {len(both)}) ===")
for pk in sorted(both)[:20]:
    pihat = plink_set[pk]
    kp = king_set[pk]
    print(f"  {pk[0]:12s} <-> {pk[1]:12s}  PI_HAT={pihat:.4f}  KING={kp['kinship']:.4f}  IBS0={kp['ibs0']:.6f}")

# Save KING data to a clean JSON for the dashboard
king_data = {
    'pairs': [{
        'iid1': p['iid1'], 'iid2': p['iid2'],
        'kinship': p['kinship'], 'nsnp': p['nsnp'],
        'ibs0': p['ibs0'], 'hethet': p['hethet']
    } for p in king_pairs],
    'summary': {
        'total_king_pairs': len(king_pairs),
        'overlap_with_plink': len(both),
        'plink_only': len(plink_only),
        'king_only': len(king_only),
        'king_only_catastrophic': king_only_catastrophic,
        'king_only_clean': king_only_clean
    }
}

# Also create the overlay data: for the 55 unexpected PLINK pairs, add KING kinship
identity_audit = d.get('identity_audit', {})
unexpected_details = identity_audit.get('unexpected_pair_details', [])
overlay = []
for up in unexpected_details:
    pk = pair_key(up['iid1'], up['iid2'])
    kp = king_set.get(pk)
    overlay.append({
        'iid1': up['iid1'], 'iid2': up['iid2'],
        'pihat': plink_set.get(pk, None),
        'king_kinship': kp['kinship'] if kp else None,
        'king_ibs0': kp['ibs0'] if kp else None,
        'king_nsnp': kp['nsnp'] if kp else None,
        'same_chip': up['same_chip'],
        'confirmed_by_king': kp is not None
    })

confirmed = sum(1 for o in overlay if o['confirmed_by_king'])
not_confirmed = sum(1 for o in overlay if not o['confirmed_by_king'])
print(f"\n=== UNEXPECTED PAIR KING CONFIRMATION ===")
print(f"  55 unexpected PLINK pairs:")
print(f"  Confirmed by KING (kinship >= 0.354): {confirmed}")
print(f"  NOT confirmed by KING: {not_confirmed}")
if not_confirmed:
    print(f"  Pairs NOT confirmed:")
    for o in overlay:
        if not o['confirmed_by_king']:
            print(f"    {o['iid1']:12s} <-> {o['iid2']:12s}  PI_HAT={o['pihat']:.4f}")

# Save to JSON
king_data['unexpected_overlay'] = overlay
king_data['unexpected_confirmed'] = confirmed
king_data['unexpected_not_confirmed'] = not_confirmed

with open('data/king_data.json', 'w') as f:
    json.dump(king_data, f)
print(f"\nSaved king_data.json ({len(king_pairs)} pairs)")
