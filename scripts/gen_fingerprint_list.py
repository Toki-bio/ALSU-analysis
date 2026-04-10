"""Generate the fingerprinting sample list:
- All 30 unverified KEEP samples + their pair partners
- New KING-only duplicates (07-04d/04-50, 20-08/08-123, 08-754/08-749)
Output: fingerprint_samples.txt (FID IID) and a summary."""
import json

d = json.load(open('data/investigation_data_v2.json'))
es_map = {s['iid']: s for s in d['enriched_samples']}

# 1. Unverified KEEP samples and their partners
verdicts = {v['iid']: v for v in d['sample_verdicts']}
unverified = [v for v in d['sample_verdicts'] if v['identity_status'] == 'unverified']

# Get their partners from identity_audit
ia = d['identity_audit']
unexpected = ia['unexpected_pair_details']

# Build partner map for unverified samples
fingerprint_set = set()
pair_info = []  # for the summary

for v in unverified:
    iid = v['iid']
    fingerprint_set.add(iid)
    # Find all partners
    for up in unexpected:
        if up['iid1'] == iid or up['iid2'] == iid:
            partner = up['iid2'] if up['iid1'] == iid else up['iid1']
            fingerprint_set.add(partner)
            pair_info.append((iid, partner, 'unverified_pair'))

# 2. New KING-only genuine duplicates
king_new = [
    ('07-04d', '04-50', 'king_new_duplicate'),
    ('20-08', '08-123', 'king_new_duplicate'),
    ('08-754', '08-749', 'king_new_duplicate'),
]
for a, b, reason in king_new:
    fingerprint_set.add(a)
    fingerprint_set.add(b)
    pair_info.append((a, b, reason))

# Sort and write
samples = sorted(fingerprint_set)

# Write FID IID format (FID = FID from .fam, which is the numeric index)
# We need the FID from enriched_samples
with open('data/fingerprint_samples.txt', 'w') as f:
    f.write("FID\tIID\tBarcode\tPosition\tReason\n")
    for iid in samples:
        s = es_map.get(iid, {})
        fid = s.get('fid', iid)
        bc = s.get('barcode', '?')
        pos = s.get('position', '?')
        # Determine reason
        reasons = []
        v = verdicts.get(iid, {})
        if v.get('identity_status') == 'unverified':
            reasons.append('KEEP_UNVERIFIED')
        if v.get('action') == 'REMOVE':
            reasons.append('REMOVE_partner')
        for a, b, r in king_new:
            if iid in (a, b):
                reasons.append('KING_new')
        if not reasons:
            reasons.append('partner_of_unverified')
        f.write(f"{fid}\t{iid}\t{bc}\t{pos}\t{'; '.join(reasons)}\n")

# Also write minimal FID IID for lab
with open('data/fingerprint_samples_lab.txt', 'w') as f:
    for iid in samples:
        s = es_map.get(iid, {})
        fid = s.get('fid', iid)
        f.write(f"{fid}\t{iid}\n")

print(f"=== FINGERPRINTING SAMPLE LIST ===")
print(f"Total unique samples: {len(samples)}")
print(f"  From unverified KEEP: {sum(1 for s in samples if verdicts.get(s,{}).get('identity_status')=='unverified')}")
print(f"  Partners of unverified: {sum(1 for s in samples if verdicts.get(s,{}).get('action')=='REMOVE' and s not in [x[0] for x in king_new]+[x[1] for x in king_new])}")
print(f"  KING new duplicates: {sum(1 for s in samples if s in {a for a,b,_ in king_new}|{b for a,b,_ in king_new})}")
print()

# Print grouped by pair
seen = set()
print("=== PAIRS TO VERIFY ===")
print(f"{'Sample A':14s} {'Sample B':14s} {'Reason':25s} {'Chip A':6s} {'Chip B':6s}")
print("-" * 70)
for a, b, reason in sorted(set(pair_info), key=lambda x: x[0]):
    key = tuple(sorted([a, b]))
    if key in seen:
        continue
    seen.add(key)
    sa = es_map.get(a, {})
    sb = es_map.get(b, {})
    bca = (sa.get('barcode','') or '')[-4:]
    bcb = (sb.get('barcode','') or '')[-4:]
    print(f"  {a:14s} {b:14s} {reason:25s} {bca:6s} {bcb:6s}")

print(f"\nTotal unique pairs: {len(seen)}")
print(f"\nFiles written:")
print(f"  data/fingerprint_samples.txt     (annotated, {len(samples)} samples)")
print(f"  data/fingerprint_samples_lab.txt  (FID+IID only, for lab submission)")
