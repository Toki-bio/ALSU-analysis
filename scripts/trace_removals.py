"""Trace which REMOVE samples are still in the final 1,047 analysis dataset."""
import json, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

d = json.load(open('data/investigation_data_v2.json'))
verdicts = d['sample_verdicts']
enriched = {s['iid']: s for s in d['enriched_samples']}

# Build verdict lookup
v_map = {s['iid']: s for s in verdicts}

# The pipeline removal stages (from step docs):
# Step 1: 92 removed (F_MISS > 0.20) -> 1,155 remain
# Step 2: 57 removed (IBD dedup PI_HAT >= 0.98) -> 1,098 remain
# Step 6: 24 het outliers removed -> 1,074 remain
# Then 27 more with --mind 0.05 -> 1,047 remain

# All REMOVE verdicts
removes = [s for s in verdicts if s['action'] == 'REMOVE']
keeps = [s for s in verdicts if s['action'] == 'KEEP']

print(f"Total verdicts: {len(verdicts)}")
print(f"  KEEP:   {len(keeps)}")
print(f"  REMOVE: {len(removes)}")
print()

# Categorize REMOVE reasons
from collections import Counter
reasons = Counter()
for s in removes:
    for r in s.get('reasons', []):
        reasons[r] += 1
print("=== REMOVE REASONS ===")
for r, c in reasons.most_common():
    print(f"  {c:4d}  {r}")
print()

# Check which removes have high missingness (would have been caught by Step 1)
print("=== REMOVE samples by fmiss category ===")
step1_caught = []  # fmiss > 0.20
step6_caught = []  # fmiss 0.05-0.20 (het outliers or mind 0.05)
still_in = []      # fmiss < 0.05 (would survive all QC)
for s in removes:
    e = enriched.get(s['iid'], {})
    fmiss = e.get('fmiss', 0)
    if fmiss > 0.20:
        step1_caught.append(s)
    elif fmiss > 0.05:
        step6_caught.append(s)
    else:
        still_in.append(s)

print(f"  Caught by Step 1 (fmiss > 0.20): {len(step1_caught)}")
print(f"  Caught by Step 6 (fmiss 0.05-0.20): {len(step6_caught)}")
print(f"  LEAKED THROUGH (fmiss < 0.05): {len(still_in)}")
print()

# The IBD dedup (Step 2) removed 57 samples
# Check which of our REMOVE samples were in the dedup removal list
# The dedup list is the 57 from duplicates_pihat098.txt
# We can check: if a REMOVE sample has a duplicate pair partner who is KEEP,
# the REMOVE one should have been the one removed in Step 2

# Check Step 2 specifically
print("=== Step 2 IBD dedup overlap ===")
ibd_removes_iids = set()
for s in removes:
    rlist = ' '.join(s.get('reasons', []))
    if 'duplicate' in rlist.lower() or 'ibd' in rlist.lower():
        ibd_removes_iids.add(s['iid'])
print(f"  REMOVE verdicts mentioning duplicate/IBD: {len(ibd_removes_iids)}")
print()

# CRITICAL: the samples that LEAKED - detail each one
print("=" * 80)
print("=== SAMPLES MARKED REMOVE BUT LIKELY STILL IN FINAL DATASET ===")
print("=" * 80)
print(f"{'IID':14s} {'fmiss':>8s} {'het':>8s} {'F_het':>8s} {'chip':>12s} {'reason'}")
for s in sorted(still_in, key=lambda x: enriched.get(x['iid'], {}).get('fmiss', 0)):
    e = enriched.get(s['iid'], {})
    print(f"{s['iid']:14s} {e.get('fmiss',0):8.4f} {e.get('het_rate',0):8.4f} {e.get('f_het',0):8.4f} {str(e.get('barcode',''))[-8:]:>12s} {'; '.join(s.get('reasons',[]))}")

print()
print(f"\nTotal REMOVE samples likely still in final 1,047: {len(still_in)}")

# Also check: KING-only new duplicates (3 pairs found)
print()
print("=== NEW KING DUPLICATES (not caught by PLINK) ===")
king = d.get('king_verification', {})
new_dups = king.get('new_duplicates', [])
for p in new_dups:
    iid1, iid2 = p.get('iid1', ''), p.get('iid2', '')
    v1 = v_map.get(iid1, {})
    v2 = v_map.get(iid2, {})
    e1 = enriched.get(iid1, {})
    e2 = enriched.get(iid2, {})
    print(f"  {iid1} ({v1.get('action','?')}, fmiss={e1.get('fmiss',0):.4f}) <-> {iid2} ({v2.get('action','?')}, fmiss={e2.get('fmiss',0):.4f})")

# Identity audit: 30 KEEP with unverified identity
print()
print("=== IDENTITY AUDIT: KEEP but UNVERIFIED ===")
unverified = [s for s in keeps if s.get('identity_status') == 'unverified']
print(f"  Count: {len(unverified)}")
for s in unverified:
    e = enriched.get(s['iid'], {})
    print(f"  {s['iid']:14s} fmiss={e.get('fmiss',0):.4f} note={s.get('identity_note','')[:60]}")
