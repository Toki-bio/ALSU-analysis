"""Trace which REMOVE samples are still in the final 1,047 analysis dataset."""
import json, io, sys
from collections import Counter
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

d = json.load(open('data/investigation_data_v2.json'))
verdicts = d['sample_verdicts']
enriched = {s['iid']: s for s in d['enriched_samples']}
v_map = {s['iid']: s for s in verdicts}

removes = [s for s in verdicts if s['action'] == 'REMOVE']
keeps = [s for s in verdicts if s['action'] == 'KEEP']

print(f"Total verdicts: {len(verdicts)}")
print(f"  KEEP:   {len(keeps)}")
print(f"  REMOVE: {len(removes)}")
print()

# Categorize by what pipeline step would catch them
# Step 1: fmiss > 0.20 (92 removed)
# Step 2: IBD dedup - PI_HAT >= 0.98 (57 removed)
# Step 6: het outliers / --mind 0.05 (24+27=51 removed)

step1 = []  # fmiss > 0.20
step2 = []  # PI_HAT mentioned in reason, fmiss <= 0.20
step6 = []  # fmiss 0.05-0.20 (caught by --mind 0.05)
leaked = [] # fmiss < 0.05 and not IBD-removed

for s in removes:
    e = enriched.get(s['iid'], {})
    fmiss = e.get('fmiss', 0)
    reason = s.get('reasons', '')
    
    if fmiss > 0.20:
        step1.append(s)
    elif 'PI_HAT' in reason or 'pihat' in reason.lower():
        step2.append(s)
    elif fmiss > 0.05:
        step6.append(s)
    else:
        leaked.append(s)

print(f"=== REMOVE breakdown by pipeline stage ===")
print(f"  Caught by Step 1 (fmiss > 0.20):     {len(step1)}")
print(f"  Caught by Step 2 (IBD PI_HAT dedup):  {len(step2)}")
print(f"  Caught by Step 6 (fmiss 0.05-0.20):   {len(step6)}")
print(f"  LEAKED INTO FINAL DATASET (fmiss<0.05, not IBD): {len(leaked)}")
print()

# Check: Step 2 removes - are these the 57 from dedup?
print(f"=== Step 2 IBD removals detail ===")
for s in step2[:5]:
    e = enriched.get(s['iid'], {})
    print(f"  {s['iid']:14s} fmiss={e.get('fmiss',0):.4f} {s['reasons'][:80]}")
print(f"  ... ({len(step2)} total)")
print()

# CRITICAL: the leaked samples
print("=" * 80)
print("SAMPLES MARKED REMOVE BUT STILL IN FINAL ~1,047 DATASET")
print("=" * 80)
print(f"{'IID':14s} {'fmiss':>7s} {'het':>7s} {'f_het':>7s} {'chip[-4:]':>8s} reason")
print("-" * 120)
for s in sorted(leaked, key=lambda x: enriched.get(x['iid'], {}).get('fmiss', 0)):
    e = enriched.get(s['iid'], {})
    chip = str(e.get('barcode', ''))[-4:]
    print(f"{s['iid']:14s} {e.get('fmiss',0):7.4f} {e.get('het_rate',0):7.4f} {e.get('f_het',0):7.4f} {chip:>8s} {s['reasons'][:80]}")

print()
print(f"TOTAL LEAKED: {len(leaked)}")

# Classify leaked by reason type
print()
print("=== LEAKED by reason category ===")
cats = Counter()
for s in leaked:
    r = s['reasons']
    if 'Catastrophic' in r or 'catastrophic' in r:
        cats['Catastrophic chip'] += 1
    elif 'Duplicate' in r or 'duplicate' in r or 'suffix' in r:
        cats['Duplicate/suffix'] += 1
    elif 'contamina' in r.lower():
        cats['Contamination'] += 1
    elif 'PI_HAT' in r or 'pihat' in r.lower():
        cats['IBD (missed by step 2)'] += 1
    else:
        cats[f'Other: {r[:50]}'] += 1
for c, n in cats.most_common():
    print(f"  {n:4d}  {c}")

# KING new duplicates
print()
print("=== NEW KING DUPLICATES (3 pairs not caught by PLINK) ===")
king = d.get('king_verification', {})
new_dups = king.get('new_duplicates', [])
if new_dups:
    for p in new_dups:
        iid1, iid2 = p.get('iid1', ''), p.get('iid2', '')
        v1 = v_map.get(iid1, {})
        v2 = v_map.get(iid2, {})
        e1 = enriched.get(iid1, {})
        e2 = enriched.get(iid2, {})
        print(f"  {iid1} ({v1.get('action','?')}, fmiss={e1.get('fmiss',0):.4f}) <-> {iid2} ({v2.get('action','?')}, fmiss={e2.get('fmiss',0):.4f})")
else:
    print("  (not in JSON - check king_data.json)")

# Identity audit: 30 KEEP with unverified identity  
print()
print("=== UNVERIFIED KEEP (identity uncertain) ===")
unverified = [s for s in keeps if s.get('identity_status') == 'unverified']
print(f"  Count: {len(unverified)}")
for s in unverified[:10]:
    e = enriched.get(s['iid'], {})
    print(f"  {s['iid']:14s} fmiss={e.get('fmiss',0):.4f} {s.get('identity_note','')[:70]}")
if len(unverified) > 10:
    print(f"  ... ({len(unverified)} total)")

# Summary
print()
print("=" * 80)
print("PIPELINE IMPACT SUMMARY")
print("=" * 80)
print(f"Pipeline final dataset: 1,047 samples")
print(f"Investigation says REMOVE: 191 total")
print(f"  Already removed by pipeline:  {len(step1) + len(step2) + len(step6)}")
print(f"    Step 1 (missingness):       {len(step1)}")
print(f"    Step 2 (IBD dedup):         {len(step2)}")
print(f"    Step 6 (het outliers):      {len(step6)}")
print(f"  STILL IN DATASET:             {len(leaked)}")
print(f"  Plus KING new dups:           {len(new_dups)} pairs")
print(f"  Plus unverified identity:     {len(unverified)} KEEP samples")
