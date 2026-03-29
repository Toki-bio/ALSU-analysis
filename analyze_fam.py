#!/usr/bin/env python3
"""Analyze V2 global FAM file to understand population composition."""

fam_lines = open("data/global/global_v2_admix.fam").read().strip().split('\n')
print(f"Total samples: {len(fam_lines)}")

iids = [line.split()[1] for line in fam_lines]

# 1000G IIDs start with HG or NA; Uzbek IIDs are numeric-prefix
ref_ids = [x for x in iids if x.startswith('HG') or x.startswith('NA')]
uzb_ids = [x for x in iids if not (x.startswith('HG') or x.startswith('NA'))]
print(f"1000G reference: {len(ref_ids)}")
print(f"Uzbek: {len(uzb_ids)}")

# The pop_mapping file has superpop — let's check what file the server used
# Since we can't SSH, let's try to use 1000G integrated_call_samples_v3 known IDs
# Actually, the population can be reconstructed from standard 1000G Phase 3 sample lists
# But we don't have that locally.

# Alternative: we know from V1 which 1000G populations were used.
# For V2 the set might be different. Let me just check unique IID prefixes.
print(f"\nFirst 5 ref IDs: {ref_ids[:5]}")
print(f"Last 5 ref IDs: {ref_ids[-5:]}")

# V2 might have the same 10 populations but with more samples each
# Let's count the ref IDs positions in the FAM
ref_positions = [i for i, iid in enumerate(iids) if iid.startswith('HG') or iid.startswith('NA')]
uzb_positions = [i for i, iid in enumerate(iids) if not (iid.startswith('HG') or iid.startswith('NA'))]
print(f"\nRef positions: {ref_positions[0]}..{ref_positions[-1]}")
print(f"Uzb positions: {uzb_positions[0]}..{uzb_positions[-1]}")

# Check if they're contiguous
if ref_positions[-1] < uzb_positions[0]:
    print("Layout: all ref samples first, then Uzbeks")
elif uzb_positions[-1] < ref_positions[0]:
    print("Layout: all Uzbek first, then ref")
else:
    print("Layout: interleaved")
