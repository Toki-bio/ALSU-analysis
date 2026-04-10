"""Embed KING verification data into investigation_data_v2.json."""
import json

d = json.load(open('data/investigation_data_v2.json'))
k = json.load(open('data/king_data.json'))

d['king_verification'] = {
    'summary': k['summary'],
    'unexpected_overlay': k['unexpected_overlay'],
    # Only keep genuinely new pairs (clean chips, not on catastrophic)
    'new_pairs': [p for p in k['pairs'] if
        tuple(sorted([p['iid1'], p['iid2']])) not in
        {tuple(sorted([pp['iid1'], pp['iid2']])) for pp in d['pairs_098']}
        and p['kinship'] >= 0.40  # only strong duplicates
    ]
}

print(f"KING summary: {k['summary']}")
print(f"New strong pairs (KING>=0.40, not in PLINK): {len(d['king_verification']['new_pairs'])}")
for p in d['king_verification']['new_pairs']:
    print(f"  {p['iid1']:12s} <-> {p['iid2']:12s}  KING={p['kinship']:.4f}  IBS0={p['ibs0']:.6f}")

with open('data/investigation_data_v2.json', 'w') as f:
    json.dump(d, f)
print("Updated investigation_data_v2.json with king_verification")
