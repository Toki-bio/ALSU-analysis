import json
from collections import Counter

D = json.load(open('data/investigation_data_v2.json', encoding='utf-8-sig'))
seven = ['08-365','08-701','08-25','08-495','08-77','08-825','12-11']
sv = {s['iid']: s for s in D['sample_verdicts']}

print("=== 7 EXTRA SAMPLES cross-ref with Step 0 ===")
for s in seven:
    v = sv.get(s, {})
    status = v.get('status', '?')
    fmiss = v.get('fmiss', '?')
    reasons = '; '.join(v.get('reasons', []))
    barcode = v.get('barcode', '?')
    print(f"  {s:10s}  status={status:16s}  F_MISS={fmiss}  chip={barcode}  reasons={reasons}")

print()

# Parse the 99 IIDs from server data
removed_iids = []
in_block = False
for line in open('data/step1_server_data.txt', encoding='utf-8'):
    if '=== REMOVED_99 ===' in line:
        in_block = True
        continue
    if in_block and line.startswith('==='):
        break
    if in_block and line.strip():
        parts = line.split()
        if len(parts) >= 3:
            removed_iids.append(parts[1])

print(f"=== All {len(removed_iids)} removed by Step 0 status ===")
statuses = []
for iid in removed_iids:
    v = sv.get(iid, {})
    statuses.append(v.get('status', 'not_in_investigation'))
c = Counter(statuses)
for k, v in c.most_common():
    print(f"  {k}: {v}")

print()
# Check which barcodes they come from
print("=== Removed 99 by chip barcode ===")
barcodes = []
for iid in removed_iids:
    v = sv.get(iid, {})
    barcodes.append(v.get('barcode', '?'))
bc = Counter(barcodes)
for k, v in bc.most_common(10):
    print(f"  {k}: {v}")
