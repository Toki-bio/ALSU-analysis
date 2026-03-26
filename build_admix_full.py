#!/usr/bin/env python3
"""Build ADMIX_DATA JSON with K=2-8 for step10.html."""
import json, os

# Load existing K=2-5 data to get ids and sort_order
with open("data/admix_v2_data.json", encoding="utf-8") as f:
    base = json.load(f)

ids = base["ids"]
sort_order = base["sort_order"]
n = len(ids)
print(f"Samples: {n}, sort_order len: {len(sort_order)}")

result = {"ids": ids, "sort_order": sort_order}

# Copy K=2-5 from existing
for k in range(2, 6):
    key = f"k{k}"
    result[key] = base[key]
    print(f"K={k}: {len(base[key])} samples")

# Parse K=6,7,8 from new Q files
for k in range(6, 9):
    qfile = f"data/UZB_v2_admix.{k}.Q"
    rows = []
    with open(qfile) as f:
        for line in f:
            vals = [round(float(x), 6) for x in line.strip().split()]
            rows.append(vals)
    assert len(rows) == n, f"K={k}: expected {n} rows, got {len(rows)}"
    result[f"k{k}"] = rows
    print(f"K={k}: {len(rows)} samples, {k} components each")

# Save full JSON
with open("data/admix_v2_data_full.json", "w") as f:
    json.dump(result, f, separators=(',', ':'))
print(f"Full JSON: {os.path.getsize('data/admix_v2_data_full.json')} bytes")

# Build JS embedding line
js_line = "const ADMIX_DATA = " + json.dumps(result, separators=(',', ':')) + ";"
with open("data/admix_v2_jsline_full.txt", "w") as f:
    f.write(js_line)
print(f"JS line: {len(js_line)} chars")
