#!/usr/bin/env python3
"""Compute V2 global ADMIXTURE superpop-level mean Q-values and generate JS data."""
import json, re
from collections import defaultdict

# Parse the mapping file
with open("fullmap_out.txt", encoding="utf-8") as f:
    raw = f.read()

# Strip SSH banner  
start = raw.find("=== FULL_MAPPING ===")
end = raw.find("=== END ===")
mapping_text = raw[start + len("=== FULL_MAPPING ==="):end].strip()

lines = mapping_text.split('\n')
print(f"Mapping lines: {len(lines)}")

# Parse: IID\tSUPERPOP
sample_pops = []
for line in lines:
    parts = line.strip().split('\t')
    if len(parts) == 2:
        sample_pops.append(parts[1])

print(f"Parsed sample populations: {len(sample_pops)}")
assert len(sample_pops) == 3595, f"Expected 3595, got {len(sample_pops)}"

# Count populations
from collections import Counter
pop_counts = Counter(sample_pops)
print(f"\nPopulation counts:")
pop_order = ['AFR', 'EUR', 'SAS', 'EAS', 'AMR', 'UZB']
for pop in pop_order:
    print(f"  {pop}: {pop_counts.get(pop, 0)}")

# For each K, compute mean Q per superpop
print("\n" + "="*60)
for k in range(2, 9):
    qfile = f"data/global/global_v2_admix.{k}.Q"
    q_rows = []
    with open(qfile) as f:
        for line in f:
            vals = [float(x) for x in line.strip().split()]
            q_rows.append(vals)
    
    print(f"\nK={k}:")
    for pop in pop_order:
        indices = [i for i, p in enumerate(sample_pops) if p == pop]
        if not indices:
            continue
        means = [sum(q_rows[i][j] for i in indices) / len(indices) for j in range(k)]
        means_str = ','.join(f'{m:.4f}' for m in means)
        print(f"  {pop:5s} (n={len(indices):4d}): [{means_str}]")

# Generate JS-ready data
print("\n\n=== JS DATA FOR STEP11 ===")
print(f"const popNames = {json.dumps(pop_order)};")
print(f"const popN = {json.dumps([pop_counts[p] for p in pop_order])};")
print(f"const popSuper = {json.dumps(pop_order)};")

# CV errors
cv_errors = [
    {"k": 2, "cv": 0.31043},
    {"k": 3, "cv": 0.29992},
    {"k": 4, "cv": 0.29703},
    {"k": 5, "cv": 0.29503},
    {"k": 6, "cv": 0.29458},
    {"k": 7, "cv": 0.29442},
    {"k": 8, "cv": 0.29422},
]
print(f"const cvErrors = {json.dumps(cv_errors)};")

# qData
print("const qData = {")
for k in range(2, 9):
    qfile = f"data/global/global_v2_admix.{k}.Q"
    q_rows = []
    with open(qfile) as f:
        for line in f:
            vals = [float(x) for x in line.strip().split()]
            q_rows.append(vals)
    
    pop_means = []
    for pop in pop_order:
        indices = [i for i, p in enumerate(sample_pops) if p == pop]
        means = [round(sum(q_rows[i][j] for i in indices) / len(indices), 4) for j in range(k)]
        pop_means.append(means)
    
    print(f"    {k}: {json.dumps(pop_means)},")
print("};")
