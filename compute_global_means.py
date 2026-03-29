#!/usr/bin/env python3
"""Compute mean Q-values per population from V2 global ADMIXTURE Q files + FAM."""
import os

# Read FAM to get sample -> population mapping
fam_lines = open("data/global/global_v2_admix.fam").read().strip().split('\n')
print(f"FAM: {len(fam_lines)} samples")

# The FAM file format: FID IID ... — FID encodes the population
# Let's check: 
for i, line in enumerate(fam_lines[:5]):
    parts = line.split()
    print(f"  Sample {i}: FID={parts[0]}, IID={parts[1]}")
print("  ...")
for i, line in enumerate(fam_lines[-5:]):
    parts = line.split()
    print(f"  Sample {len(fam_lines)-5+i}: FID={parts[0]}, IID={parts[1]}")

# Count populations
from collections import Counter
pops = [line.split()[0] for line in fam_lines]
pop_counts = Counter(pops)
print(f"\nPopulations ({len(pop_counts)}):")
for pop, n in sorted(pop_counts.items(), key=lambda x: -x[1]):
    print(f"  {pop}: {n}")

# Define population order matching V1 page
pop_order = ['YRI','CEU','FIN','GBR','IBS','TSI','GIH','PJL','CHB','JPT','UZBEK']

# For each K, compute mean Q per population
for k in range(2, 9):
    qfile = f"data/global/global_v2_admix.{k}.Q"
    q_rows = []
    with open(qfile) as f:
        for line in f:
            vals = [float(x) for x in line.strip().split()]
            q_rows.append(vals)
    assert len(q_rows) == len(fam_lines), f"K={k}: Q rows {len(q_rows)} != FAM {len(fam_lines)}"
    
    # Group by population
    pop_qs = {pop: [] for pop in pop_order}
    for i, (pop, qvals) in enumerate(zip(pops, q_rows)):
        # Map FID to our pop names
        if pop in pop_qs:
            pop_qs[pop].append(qvals)
        else:
            # Check if it's the Uzbek cohort with a numeric FID
            # Let's see what unmapped looks like
            pass
    
    # Check if UZBEK mapped
    mapped = sum(len(v) for v in pop_qs.values())
    if mapped < len(fam_lines):
        # Some didn't map — print unmapped FIDs
        unmapped_fids = set(pops) - set(pop_order)
        print(f"\nK={k}: {mapped} mapped, {len(fam_lines)-mapped} unmapped. FIDs: {unmapped_fids}")
        break

    print(f"\nK={k}: Mean Q per population:")
    for pop in pop_order:
        qs = pop_qs[pop]
        if not qs:
            print(f"  {pop}: NO DATA")
            continue
        n = len(qs)
        means = [sum(row[j] for row in qs) / n for j in range(k)]
        mean_str = ','.join(f'{m:.4f}' for m in means)
        print(f"  {pop} (n={n}): [{mean_str}]")
