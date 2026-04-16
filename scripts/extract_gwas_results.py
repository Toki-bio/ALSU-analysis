#!/usr/bin/env python3
"""Extract and summarize GWAS results from PLINK2 output."""
import sys

WORKDIR = "/staging/ALSU-analysis/spring2026/gwas"

for tier, fname in [("TIER 3 (Genome-wide)", "tier3_gwas.RPL.glm.logistic.hybrid")]:
    path = f"{WORKDIR}/{fname}"
    try:
        lines = open(path).readlines()
    except FileNotFoundError:
        print(f"{tier}: file not found")
        continue

    header = lines[0].strip().split("\t")
    p_idx = header.index("P")

    valid = []
    for l in lines[1:]:
        parts = l.strip().split("\t")
        try:
            p = float(parts[p_idx])
            if 0 < p < 1:
                valid.append((p, l.strip()))
        except (ValueError, IndexError):
            pass

    valid.sort()
    n_total = len(valid)
    n_gw = sum(1 for p, _ in valid if p < 5e-8)
    n_sug = sum(1 for p, _ in valid if p < 1e-5)

    print(f"=== {tier} ===")
    print(f"Valid SNPs tested:              {n_total}")
    print(f"Genome-wide significant (5e-8): {n_gw}")
    print(f"Suggestive (1e-5):              {n_sug}")
    print()
    print("Top 20 hits:")
    print("\t".join(header))
    for p, line in valid[:20]:
        print(line)
    print()

# Lambda GC
try:
    import math
    path = f"{WORKDIR}/tier3_gwas.RPL.glm.logistic.hybrid"
    lines = open(path).readlines()
    header = lines[0].strip().split("\t")
    p_idx = header.index("P")
    chi2_vals = []
    for l in lines[1:]:
        parts = l.strip().split("\t")
        try:
            p = float(parts[p_idx])
            if 0 < p < 1:
                z_idx = header.index("Z_STAT")
                z = float(parts[z_idx])
                chi2_vals.append(z * z)
        except (ValueError, IndexError):
            pass
    if chi2_vals:
        chi2_vals.sort()
        median_chi2 = chi2_vals[len(chi2_vals) // 2]
        lambda_gc = median_chi2 / 0.4549364
        print(f"Lambda_GC = {lambda_gc:.4f}")
        if lambda_gc > 1.10:
            print("  WARNING: Lambda > 1.10")
        elif lambda_gc < 0.95:
            print("  NOTE: Lambda < 0.95 (conservative)")
        else:
            print("  OK: Lambda within acceptable range")
except Exception as e:
    print(f"Lambda calculation failed: {e}")
