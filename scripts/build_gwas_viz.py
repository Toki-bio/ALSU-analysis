"""
Build GWAS visualization data from downloaded PLINK2 results.

Reads:
  data/tier3_plot.tsv.gz   — genome-wide results (CHR, POS, ID, OR, SE, P)
  data/tier1_pbs.tsv       — PBS candidate results
  data/tier2_rpl_known.tsv — known RPL variant results

Writes:
  data/gwas_results.json   — compact JSON for step16.html rendering
"""

import gzip, json, math, os, random

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ── Load tier3 genome-wide results ─────────────────────────────
print("Loading tier3 genome-wide results...")

chr_offsets = {}
chr_lengths_hg38 = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
    6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
    11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
    16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
    21: 46709983, 22: 50818468,
}
offset = 0
for c in range(1, 23):
    chr_offsets[c] = offset
    offset += chr_lengths_hg38[c]

all_p = []       # for QQ plot
manhattan = []   # thinned for Manhattan: [genome_pos, -log10p, chr]
top_hits = []    # top 50 SNPs by p-value

random.seed(42)
THIN_P_THRESHOLD = 0.05  # keep all SNPs with p < 0.05; thin the rest

n_total = 0
n_valid = 0

with gzip.open(os.path.join(ROOT, "data", "tier3_plot.tsv.gz"), "rt") as f:
    header = f.readline()
    for line in f:
        parts = line.rstrip("\n").split("\t")
        n_total += 1
        try:
            chrom = int(parts[0])
            pos = int(parts[1])
            snp_id = parts[2]
            or_val = float(parts[3])
            se_val = float(parts[4])
            p = float(parts[5])
        except (ValueError, IndexError):
            continue
        if p <= 0 or p >= 1:
            continue
        if chrom < 1 or chrom > 22:
            continue
        n_valid += 1
        neglog10p = -math.log10(p)
        all_p.append(p)

        # Top hits
        if len(top_hits) < 50 or p < top_hits[-1][0]:
            ci_lo = round(math.exp(math.log(or_val) - 1.96 * se_val), 3) if or_val > 0 else None
            ci_hi = round(math.exp(math.log(or_val) + 1.96 * se_val), 3) if or_val > 0 else None
            top_hits.append((p, chrom, pos, snp_id, round(or_val, 4), ci_lo, ci_hi, round(neglog10p, 3)))
            top_hits.sort()
            top_hits = top_hits[:50]

        # Manhattan thinning: keep all p < 0.001, 10% of 0.001-0.05, 0.5% of rest
        gpos = chr_offsets[chrom] + pos
        if p < 0.001:
            manhattan.append([gpos, round(neglog10p, 2), chrom])
        elif p < 0.05 and random.random() < 0.10:
            manhattan.append([gpos, round(neglog10p, 2), chrom])
        elif random.random() < 0.005:
            manhattan.append([gpos, round(neglog10p, 2), chrom])

print(f"  Total rows: {n_total}, valid: {n_valid}, Manhattan points: {len(manhattan)}")

# ── QQ plot data ──────────────────────────────────────────────
print("Building QQ plot data...")
all_p.sort()
n_qq = len(all_p)

# Thin QQ: keep 2000 evenly spaced + all p < 1e-3
qq_data = []
keep_indices = set(range(0, n_qq, max(1, n_qq // 2000)))
for i, p in enumerate(all_p):
    if i in keep_indices or p < 1e-3:
        expected = -math.log10((i + 0.5) / n_qq)
        observed = -math.log10(p)
        qq_data.append([round(expected, 3), round(observed, 3)])

# Lambda GC — use chi2 approximation without scipy
median_idx = n_qq // 2
median_p = all_p[median_idx]
# For p → chi2(1df): chi2 ≈ (norm.ppf(1-p/2))^2
# Use Wichura's rational approximation for inverse normal
def inv_norm(p_val):
    """Approximate inverse normal CDF for 0 < p < 0.5 (upper tail)."""
    import math as m
    t = m.sqrt(-2 * m.log(p_val))
    c0, c1, c2 = 2.515517, 0.802853, 0.010328
    d1, d2, d3 = 1.432788, 0.189269, 0.001308
    return t - (c0 + c1*t + c2*t*t) / (1 + d1*t + d2*t*t + d3*t*t*t)

z_median = inv_norm(median_p / 2)
lambda_gc = round((z_median ** 2) / 0.4549364, 4)
print(f"  Lambda_GC = {lambda_gc}")

# ── Load tier1 PBS results ─────────────────────────────────────
print("Loading tier1 PBS results...")
tier1 = []
with open(os.path.join(ROOT, "data", "tier1_pbs.tsv")) as f:
    hdr = f.readline().strip().split("\t")
    for line in f:
        parts = line.strip().split("\t")
        d = dict(zip(hdr, parts))
        or_val = float(d["OR"])
        se = float(d["LOG(OR)_SE"])
        tier1.append({
            "chrom": int(d["#CHROM"]),
            "pos": int(d["POS"]),
            "id": d["ID"],
            "ref": d["REF"],
            "alt": d["ALT"],
            "a1_freq": round(float(d["A1_FREQ"]), 4),
            "or": round(or_val, 4),
            "ci_lo": round(math.exp(math.log(or_val) - 1.96 * se), 3),
            "ci_hi": round(math.exp(math.log(or_val) + 1.96 * se), 3),
            "p": float(d["P"]),
            "n": int(d["OBS_CT"]),
        })

# ── Load tier2 known RPL results ──────────────────────────────
print("Loading tier2 known RPL results...")
tier2 = []
with open(os.path.join(ROOT, "data", "tier2_rpl_known.tsv")) as f:
    hdr = f.readline().strip().split("\t")
    for line in f:
        parts = line.strip().split("\t")
        d = dict(zip(hdr, parts))
        or_val = float(d["OR"])
        se = float(d["LOG(OR)_SE"])
        tier2.append({
            "chrom": int(d["#CHROM"]),
            "pos": int(d["POS"]),
            "id": d["ID"],
            "ref": d["REF"],
            "alt": d["ALT"],
            "a1_freq": round(float(d["A1_FREQ"]), 4),
            "or": round(or_val, 4),
            "ci_lo": round(math.exp(math.log(or_val) - 1.96 * se), 3),
            "ci_hi": round(math.exp(math.log(or_val) + 1.96 * se), 3),
            "p": float(d["P"]),
            "n": int(d["OBS_CT"]),
        })

# ── Build output JSON ─────────────────────────────────────────
results = {
    "meta": {
        "n_cases": 273,
        "n_controls": 509,
        "n_snps_tested": n_valid,
        "n_snps_total": n_total,
        "lambda_gc": lambda_gc,
        "covariates": "AGE, BMI, PC1, PC2, PC3",
        "bfile": "UZB_imputed_HQ_clean",
        "date": "2026-04-13",
    },
    "tier1_pbs": tier1,
    "tier2_known": tier2,
    "top_hits": [
        {"p": h[0], "chr": h[1], "pos": h[2], "id": h[3], "or": h[4],
         "ci_lo": h[5], "ci_hi": h[6], "neglog10p": h[7]}
        for h in top_hits
    ],
    "manhattan": manhattan,
    "qq": qq_data,
    "chr_offsets": {str(c): chr_offsets[c] for c in range(1, 23)},
    "genome_length": offset,
    "gw_sig": sum(1 for p, *_ in top_hits if p < 5e-8),
    "suggestive": sum(1 for h in top_hits if h[0] < 1e-5),
    "n_suggestive_total": sum(1 for p in all_p if p < 1e-5),
}

out_path = os.path.join(ROOT, "data", "gwas_results.json")
with open(out_path, "w") as f:
    json.dump(results, f, separators=(",", ":"))
size_mb = os.path.getsize(out_path) / 1024 / 1024
print(f"\nWrote {out_path} ({size_mb:.1f} MB)")
print(f"  Manhattan points: {len(manhattan)}")
print(f"  QQ points: {len(qq_data)}")
print(f"  Top hits: {len(top_hits)}")
print(f"  GW significant: {results['gw_sig']}")
print(f"  Suggestive: {results['n_suggestive_total']}")
print(f"  Lambda_GC: {lambda_gc}")
