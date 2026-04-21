#!/usr/bin/env python3
"""Process plink2 logistic.hybrid output: count, median p, top-10."""
import sys, math, heapq, os

def process(path, label):
    if not os.path.exists(path):
        print(f"  [{label}] missing: {path}")
        return
    n = 0; gw = 0; sug = 0
    ps = []
    top = []  # min-heap of (-p, line)  -> retain 10 smallest p
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            pi = header.index("P")
        except ValueError:
            print(f"  [{label}] no P column"); return
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= pi: continue
            try:
                p = float(parts[pi])
            except ValueError:
                continue
            if not (0 < p < 1): continue
            n += 1
            ps.append(p)
            if p < 5e-8: gw += 1
            if p < 1e-5: sug += 1
            if len(top) < 10:
                heapq.heappush(top, (-p, line.rstrip("\n")))
            elif -top[0][0] > p:
                heapq.heapreplace(top, (-p, line.rstrip("\n")))
    if not n:
        print(f"  [{label}] no valid p"); return
    ps.sort()
    median_p = ps[n//2]
    chi_med = -2 * math.log(median_p)  # NOT lambda — see below
    # Proper lambda_GC for 1df chi2: lambda = chi2_inv(1-median_p, 1) / 0.4549
    # Use approximation via inverse error function: chi2_1df_inv(1-p) = 2*erfcinv(p)^2 ... need scipy or approx
    # Beasley-Springer-Moro approximation
    def erfcinv(x):
        # Acklam's inverse normal then convert: erfcinv(x) = -inv_norm(x/2)/sqrt(2)
        from math import sqrt, log
        # Use rational approximation
        # For 0<x<2
        if x <= 0 or x >= 2:
            return float('nan')
        # Abramowitz & Stegun 26.2.23 via inverse normal
        p = x / 2.0
        if p < 0.5:
            t = math.sqrt(-2.0*math.log(p))
            num = 2.515517 + 0.802853*t + 0.010328*t*t
            den = 1 + 1.432788*t + 0.189269*t*t + 0.001308*t*t*t
            z = -(t - num/den)
        else:
            t = math.sqrt(-2.0*math.log(1-p))
            num = 2.515517 + 0.802853*t + 0.010328*t*t
            den = 1 + 1.432788*t + 0.189269*t*t + 0.001308*t*t*t
            z = (t - num/den)
        return -z / math.sqrt(2.0)
    chi2_median = 2.0 * (erfcinv(median_p) ** 2)
    lambda_gc = chi2_median / 0.4549364
    print(f"  [{label}] tested={n}  GW(<5e-8)={gw}  suggestive(<1e-5)={sug}")
    print(f"  [{label}] median_p={median_p:.6g}  lambda_GC={lambda_gc:.4f}")
    # Write top 10
    out = path.replace(".RPL.glm.logistic.hybrid", "_top10.tsv")
    with open(out, "w") as o:
        o.write("\t".join(header) + "\n")
        for negp, ln in sorted(top, key=lambda x: -x[0]):
            o.write(ln + "\n")
    print(f"  [{label}] top10 -> {out}")

if __name__ == "__main__":
    base = "/staging/ALSU-analysis/spring2026/gwas/sensitivity/results"
    for name in ["S1_strict3", "S2_review_as_case", "S3_gray_as_case", "S4_probable"]:
        process(f"{base}/{name}.RPL.glm.logistic.hybrid", name)
