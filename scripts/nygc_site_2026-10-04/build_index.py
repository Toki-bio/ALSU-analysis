from sg import *

t = read("index.html")
P = R["pairs"]; S_ = R["pbs"]
reps = [
    ("Input: UZB + 1000 Genomes reference", "Input: UZB + NYGC 30x reference (EUR/EAS/SAS/AFR)"),
    ("Output: Global PCA, ancestry inference | Jan 4", "Output: Global PCA, 3,968 samples, 64,959 SNPs | rebuilt Oct 2026"),
    ("Step 9: Fst Analysis (UZB vs EUR)", "Step 9: Genome-wide Fst (10 pairs)"),
    ("Input: UZB + 1000G EUR (376K SNPs)", "Input: UZB + NYGC reference (82,744 SNPs)"),
    ("Output: Genome-wide Fst = 0.020 | Oct–Nov", f"Output: UZB–EUR {P['EUR_UZB']['weighted_all']:.4f}, UZB–SAS {P['SAS_UZB']['weighted_all']:.4f} | Oct 2026"),
    ("Input: UZB + EUR/EAS/SAS/AFR (376K SNPs)", f"Input: UZB + NYGC EUR/EAS/SAS/AFR ({S_['n']:,} SNPs)"),
    ("Output: PBS (old run DEPRECATED)", f"Output: PBS, no candidates (max {S_['max']:.3f}) | Oct 2026"),
    ("Input: UZB + 1000G (376K SNPs)", "Input: UZB + NYGC reference (64,959 SNPs)"),
    ("Output: K=5 old run DEPRECATED; NYGC rerun current", "Output: K=4: ~50% EUR / 30% EAS / 19% SAS-like | Oct 2026"),
    ("Step 12: PBS SNP Annotation", "Step 12: Annotation of top PBS / Fst SNPs"),
    ("Input: old PBS SNPs (DEPRECATED)", "Input: top PBS and Fst SNPs (NYGC)"),
    ("Output: VEP + GWAS + GTEx results | Mar 2026", "Output: gene, consequence, GWAS Catalog lookups | Oct 2026"),
    ("Output: 401 independent loci | Mar 2026", "Output: LD decay, UZB vs NYGC groups | Oct 2026"),
    ("Output: 5×5 FST matrix + MDS | Mar 2026", "Output: 5×5 FST matrix + MDS (NYGC) | Oct 2026"),
    ("Input: UZB BED (1,047 × 5.41M SNPs)", "Input: UZB BED (1,256 × 5.36M SNPs)"),
    ("Output: 36.7K ROH | F_ROH | IBD | Mar 2026", "Output: 44.9K ROH | F_ROH | IBD | Aug 2026"),
    ("K=5 map (DEPRECATED)", "K=4 map (NYGC)"),
]
for a, b in reps:
    t = wsub(t, a, b)
write("index.html", t)
print("index ok")
