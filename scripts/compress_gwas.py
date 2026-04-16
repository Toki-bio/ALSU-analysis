#!/usr/bin/env python3
"""Extract plotting columns from GWAS results and gzip."""
import gzip

INFILE = "/staging/ALSU-analysis/spring2026/gwas/tier3_gwas.RPL.glm.logistic.hybrid"
OUTFILE = "/staging/ALSU-analysis/spring2026/gwas/tier3_plot.tsv.gz"

with open(INFILE) as fin, gzip.open(OUTFILE, "wt") as fout:
    header = fin.readline().strip().split("\t")
    ci = header.index("#CHROM")
    pi = header.index("P")
    idi = header.index("ID")
    ori = header.index("OR")
    sei = header.index("LOG(OR)_SE")
    fout.write("CHR\tPOS\tID\tOR\tSE\tP\n")
    n = 0
    for line in fin:
        parts = line.strip().split("\t")
        try:
            p = parts[pi]
            if p == "NA":
                continue
            fout.write(f"{parts[ci]}\t{parts[1]}\t{parts[2]}\t{parts[ori]}\t{parts[sei]}\t{p}\n")
            n += 1
        except Exception:
            pass
    print(f"Wrote {n} rows to {OUTFILE}")
