# Site rebuild on the NYGC 30x reference, 2026-10-04

Scripts that produced the numbers, tables and figures of the rebuilt pages (steps 8-15, relatedness page, manuscripts). Page builders (`build_*.py`) read `results.json` through `sg.py`; they contain no hand-typed result numbers except where a value is quoted from a log and says so on the page.

## Data analysis (run on a workstation from files copied from DRAGEN)
| Script | What it does | Feeds |
|---|---|---|
| `analysis.py` | ten FST pairs (+ leave-one-chromosome-out jackknife), UZB-EUR per-SNP distribution, top-30 SNPs with genes (Ensembl), PBS statistics, top-15 PBS, non-negative mixture fit | steps 9, 10, 14 |
| `figs.py` | mixture-fit test, FST histogram/Manhattan, PBS Manhattan, FST heatmap + MDS | steps 9, 10, 14 |
| `pca_analysis.py` | global PCA groups, nearest-reference labels, PC4 offset, AFR-cluster sample | step 8 |
| `eth_analysis.py` | ancestry by nationality/birthplace, ID-keyed, with permutation and positive controls | step 11 |
| `roh_analysis.py` | FROH versus ancestry components | step 15 |
| `rel_q.py` | ancestry distances for the 1,193 relatedness-set samples on projected K=4, permutation test, sample reuse | relatedness page |
| `ld_fig.py` | LD-decay figure | step 13 |
| `annotate.py`, `annotate_gwas.py` | Ensembl / GWAS Catalog lookups (GWAS helper has a positive-control gate and records failures separately from "none") | step 12 |
| `old8c.py` | NYGC allele frequencies of the eight withdrawn PBS candidates (run on DRAGEN) | PBS candidates page |
| `linkgraph.py` | link audit: pages reachable from index.html, broken internal links | archive index |

## Server jobs (DRAGEN; all write under /staging, never /home)
`ld13.sh` (LD decay UZB vs NYGC groups), `proj_k4.sh` (project samples onto fixed NYGC K=4 frequencies), `fst22.sh` + `fst22b.py` (unascertained-SNP FST check on chr22; Hudson estimator validated against PLINK), `hub.sh` (per-sample missingness/heterozygosity in the relatedness set), `v2_move.sh` (move of the 12 GB V2 data out of /home), `dl_kaz.sh`/`dl2.sh`/`kaz_pipeline.sh` (Kazakh GSA reference), `ibdne.sh` (IBD-based recent Ne).

## Oracle checks used (rule 11 of the project CLAUDE.md)
- ID-keyed phenotype join reproduces 1,057 matches; label permutation gives median p about 0.5; reference-label positive control p < 1e-300.
- ADMIXTURE K=4 Uzbek means (0.50/0.30/0.19) reproduced by an allele-frequency mixture fit (0.49/0.28/0.23, r = 0.994).
- FROH recomputed from `.hom.indiv` matches the page (median 0.0147, 35 above 0.0625).
- Hudson FST from allele frequencies matches PLINK on reference pairs (max difference 0.0044) before it is used for Uzbek pairs.
- Projection onto fixed frequencies reproduces the full-run ancestry (r = 0.94-0.99) before it is used for samples absent from the run.
- GWAS Catalog helper must return associations for a known SNP (rs1800414) before any "none listed" is trusted; first pass had silently hit HTTP 429.
- Every copy to a new location (V2 move) verified by file bytes, file/dir/symlink counts and a checksum diff before the source was removed.
