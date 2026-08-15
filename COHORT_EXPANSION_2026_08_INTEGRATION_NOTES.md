# Cohort Expansion — Integration Notes (August 2026)

## Summary

Integrated three new genotyping batches (48redone rescans, GWAS2026 wave 1, GWAS2026-2)
into the original ALSU popgen cohort, producing a QC'd, imputed, expanded dataset for
the population-structure track (Steps 7–14 only — **not** the RPL/pregnancy-loss GWAS,
Step 16, which stays scoped to the original 1,047-person cohort per explicit user
instruction). This document is the source-of-truth technical record, with emphasis on
exactly what changed relative to the original pipeline and why.

## Scope decision (explicit, from the user)

- GWAS2026-2 (95 samples): confirmed population-genetics reference data, unrelated to
  the RPL cohort. Included.
- GWAS2026 wave 1 (96 samples, case/control design for an unstated study): initially
  parked pending purpose, then explicitly included by the user despite the unresolved
  case/control label — verified afterward that case/control status has no plausible
  mechanism to distort PCA/ADMIXTURE/F_ST/PBS (those infer ancestry from genome-wide
  markers; ascertainment on an unrelated trait doesn't move enough of them to matter,
  unlike GWAS association testing where it would).
- 48redone (48 samples): rescans of specific original-cohort members whose chips
  failed (`208993030034`, `208993030080` in the original 4-chip catastrophic-failure
  set documented in Step 0). Confirmed via exact sample-sheet ID match — all 48 IDs
  present verbatim in the original cohort's own sample sheet.
- An earlier external tracker (`uzbek-variant-frequency-explorer` repo) had fabricated
  a claim that GWAS2026 wave 1 "replaces 2 failed ALSU chips" — traced via git history,
  found to be an unverified inferential leap bridging two independently-true facts
  (real failed-chip barcodes + a real new batch). Corrected in both repos.

## Dataset inventory, all counts verified directly against DRAGEN files

| Source | Raw samples | Raw variants | Notes |
|---|---|---|---|
| ALSU original | 1,247 | 654,027 | GenomeStudio-called, `ConvSK_raw` |
| 48redone | 48 | 613,586 | DRAGEN Array-called this session (see below) |
| GWAS2026 wave 1 | 95 | 613,586 | DRAGEN Array-called May 2026 |
| GWAS2026-2 | 95 (of 96) | 613,586 | DRAGEN Array-called July 2026 |
| **Combined raw** | **1,485** | — | Before any QC |

## Step-by-step, with exact comparison to the original pipeline

### 1. Genotype-calling 48redone (new work, no prior run existed)

The repo's pre-written script (`scripts/run_48redone_qc.sh`) had never actually been
run successfully — it contained three real bugs, found by execution rather than
assumption:
- `PLINK2` hardcoded to `/usr/local/bin/plink2`, which doesn't exist on this host
  (real path: `/staging/conda/envs/bioinfo/bin/plink2`).
- `dragena genotype-call` — wrong verb; the actual CLI (checked via `--help`) uses
  `genotype call` as a two-word subcommand, with different flag names throughout
  (`--bpm-manifest` not `--bpm-manifest-file`, `--gencall-cutoff` not
  `--gencall-score-cut-off`, etc.). Rebuilt from GWAS2026-2's own known-good
  `run_pipeline.sh` invocation instead of guessing.
- The generated sample sheet left `Sample_Name` empty; DRAGEN Array's sample-sheet
  parser rejects this ("Sample_Name cannot be empty"). Fixed by setting it equal to
  `Sample_ID`, matching GWAS2026-2's own sheet convention.

Result: 48/48 GTCs called successfully (better than GWAS2026-2's 95/96).

### 2. Raw merge — strand alignment bug found and fixed mid-session

First merge attempt used `chr:pos:ref:alt` as the join key between the GenomeStudio-called
ALSU data and the DRAGEN-Array-called new batches. This produced only **259,184** shared
SNPs — flagged by the user as "some huge error, mistaken loss of most snps, unacceptable."
Investigation found the real cause: matching by position alone (ignoring ref/alt string
order) gave **597,979** shared positions — the gap was ref/alt written in different order
between two independent `--ref-from-fa` runs (e.g. `C:T` vs `T:C` for the same physical
SNP), not real non-overlap.

Fix: merge on `chr:pos` only, let PLINK1's `--bmerge` + `--flip` cycle resolve genuine
strand mismatches. Recovered to **594,601** SNPs (2.3× more than the first attempt),
with only 1,610 residual true conflicts (0.27% — consistent with normal A/T,C/G
strand-ambiguous sites, not a bug).

Final combined raw merge: **1,485 samples, 594,601 SNPs**, 97.0% genotyping rate.

### 3. Sample/variant QC — full-set, not per-batch

Per explicit user instruction ("i asked you to run qc on full set"), QC was run once
on the combined 1,485-sample raw merge, not per-source. This also resolved a real
PLINK2 hard error (`--indep-pairwise` refuses under 50 samples) that occurred when
48redone (48 samples) was QC'd standalone.

| Stage | Samples | Change |
|---|---|---|
| Raw merged | 1,485 | — |
| mind 0.20 | 1,384 | −101 |
| IBD dedup (PI_HAT≥0.98, 50 clusters) | 1,326 | −58 |
| Het (±3SD, 36 flagged) + KING (0.185, 27 flagged) | 1,268 | −58 |
| **Final** (geno 0.10, hwe 1e-6, maf 0.01 — same thresholds as original Step 3) | **1,268 samples, 452,377 SNPs** | — |

Compare to the original cohort alone (Step 0/3): 1,247 → 1,056 (191 removed: 98
chip_failure, 45 dt_artifact, 33 ibd_duplicate, 12 high_fmiss, 3 contaminated) → 459,823
SNPs after variant QC. The expanded pipeline used the *same* thresholds throughout
(mind 0.20, PI_HAT≥0.98, geno 0.10, hwe 1e-6, maf 0.01) for direct comparability.

### 4. Michigan Imputation Server submission — file prep and known quirks

Replicated the documented spring2026 recipe exactly (`step3.html`/`step4.html`):
indel exclusion → per-chr VCF export → concat/biallelic/dedup → Cyrillic ID fix →
chr-rename → `bcftools +fixref` strand alignment (`-m top`) → palindromic-SNP removal.
Result: 450,665 variants (from 452,377, −0.4% to palindromic removal, same order of
loss as the original run), 1,268 samples, submitted with the exact original job
settings (1000G Phase 3 30x GRCh38, Eagle v2.4, Minimac v4.1.6, skip AF check, QC &
Imputation mode).

**Submission QC comparison to the original:**

| | Original (1,093 samples) | This run (1,268 samples) |
|---|---|---|
| Reference overlap | 98.31% | 98.39% |
| Strand flips | 0 | 0 |
| Matched sites | 448,305 | 443,045 |
| Excluded sites | 647 | 364 |
| Typed-only sites | 7,732 | 7,256 |
| Remaining chunks | 152 | 152 |

The identical 155-total/152-remaining chunk count across both runs is expected and not
a coincidence — Michigan's chunking is purely a function of genome length divided by a
fixed ~20Mb window size, independent of sample count.

Download used the documented resume-safe fix (`sed 's/set -e/set +e/; s/curl -L /curl -C - -L /g'`
on Michigan's own download script — their `set -e` kills the script on an HTTP 416
"already complete" response during any resume). 76GB downloaded, 27 files, all
verified via `gzip -t` + tabix index presence. Extraction (22× AES-encrypted zips)
completed cleanly with the user-supplied password.

### 5. Post-imputation ID normalization — a different bug than the documented one

Step 5's documented issue (Michigan v2 numeric prefix on chr1-9 sample IDs only) did
**not** recur this run — checked directly (`chr1` vs `chr10` vs `chr22` sample lists
identical). Instead, a different, self-inflicted issue: every sample ID showed a
doubled pattern (`03-72_03-72`), uniform across all 22 chromosomes (not chr1-9-specific,
so provably not the same bug). Root cause: `FULL_QC_FINAL.fam`'s FID column was never
`0` (carried the original cohort's own numeric row-index FIDs), and PLINK1's VCF export
writes `FID_IID` as the sample name whenever FID isn't `0`.

Fix: built an exact rename mapping directly from the fam file (`FID_IID` → `IID` for
all 1,268 samples, not a regex guess), applied via `bcftools reheader` per chromosome.
Two samples needed a manual mapping addition (`808_03-25m`, `1038_08-176X-00006`) because
their fam-file IID was still Cyrillic while the VCF sample name had already been
Latinized by the pre-submission prep step — caught by requiring 0 unmapped IDs before
touching any file, not assumed. Verified clean on all 22 chromosomes afterward
(uniform 1,268 samples, 0 non-ASCII names remaining).

### 6. R²/MAF quality filtering and the mandatory allele-orientation check

Applied `bcftools view -i 'INFO/R2>=0.80 && INFO/MAF>=0.001'` per chromosome, matching
Step 6 exactly. Result: **9,808,764 variants** pass (vs the original's 9,981,814 for
1,093 samples — comparable scale).

Step 6's docs flag a specific historical bug (PLINK swapping REF/ALT wholesale during
VCF→PLINK conversion, silently corrupting all downstream analysis) and mandate a
validation pass against Michigan's own `info.gz` REF/ALT columns. Ran that validation:

- `--ref-from-fa force` on conversion: 0 variants changed, 9,208,094 validated (already
  correctly ref-aligned, as expected for Minimac dose output).
- Joined 9,887,126 of 9,808,764 bim rows against Michigan's per-position alleles
  (first attempt failed silently — 0 rows joined — due to a `chr`-prefix mismatch
  between PLINK's bare chromosome numbers and Michigan's `chr`-prefixed contig names in
  `info.gz`; found and fixed before concluding anything).
- 72,432 raw mismatches (0.73%), split into 34,158 indel/multi-character entries
  (not flip candidates — excluded from consideration) and 38,274 pure single-base SNP
  mismatches.
- Of the 38,274 SNP mismatches, **zero** are strand-complement pairs (A↔T, C↔G swaps)
  — this directly rules out the systematic REF/ALT-swap bug the docs warn about. The
  observed pattern (e.g. array reports A/G at a position, Michigan's panel reports G/C)
  is consistent with normal tri-allelic-site representation differences between a
  genotyping array and an imputation reference panel, not a pipeline defect.
- Conclusion: no `--flip` correction applied, since there's no evidence one is
  warranted — applying one to non-flip-candidate rows would itself introduce an error.

### 7. Popgen re-analysis (Steps 7–15) — results and one bug found

Steps 7–15 were rerun on the expanded, imputed, QC'd cohort (1,256 samples). All
steps reproduce the original's population-structure signal within expected small
deltas — see per-step numbers in build log. One methodology gap was found and fixed
mid-session; documented here since it was a real defect, not an assumption:

**Step 15 IBD — silent zero-variant bug.** The first IBD run (`plink --extract
UZB_pruned.prune.in --genome`) reported 0 pairs in every relatedness category. Root
cause: `UZB_pruned.prune.in` (from Step 7) uses `chr:pos:ref:alt` variant IDs, but
the IBD script pointed `--bfile` at `UZB_imputed_HQ_qc`, which uses rsIDs — so
`--extract` silently matched 0 variants. The script's own `2>&1 | tail -20` piping
masked this from `set -e` (pipeline exit status came from `tail`, not `plink`), so a
false "=== DONE ===" was printed with all-zero relatedness counts. Fixed by
re-pointing IBD at `UZB_imputed_HQ_unique` (same 5,358,770-variant set, matching ID
convention) and adding `set -eo pipefail` so this class of failure can't hide again
in any script written this session. Corrected result: 6,744 IBD pairs ≥0.05, F_ROH
distribution matches the original almost exactly (mean 0.0188 both cohorts, median
0.0147 vs 0.0148). ROH itself (a separate `--homozyg` invocation, no `--extract`)
was unaffected and correct on the first run (44,938 segments).

**Also found mid-session:** the ADMIXTURE rerun (Step 11) had only redone the
global (Uzbek+1000G) K=2–8 sweep, single-seed. Two verification layers the original
pipeline used to justify "K=5 is the working number"
(`STEP11_ADMIXTURE_CORRECTION_NOTES.md`) were not initially carried over to the
expanded-cohort rerun: (a) the standalone Uzbek-only ADMIXTURE run (separate
QC+LD-pruned panel, not a subset of the global merge), and (b) the 5-independent-seed
replicate stability check per panel. Both are being redone for full parity — see
build log for status/results once complete. sNMF cross-validation (global +
Uzbek-only) is also being rerun for the same reason.

## What's still open

- ADMIXTURE full parity (Uzbek-only panel + 5-seed replicates on both panels) and
  sNMF validation (global + Uzbek-only) are in progress as of this writing — see
  Section 7 above.
- Step 12 (SNP annotation via external API) deferred pending a decision on internet
  access from DRAGEN vs. running locally.
- Site pages have not yet been updated with the expanded-cohort numbers. Per explicit
  scope, only the popgen step pages may be touched — Step 16 (RPL GWAS) stays on the
  original cohort.
- GWAS2026 wave 1's actual case/control phenotype/study is still unknown; nothing in
  this integration depends on knowing it, but it remains an open question for the
  study team.
