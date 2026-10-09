# ALSU: where the GSA genotyping data is, and where it is described

Written 2026-10-09 after a session asked the user where the "new male samples" were, which showed that no single
file said where each batch lives. Everything marked **verified** was read on DRAGEN (copilot@100.104.25.22, see
skill `server-connections`) on 2026-10-09. Everything marked **stated** was taken from a doc and not re-checked.
Update this file in the same turn whenever the user tells you where something is (global rule).

Participant sample IDs are in some of the files named below. Keep this file local (C:\work\alsu), not in the public
GitHub repo.

## 1. One-screen answer

| I need... | Go to |
|---|---|
| Males (Y, mtDNA, X) | Section 3: batches GWAS2026 (wave 1) and GWAS2026-2, as **VCF**, not PLINK |
| All chromosomes for the original 1,247 (all female) | `ConvSK_raw` PLINK (section 2) |
| Autosome-only analysis sets (PCA, ADMIXTURE, PBS, FST) | Section 4 |
| Imputed genotypes | Section 5 |
| Which doc explains X | Section 6 |
| Sex calls, Y haplogroups already computed | Section 8 |
| Doc contradictions, traps that cost time | Sections 9 and 10 |

## 1a. How the genotyping documentation is organised (so nobody has to ask)

One master file, copies in the places a session will look first (status 2026-10-09: the repo copy was swept into GitHub commit
b83cd58 by another session; the DRAGEN copy was blocked by the permission classifier and is NOT placed). If you change the master, refresh the copies:

| Copy | Why there |
|---|---|
| `C:\work\alsu\DATA_LOCATIONS.md` (this file, the master) | Local project folder; `CLAUDE.md` there links to it |
| `/staging/tmp/ALSU_DATA_LOCATIONS.md` on DRAGEN | For anyone who logs in. (`/staging/ALSU-analysis/` itself is root-owned and refused the write on 2026-10-09; if an admin can place it there, that is the better spot.) |
| `C:\work\ALSU-analysis\DATA_LOCATIONS.md` (GitHub Toki-bio/ALSU-analysis, linked from README) | Public repo, no participant IDs in it |
| Pointers: skill `alsu-project`, memory `alsu-genotype-data-locations` | Loaded at session start |

Rules for keeping it findable:
1. A new batch, file or location gets a row in section 2 or 4 **in the same turn** the user mentions it.
2. Every row says the format (PLINK / VCF / IDAT / GTC), the build, and whether X/Y/MT are present.
3. "verified" means read from the server on a stated date; "stated" means copied from a doc. Never blur them.
4. Results derived from genotypes (sex calls, haplogroups) go in section 8 with their output path, not only in chat.
5. A doc that contradicts this file is listed in section 9 until it is fixed.

## 2. Platform and batches (all on DRAGEN)

Chip: Illumina Global Screening Array **GSA-24v3-0, manifest A2** for all batches. Manifest:
`/staging/ALSU-analysis/Conversion/GSA-24v3-0_A2.bpm` (+ `.csv`). Source: `BpmManifestFilePath` in the DRAGEN
`genotype_call.log` files, quoted in `DRAFT_SECTION_2.1_2.4_2.6_GENOTYPING.md`. **stated**
Genome build of the called data: GRCh38 (rs12913832 at 28,120,472 in `ConvSK_raw.bim`; **verified**).

| Batch | N | Caller | Raw variants | Files (DRAGEN) | Sex in data |
|---|---|---|---|---|---|
| **Original ALSU cohort** ("ConvSK") | 1,247 | Illumina GenomeStudio (not DRAGEN Array; the old report text is wrong on this) | 654,027 | PLINK `/staging/ALSU-analysis/winter2025/PLINK_301125_0312/ConvSK_raw.{bed,bim,fam}`; GenomeStudio export `/staging/ALSU-analysis/winter2025/1_GenomeStudio_qc/FullDataTable.txt` (26 GB); sample sheet `/staging/ALSU-analysis/Conversion/1248_merged_sample_sheet_12022025.csv.csv` | **verified**: `plink --check-sex` gives 1,243 female-like, 0 male-like, 4 ambiguous. Probes: X 27,171, Y 4,138, XY 878, MT 1,138. `.fam` sex column is all 0. Output: `/staging/tmp/uniparental_probe/sexchk.sexcheck` |
| **GWAS2026 wave 1** | 95 of 96 | DRAGEN Array, May 2026 | 613,586 | IDATs `/staging/GWAS2026/GWAS_IDAT/{208993030058,208993030083,209422280033,209422280035}`; sample sheet `/staging/GWAS2026/96_merged_sample_sheet_08052026.xlsx`; working dir `/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/gwas2026/` with `gtc/`, `vcf/`, `gwas2026_merged.vcf.gz` (+`.tbi`), `plink/gwas2026_raw.*`, `logs/`, `qc/`, `GWAS2026_dragen_sample_sheet_95_existing_idats.csv` | **verified** 2026-10-09: VCF has Y 3,822, MT 987, X 27,516 records. Per-sample test (Y call rate >=0.80 with chrX heterozygosity <0.05): **37 male-like**, 58 female-like, no in-between. Case/control labels exist (`Sample_Group`), the study is not identified |
| **GWAS2026-2 (wave 2)** | 95 of 96 | DRAGEN Array, July 2026 | 613,586 | `/staging/GWAS2026-2/` : `gwas2026_2_merged.vcf.gz` (+`.tbi`, 420 MB), `gtc/`, IDAT chip dirs `209422280026, 209422280029, 209422280044, 209422280066`, `GWAS_192_merged_sample_sheet_31072026.xlsx`, `GWAS2026-2_dragen_sample_sheet_96.csv`, `..._95_existing_gtc.csv`, `plink/`, `qc/`, `logs/` | **verified** 2026-10-09: same contigs as above; **26 male-like**, 68 female-like (one sample between: Y call rate 0.1) |
| **48redone** (rescans of 48 original-cohort members whose chips failed) | 48 | DRAGEN Array, Aug 2026 | 613,586 | IDATs `/staging/GWAS2026/48redone/ALSU rescan/{209422280005,209422280032}`; sheet `/staging/GWAS2026/48redone/48_rescan_sample_sheet_18052026.xlsx`; working dir `.../dragen_array_test_20260513/gwas48redone/` (`gtc/`, `gwas48redone_merged.vcf.gz`, `plink/gwas48redone_raw.*`) | **verified**: 0 male-like, 48 female-like |

So the **known male samples are 37 (GWAS2026) + 26 (GWAS2026-2) = 63**, from per-sample chrY call rate plus chrX
heterozygosity in the VCFs. Both signals agree (bimodal). This is my inference from genotype data; the sample sheets
have not been read for a sex column. Per-sample tables: `/staging/tmp/uniparental_probe/{gwas2026_merged,gwas2026_2_merged,gwas48redone_merged}_sexsummary.tsv`
(columns: sample, Y call rate, X heterozygosity). The user said on 2026-10-09 that "males are in new samples"; this
is consistent with it.

Other GTC pools: `.../dragen_array_test_20260513/gtc/` has 1,251 files (original cohort re-called as a test; see
`reports/caller_concordance_48_report_2026-08-27.md` there for GenomeStudio vs DRAGEN Array concordance on 48 samples,
**stated**, not read by me).

## 3. What each file type can and cannot give

- **VCFs of the new batches** (`gwas2026_merged.vcf.gz`, `/staging/GWAS2026-2/gwas2026_2_merged.vcf.gz`): keep X, Y and MT.
  Use these for Y haplogroups and mtDNA. **verified**
- **PLINK sets of the new batches and everything derived from them** (`gwas*/plink/*`, `FULL_MERGED`, `dragen3_*`,
  `cross_gwas2026_alsu_winter/merged.*`): probed 2026-10-09 for chromosome codes above 22; none found. Autosomes only.
  **verified** for `gwas2026/plink/*`, `gwas48redone_raw`, `FULL_MERGED`, `dragen3_step1`, `dragen3_merged`,
  `cross_gwas2026_alsu_winter/merged.bim`, and not checked for the other files.
- **All fam files have sex = 0** (checked for ConvSK_raw, FULL_MERGED, dragen3_*, the gwas* raw sets).
- Coordinates: GRCh38 everywhere. Yhaplo wants GRCh37 (liftover needed). SNAPPY build not verified.

## 4. Derived analysis sets (autosomes)

| Set | N | Path |
|---|---|---|
| Raw merge of all batches | 1,485 | `/staging/ALSU-analysis/spring2026/full_expanded_cohort/FULL_MERGED.*` |
| After QC | 1,268 | `.../full_expanded_cohort/FULL_QC_FINAL.*` |
| New batches only | 238 (143 + 95) | `.../full_expanded_cohort/dragen3_merged.*`, `dragen3_step1.*` (143 = GWAS2026-2 + 48redone) |
| Original-cohort QC chain | 1,247 to 1,056/1,047 | `/staging/ALSU-analysis/winter2025/PLINK_301125_0312/ConvSK_mind20_dedup_snpqc.*`, `/staging/ALSU-analysis/Conversion/OUTPUT/ConvSK/PLINK_*` |
| Cross merge for relatedness | 1,193 (1,098 + 95) | `.../dragen_array_test_20260513/cross_gwas2026_alsu_winter/merged.*` |
| NYGC reference merges, PBS, ADMIXTURE | 1,047 / 1,256 | `.../full_expanded_cohort/{nygc30x,nygc30x_1256,nygc30x_1256_filt,nygc30x_adm}` |

All sample-count definitions: `steps/cohorts_and_sample_sets.html` in `C:\work\ALSU-analysis`.

## 5. Imputation (Michigan Imputation Server, 1000G Phase 3 deep GRCh38, Minimac 4.1.6)

Run 1 winter2025 (1,098): `.../winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz/`.
Run 2 April 2026 (1,093): `/staging/ALSU-analysis/spring2026/imputation/`.
Run 3 August 2026 (1,268): `/staging/ALSU-analysis/spring2026/full_expanded_cohort/imputation_results/`.
Imputation is autosomal; it gives nothing for Y or MT. Details: `IMPUTATION_HISTORY_COMPARISON.md`.

## 6. Where the descriptions live (read these, not me)

Local, `C:\work\alsu\`: `CLAUDE.md` (rules, current state), `DRAFT_SECTION_2.1_2.4_2.6_GENOTYPING.md` (chip, manifest,
which caller for which batch), `IMPUTATION_HISTORY_COMPARISON.md`, `CURRENT_DATA_SUMMARY_FOR_REPORT_UPDATE.md`,
`presentation/DATA_AND_METHODS.md`, `report_output/REPORT_2026_SECTIONS_READY.md`, phenotype sheet
`GWAS от 27.08 - For_Plink.csv` (female pregnancy questionnaire; the 2026 batches mostly lack rows).
Site repo `C:\work\ALSU-analysis`: `COHORT_EXPANSION_2026_08_INTEGRATION_NOTES.md` (batch integration, QC funnel,
bugs), `data-sources.html` (GWAS2026 vs GWAS2026-2 table), `steps/cohorts_and_sample_sets.html`, `steps/step0.html`,
`steps/step_pre.html`, `steps/step1.html`, `old_logs/ConvSK_QC_and_VCF_report.md` and
`old_logs/alsu_genotype_imputation_pipeline_full_technical_log.md` (original conversion). The files in this
paragraph were read only as far as the passages quoted above; the step pages and old logs I did not read in full.
On DRAGEN: `.../dragen_array_test_20260513/reports/caller_concordance_48_report_2026-08-27.md` (**not read**).

## 8. Results derived from the genotypes (uniparental markers)

- **Sex calls** for the new batches: `/staging/tmp/uniparental_probe/*_sexsummary.tsv` (section 2). Original cohort: `sexchk.sexcheck` there.
- **Y haplogroups, 63 males**, 2026-10-09, Yhaplo (commit f72f88f, ISOGG snapshot **2016-01-04**, so deep branches use old names):
  work dir `/staging/tmp/yhaplo_males/` (`venv/`, `yhaplo_src/`, `build_genos.py`, `run/`); calls `run/out/haplogroups.males.txt`;
  local copy `C:\work\alsu\uniparental\y_haplogroups_63males_2026-10-09.txt` (has sample IDs, keep local); method and numbers in
  `C:\work\alsu\uniparental\Y_HAPLOGROUPS.md`. Input was built from the GRCh38 VCF genotypes by liftOver to GRCh37
  (`run/hg38ToHg19.over.chain.gz`).
- **SNAPPY 2.2 on the same 63 males:** `/staging/tmp/snappy_males/real/males_snappy.out` (tool and venv in `/staging/tmp/snappy_males/`);
  58 of 63 agree with Yhaplo at 3-character level, see `Y_HAPLOGROUPS.md`.
- **mtDNA haplogroups, 2026-10-09**, HaploGrep3 3.3.2 (`/staging/tmp/mt_haplo/hg3/`, tree phylotree-rcrs@17.3, `--chip`): 238 samples of the three 2026
  batches `/staging/tmp/mt_haplo/new238_haplo.txt` (MT VCF `mt_new238.vcf.gz`), original 1,247 `/staging/tmp/mt_haplo/convsk_haplo.txt` (ConvSK MT genotypes
  converted by `mt_convsk.py`, 883 sites after strand alignment to the DRAGEN VCF). Coarse (495 of 1,247 calls shorter than 3 characters). Not yet related to ancestry.
- **Public summary page:** `steps/uniparental_markers.html` in the ALSU-analysis repo (no sample IDs); linked from index and the cohorts page.
- HLA: not done (needs imputation; the array has about 7,760 probes in the MHC window).

## 9. Known contradictions in other docs (fix or keep in mind)

- Old report text says the original cohort was converted with DRAGEN Array v4.3. The server holds a GenomeStudio export; DRAGEN Array was used only for the 2026 batches (`DRAFT_SECTION_2.1_2.4_2.6_GENOTYPING.md`).
- `data-sources.html` calls GWAS2026 "96 samples"; the analysed set is 95 (sample sheet `..._95_existing_idats.csv`), and GWAS2026-2 is also 95 of 96.
- The public pages and `cohorts_and_sample_sets.html` give sample counts for autosomal analysis sets only; none of them mentions that Y/MT exist in the VCFs.
- `.fam` sex columns are all 0, so any tool that reads sex from `.fam` sees no males.

## 10. Traps

- `/staging` top level is execute-only. `ls /staging/*GWAS*` fails ("cannot access") even though `/staging/GWAS2026-2`
  exists. A failed glob there is not evidence of absence; `ls` the exact path.
- Short names mislead: "ConvSK" is the original cohort; "dragen3_*" are the 2026 batches; "GWAS2026" is wave 1 and
  "GWAS2026-2" is wave 2; "48redone" is not new people.
- Sample sex must come from genotypes (Y call rate, X heterozygosity, `--check-sex`) or a sample sheet, never from `.fam`.
- Do not run Y/MT tools on the PLINK sets; use the VCFs.
- Read-only rule for these raw directories; write scratch to `/staging/tmp/...` and begin scripts with
  `cd /staging/... || exit 1` (the /home incident of 2026-10-01).
