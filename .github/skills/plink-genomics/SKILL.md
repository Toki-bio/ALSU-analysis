---
name: plink-genomics
description: >
  Population genetics and genomics analysis using PLINK 1.9 and PLINK2. Use this skill
  whenever the user works with genotype data in PLINK binary format (.bed/.bim/.fam),
  VCF files, or needs to run: quality control (missingness, IBD, HWE, MAF), PCA,
  FST analysis, ADMIXTURE preparation, LD pruning, liftover between genome builds,
  merging with reference panels (1000 Genomes), or any GWAS-related preprocessing.
  Also trigger for VCF↔PLINK conversion, population assignment files, or
  interpreting PLINK output logs. This skill incorporates principles derived from
  the ALSU population genomics pipeline developed for the Uzbek cohort, but is
  designed to generalize to any dataset — parameters are always chosen dynamically
  based on data type, sample size, and analysis goal.
---

# PLINK Genomics Skill

## CRITICAL PRINCIPLE: No Hard-Coded Parameters

**Never copy parameters directly from this skill into commands without checking context.**
Every numeric threshold has a decision rule. Always:
1. Read the context (dataset size, data type, analysis goal)
2. Apply the decision rule from the relevant section below
3. If context is missing — ask before running

This applies even to "standard" values. The ALSU pipeline used specific thresholds for
specific reasons; those reasons may not apply to the current task.

---

## Quick Reference: Decision Rules

Before running any command, determine:

| Question | Why it matters |
|---|---|
| How many samples? | Affects IBD computation time; < 200 = fast; large datasets (>5000 samples or >1M SNPs) may benefit from PLINK2 speed, but PLINK 1.9 handles them correctly |
| Genotyped or imputed? | Different QC thresholds (see QC section) |
| hg38 or hg19? | Liftover needed before merging with 1000G Phase 3 |
| Goal: PCA / FST / ADMIXTURE / GWAS? | Different LD pruning parameters |
| PLINK 1.9 or PLINK2 available? | Some commands differ (esp. VCF import, PCA) |

---

## Module 1: Sample QC

Read `references/qc.md` for full parameter decision tables.

### 1.1 Missingness filter

```bash
# Step 1: calculate per-sample missingness
plink --bfile {PREFIX} --missing --out {PREFIX}_miss

# Step 2: identify samples to remove
awk 'NR>1 && $6+0 > {MIND_THRESHOLD} {print $1"\t"$2}' \
  {PREFIX}_miss.imiss > remove_miss.txt

# Step 3: apply filter
plink --bfile {PREFIX} --remove remove_miss.txt --make-bed --out {PREFIX}_mind
```

**MIND_THRESHOLD decision rule:**
- Raw genotyping array (pre-imputation): `0.20` — tolerant, removes only truly bad samples
- Post-imputation QC: `0.05` — stricter, imputed data should be near-complete
- GWAS-ready dataset: `0.02` — very strict
- Unknown/first look: start with `0.20`, inspect histogram, adjust

### 1.2 IBD / duplicate removal

```bash
# LD-prune first (required for IBD)
plink --bfile {PREFIX} \
  --indep-pairwise {WINDOW} {STEP} {R2} \
  --out {PREFIX}_ibd_prune

# Calculate IBD
plink --bfile {PREFIX} \
  --extract {PREFIX}_ibd_prune.prune.in \
  --genome --out {PREFIX}_ibd

# Extract duplicates/relatives
awk 'NR>1 && $10 >= {PIHAT} {print $3, $4}' \
  {PREFIX}_ibd.genome > remove_dups.txt

plink --bfile {PREFIX} --remove remove_dups.txt --make-bed --out {PREFIX}_dedup
```

**Parameter decision rules:**

| Parameter | Duplicates only | First-degree relatives | Standard GWAS | All relateds |
|---|---|---|---|---|
| PI_HAT | ≥ 0.98 | ≥ 0.45 | **≥ 0.185** | ≥ 0.125 |
| Use case | Remove exact duplicates (ALSU default) | Remove sibs/parent-child | **Most common threshold** | Remove all cousins+ |

**Standard GWAS threshold:** PI_HAT ≥ 0.185 removes first- and second-degree relatives
(corresponds to ~halfway between 2nd-degree (0.25) and 3rd-degree (0.125) relatives).
This is the most widely used threshold in population genomics studies.

**LD pruning for IBD** (not PCA — different parameters):
- Window: `50`, Step: `5`, r²: `0.1` — standard for IBD
- ALSU used: `50 5 0.1` → retained ~44,782 SNPs from 654K array

---

## Module 2: SNP QC

Read `references/qc.md` for full SNP filter decision tables.

```bash
plink --bfile {PREFIX} \
  --maf {MAF} \
  --hwe {HWE} \
  --geno {GENO} \
  --chr 1-22 \
  --make-bed \
  --out {PREFIX}_snpqc
```

**Parameter decision rules:**

| Parameter | Conservative | Standard | Lenient | When |
|---|---|---|---|---|
| `--maf` | 0.05 | 0.01 | 0.001 | 0.01 for GWAS prep; 0.001 for rare variant studies |
| `--hwe` | 1e-10 | 1e-6 | 1e-4 | 1e-6 for controls; relax for cases; always apply to controls only if case-control |
| `--geno` | 0.01 | 0.02 | 0.10 | 0.02 standard (Anderson 2010); 0.10 used in ALSU pipeline (forensically recovered) |

**ALSU pipeline used:** `--maf 0.01 --hwe 1e-6 --geno 0.10` for pre-imputation VCF export

---

## Module 3: VCF Conversion

### PLINK binary → VCF (for imputation servers)

```bash
# PLINK 1.9
plink --bfile {PREFIX} --recode vcf --out {PREFIX}

# PLINK2 (preferred for large datasets)
plink2 --bfile {PREFIX} --export vcf --out {PREFIX}

# Split by chromosome (required by Michigan/TOPMed imputation servers)
for chr in {1..22}; do
  bcftools view {PREFIX}.vcf --regions $chr \
    -o chr${chr}.dose.vcf.gz -O z
  bcftools index chr${chr}.dose.vcf.gz
done
```

### VCF (imputed) → PLINK binary

```bash
# PLINK2 — handles dosage field (DS) from imputation servers
plink2 --vcf {VCF_FILE} 'dosage=DS' \
  --double-id \
  --vcf-half-call missing \
  --fa {REFERENCE_FASTA} --ref-from-fa force \
  --make-bed \
  --out {PREFIX} \
  --threads {THREADS} \
  --memory {MEMORY_MB}

# Validate allele order after conversion (always check rs429358 / APOE ε4 as sanity check)
awk '$2=="rs429358" {print "A1:"$5, "A2:"$6}' {PREFIX}.bim
plink2 --bfile {PREFIX} --freq --out {PREFIX}_freq
grep rs429358 {PREFIX}_freq.afreq
# Expected: A1=C (ALT), frequency ~10-15% for Central Asian / European populations
```

**Memory/threads guidelines:**
- Threads: use `nproc` to check available; safe default: 8
- Memory: leave 20% headroom; 200000 MB (200 GB) for large imputed datasets; 8000 for small arrays

---

## Module 4: Post-Imputation QC + Unique IDs

```bash
# Full post-imputation QC
plink2 --bfile {PREFIX} \
  --geno 0.05 --mind 0.05 --maf 0.01 \
  --make-bed --out {PREFIX}_qc \
  --threads 8

# Assign unique CHR:POS:REF:ALT IDs (prevents duplicate ID errors downstream)
plink2 --bfile {PREFIX}_qc \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 50 \
  --make-bed --out {PREFIX}_unique \
  --threads 8
```

**Why unique IDs matter:** Imputed datasets often contain duplicate rsIDs or missing IDs.
CHR:POS:REF:ALT format guarantees uniqueness and enables cross-dataset matching by position.

---

## Module 5: LD Pruning

**Different goals require different pruning parameters:**

```bash
plink2 --bfile {PREFIX} \
  --indep-pairwise {WINDOW} {STEP} {R2} \
  --out {PREFIX}_pruned \
  --threads 8

plink2 --bfile {PREFIX} \
  --extract {PREFIX}_pruned.prune.in \
  --make-bed --out {PREFIX}_for_{GOAL}
```

**Parameter decision table:**

| Goal | Window | Step | r² | Notes |
|---|---|---|---|---|
| PCA (same population, imputed) | 1000kb | 1 | 0.05 | ALSU Step 7: ~88K from 5.4M SNPs |
| PCA (same population, array) | 200kb | 50 | 0.2 | Common alternative for array data |
| PCA (global, multi-population) | 200kb | 1 | 0.05 | ALSU Step 8: ~100-150K |
| ADMIXTURE (imputed) | 50 | 5 | 0.3 | ALSU Step 11: ~380K from 5.4M |
| ADMIXTURE (array) | 50 | 10 | 0.1 | Stricter for lower SNP density |
| IBD calculation | 50 | 5 | 0.1 | ALSU Step 2: ~44K from 654K |
| GWAS covariates | 50 | 5 | 0.2 | Standard for array-based GWAS |

**Note on reviewer suggestion:** Window `200 50 0.2` and `50 10 0.1` are valid alternatives
widely used in the literature. ALSU used `1000kb 1 0.05` for PCA on imputed data
(very aggressive pruning justified by high SNP density). For array data (~650K SNPs),
`200 50 0.2` is more appropriate and faster.

**Note:** Aggressive pruning (r²=0.05) is correct for PCA — removes LD structure so PCs
reflect ancestry, not local haplotypes. ADMIXTURE needs more SNPs (r²=0.3) for stable K estimation.

---

## Module 6: PCA

Read `references/pca.md` for visualization code (R/Python).

### Local PCA (single population)

```bash
plink2 --bfile {PREFIX}_unique \
  --extract {PREFIX}_pruned.prune.in \
  --pca 10 \
  --out {PREFIX}_pca \
  --threads 8
```

**Output:** `{PREFIX}_pca.eigenvec` (sample coordinates), `{PREFIX}_pca.eigenval` (variance explained)

**Interpret eigenvalues:** PC1 variance = eigenval[1] / sum(eigenval) × 100
- ALSU Uzbek-only: PC1=28.93%, PC2=11.51% → strong structure
- If PC1 < 5%: low internal structure (expected for homogeneous cohort)

**GWAS covariates:** Include PC1–PC10 in association models to control for stratification.
For Uzbek data: PC1+PC2 capture ~40% variance; use PC1–PC5 as minimum covariate set.

### Global PCA (with 1000G reference)

See `references/merge_1000g.md` for full workflow including coordinate liftover.

---

## Module 7: FST Analysis

Read `references/fst.md` for genome-wide FST pipeline.

```bash
# Requires merged dataset + population assignment file
# populations.txt format: FID  IID  GROUP (tab-separated)

plink --bfile {MERGED_PREFIX} \
  --fst \
  --within populations.txt \
  --out {PREFIX}_fst
```

**Interpret results:**
- Mean FST Uzbek vs EUR ≈ 0.016–0.020 (ALSU result: weighted FST = 0.0204)
- FST > 0.30: flag as ancestry-informative marker; caution in GWAS
- FST > 0.50: potential selection signal or artifact — annotate with nearest gene

**FST reference values (for context):**

| Pair | Typical FST |
|---|---|
| Within European | 0.002–0.006 |
| Uzbek vs EUR | ~0.020 |
| EUR vs SAS | ~0.020–0.040 |
| EUR vs EAS | ~0.10–0.11 |
| EUR vs AFR | ~0.12–0.15 |

---

## Module 8: Merging with 1000 Genomes

**This is complex — read `references/merge_1000g.md` before running.**

Key steps:
1. Liftover Uzbek data hg38→hg19 (if needed) using UCSC liftOver
2. Per-chromosome: extract EUR/EAS/SAS samples from 1000G VCFs
3. Rename SNPs to CHR:POS format for cross-dataset matching
4. Find overlapping positions via `comm -12`
5. Extract + merge per chromosome (handle 3-allele conflicts with `--exclude missnp`)
6. Concatenate all chromosomes

---

## Module 9: ADMIXTURE Preparation

```bash
# After LD pruning (r²=0.3, see Module 5):
# ADMIXTURE requires PLINK1 bed format and no chromosome code issues

# Run ADMIXTURE for K=2 to K=8 with cross-validation
for K in {2..8}; do
  admixture --cv -j{THREADS} {PREFIX}_for_admixture.bed $K \
    | tee log_K${K}.out
done

# Extract CV errors
grep "CV error" log_K*.out | sort -t: -k2 -n
```

**K selection rules (from ALSU Step 11):**
- Use CV error minimum as primary criterion (not Evanno ΔK — inappropriate for ADMIXTURE)
- Validate with sNMF (LEA package in R) — run 10 replicates, compare cross-entropy
- ALSU global result (Uzbek + 1000G): optimal K=7
- ALSU Uzbek-only result: optimal K=2 (EUR↔EAS admixture cline)
- If CV curve is flat (< 0.1% improvement per K): choose most parsimonious K

---


## Module 4b (run after Module 4): Imputed Data — INFO Score Filter

**This step is required before any analysis on imputed data.**
Imputed variants with low INFO/R² score are unreliable and must be removed first.

```bash
# Michigan/TOPMed imputation servers provide INFO scores in the VCF INFO field
# Filter by INFO score BEFORE converting to PLINK format

# Method 1: bcftools filter (recommended — operates on VCF directly)
bcftools filter -i 'INFO/R2 > {INFO_THRESHOLD}' {VCF_FILE} -Oz -o {VCF_FILE}_filtered.vcf.gz
bcftools index {VCF_FILE}_filtered.vcf.gz

# Method 2: using the .info.gz file from Michigan server
# File format: SNP  REF(0)  ALT(1)  ALT_Frq  MAF  AvgCall  Rsq  ...
awk 'NR>1 && $7 > {INFO_THRESHOLD} {print $1}' chr{N}.info.gz > snps_pass_info.txt
```

**INFO threshold decision rules:**

| Threshold | Use case |
|---|---|
| R² > 0.80 | **Standard — use for most analyses** (GWAS, PCA, FST) |
| R² > 0.30 | Rare variant studies (MAF < 1%) — imputation less reliable for rare variants |
| R² > 0.90 | High-confidence only (fine-mapping, clinical annotation) |

**Per-chromosome loop (Michigan server output):**
```bash
for chr in {1..22}; do
  # Extract passing SNPs from info file
  zcat chr${chr}.info.gz | awk 'NR>1 && $7 > 0.8 {print $1}'     > chr${chr}_pass_info.txt
  
  # Filter VCF
  bcftools view chr${chr}.dose.vcf.gz     --include ID=@chr${chr}_pass_info.txt     -Oz -o chr${chr}_filtered.vcf.gz
  bcftools index chr${chr}_filtered.vcf.gz
  
  echo "Chr${chr}: $(wc -l < chr${chr}_pass_info.txt) variants pass INFO>0.8"
done
```

**After INFO filtering, proceed to PLINK conversion (Module 3) and QC (Module 4).**

---

## Module 4c (run after 4b): Duplicate Variant Removal

Imputed and merged datasets often contain duplicate SNP positions. Must be resolved before
merge with reference panels and before ADMIXTURE/PCA.

```bash
# PLINK2 — identify and remove duplicates
plink2 --bfile {PREFIX}   --rm-dup exclude-mismatch   --make-bed --out {PREFIX}_nodup   --threads 8

# PLINK 1.9 — list duplicates first, then remove
plink --bfile {PREFIX}   --list-duplicate-vars suppress-first   --out {PREFIX}_dupcheck

# Review: {PREFIX}_dupcheck.dupvar
# Remove listed duplicates:
awk '{print $4}' {PREFIX}_dupcheck.dupvar > dup_snps_to_remove.txt
plink --bfile {PREFIX}   --exclude dup_snps_to_remove.txt   --make-bed --out {PREFIX}_nodup
```

**When to run:**
- After imputation + VCF→PLINK conversion (Module 3)
- Before any merge operation (Module 8)
- If `--bmerge` fails with "multiple markers with same name" error

---

## Module 4d (run before any merge): Strand Alignment for Cross-Dataset Merging

**Critical step before any merge with 1000 Genomes or other reference panels.**
Genotyping arrays from different manufacturers may have opposite strand orientations.

### What can go wrong
- **Strand flip:** Array typed A on + strand, reference typed A on − strand → apparent allele mismatch
- **A/T and G/C SNPs:** Cannot be strand-resolved by allele codes alone → must be removed
- **Allele swap:** REF/ALT assignment differs between datasets → flip required

### Step 1: Remove ambiguous A/T and G/C SNPs

```bash
# These SNPs cannot be strand-resolved reliably → remove before merging
python3 << 'PYEOF'
import pandas as pd

bim = pd.read_csv('{PREFIX}.bim', sep='	', header=None,
                  names=['CHR','SNP','CM','POS','A1','A2'])

ambiguous = bim[
    ((bim['A1']=='A') & (bim['A2']=='T')) |
    ((bim['A1']=='T') & (bim['A2']=='A')) |
    ((bim['A1']=='G') & (bim['A2']=='C')) |
    ((bim['A1']=='C') & (bim['A2']=='G'))
]
ambiguous[['SNP']].to_csv('ambiguous_AT_GC_snps.txt', index=False, header=False)
print(f"Ambiguous A/T or G/C SNPs to remove: {len(ambiguous)}")
PYEOF

plink --bfile {PREFIX}   --exclude ambiguous_AT_GC_snps.txt   --make-bed --out {PREFIX}_noAT_GC
```

### Step 2: First merge attempt — identify conflicts

```bash
plink --bfile {YOUR_DATA}_noAT_GC   --bmerge {REF_DATA}_noAT_GC   --make-bed --out merged_test 2>/dev/null

# If merge.missnp file is created → allele conflicts exist
```

### Step 3: Flip strands for conflicting SNPs

```bash
# Attempt strand flip for conflicting SNPs
plink --bfile {YOUR_DATA}_noAT_GC   --flip merged_test-merge.missnp   --make-bed --out {YOUR_DATA}_flipped

# Retry merge
plink --bfile {YOUR_DATA}_flipped   --bmerge {REF_DATA}_noAT_GC   --make-bed --out merged_test2 2>/dev/null
```

### Step 4: Exclude remaining mismatches

```bash
# SNPs that still conflict after flipping → true allele mismatches → exclude
if [ -f merged_test2-merge.missnp ]; then
  plink --bfile {YOUR_DATA}_flipped     --exclude merged_test2-merge.missnp     --make-bed --out {YOUR_DATA}_strand_clean
  
  plink --bfile {REF_DATA}_noAT_GC     --exclude merged_test2-merge.missnp     --make-bed --out {REF_DATA}_strand_clean
  
  # Final merge
  plink --bfile {YOUR_DATA}_strand_clean     --bmerge {REF_DATA}_strand_clean     --make-bed --out {MERGED_PREFIX}
fi
```

**Summary of strand alignment workflow:**
```
Remove A/T + G/C SNPs → First merge attempt → Flip mismatches → 
Retry merge → Exclude remaining conflicts → Final merge
```

**ALSU approach:** Used CHR:POS renaming strategy (Module 8) which avoids strand issues
by matching on genomic position rather than allele codes. Both approaches are valid;
position-based matching is more robust for cross-build datasets.

---

## Module 10: Common Troubleshooting

| Error | Cause | Fix |
|---|---|---|
| "3+ alleles" on merge | Multiallelic or strand-flip conflict | Run Module 4d strand alignment workflow |
| A/T or G/C merge failures | Ambiguous SNPs can't be strand-resolved | Remove them (Module 4d Step 1) before merging |
| "Multiple markers with same name" | Duplicate variant IDs | Run Module 4c duplicate removal |
| Negative FST values | Sampling noise, closely related populations | Normal; treat as ~0 |
| "Base-pair positions unsorted" after liftover | Expected after updating coordinates | Add `--sort-vars` or re-sort bim |
| Very high PC1 variance (>50%) | LD not pruned, or strong population structure | Check pruning was applied |
| ADMIXTURE fails or slow | Too many SNPs | Ensure LD pruning r²≤0.3 applied |
| Duplicate variant IDs | Imputed data without unique ID step | Run Module 4 unique ID assignment first |

---

## File Naming Convention (ALSU style)

```
{COHORT}_{STEP}_{FILTER}.{ext}

Examples:
ConvSK_raw                    ← raw array data
ConvSK_mind20                 ← after missingness filter (mind=0.20)
ConvSK_mind20_dedup           ← after IBD deduplication
ConvSK_mind20_dedup_snpqc     ← after SNP QC
UZB_imputed_HQ_clean          ← post-imputation full dataset
UZB_imputed_HQ_qc             ← post-imputation QC filtered
UZB_imputed_HQ_unique         ← with CHR:POS:REF:ALT IDs
UZB_final_pca                 ← PCA output prefix
```

---

## Reference Files

- `references/qc.md` — Full QC parameter tables, INFO score thresholds, sex-check, heterozygosity
- `references/pca.md` — PCA visualization in R and Python (ggplot2, matplotlib)
- `references/fst.md` — Genome-wide FST pipeline (per-chromosome + merge + Manhattan plot)
- `references/merge_1000g.md` — Full 1000 Genomes merge workflow with liftover
- `references/admixture.md` — ADMIXTURE + sNMF workflow, K selection, visualization
