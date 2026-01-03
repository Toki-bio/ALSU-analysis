# ALSU Genotype Imputation Pipeline - Complete Summary

**Project**: ALSU Uzbek Pregnancy Loss Cohort  
**Date**: December 2025  
**Status**: ✅ Complete - Ready for Population Genetics Analyses

---

## Executive Summary

Successfully processed 1,247 Uzbek individuals from raw Illumina array data through comprehensive QC, imputation, and post-imputation filtering to create a high-quality dataset of **1,074 samples** with **10.8 million variants**.

---

## Pipeline Overview

```
Raw Data (Nov 2025)
├─ 1,247 samples
├─ 654,027 variants
└─ Illumina Array (hg38)

    ↓ PHASE 1: PRE-IMPUTATION QC

Step 1: Sample Missingness (−92 samples)
Step 2: IBD Deduplication (−57 samples)
Step 3: SNP QC for Imputation (−180K variants)
Step 4: Michigan Imputation Server (+58M variants)

    ↓ PHASE 2: POST-IMPUTATION QC

Step 5: Heterozygosity Filter (−24 samples)
Step 6: Variant Quality Filter + Merging

    ↓ FINAL DATASET

UZB_imputed_HQ_clean
├─ 1,074 samples (13.9% removed via QC)
├─ 10,846,569 variants (23× expansion)
└─ 100% genotyping rate
```

---

## Final Dataset Specifications

### Primary Dataset

**Name**: `UZB_imputed_HQ_clean`  
**Format**: PLINK binary (.bed/.bim/.fam)  
**Location**: `/staging/ALSU-analysis/winter2025/PLINK_301125_0312/michigan_ready_chr/imputation_results/unz/filtered_clean/`

**Dimensions**:
- Samples: 1,074
- Variants: 10,846,569
- Chromosomes: 1-22 (autosomes)
- Genotyping Rate: 100%

**Quality Criteria Met**:
- ✅ R² ≥ 0.8 (imputation quality)
- ✅ MAF ≥ 0.001 (Uzbek cohort)
- ✅ |F| ≤ 0.2 (heterozygosity)
- ✅ PI_HAT < 0.98 (no duplicates)
- ✅ 100% genotyping (imputed)

---

## Quality Metrics

### Sample Quality

| Metric | Value | Notes |
|--------|-------|-------|
| Final Sample Count | 1,074 | From 1,247 original |
| Samples Removed | 173 (13.9%) | 92 missingness + 57 duplicates + 24 het outliers |
| Mean Heterozygosity F | 0.005 | Near zero (expected for outbred) |
| F Coefficient Range | -0.198 to +0.146 | Tight distribution after filtering |
| Genotyping Rate | 100% | No missingness (imputed data) |

### Variant Quality

| Metric | Value | Notes |
|--------|-------|-------|
| Total Variants | 10,846,569 | Typed + imputed |
| R² ≥ 0.9 | 6,910,161 (69%) | Excellent imputation quality |
| R² 0.8-0.9 | 3,099,369 (31%) | Good imputation quality |
| Common (MAF ≥ 5%) | ~6.0M (60%) | Standard GWAS variants |
| Low-freq (1-5%) | ~2.0M (20%) | Population-specific signals |
| Rare (0.1-1%) | ~2.0M (20%) | Founder effects, recessive |

---

## Detailed QC Steps

### Step 1: Sample Missingness Filter

**Threshold**: F_MISS > 0.20  
**Removed**: 92 samples  
**Result**: 1,247 → 1,155 samples

### Step 2: IBD-Based Deduplication

**Method**: PI_HAT ≥ 0.98 (genetic duplicates)  
**LD Pruning**: 44,782 independent SNPs  
**Duplicate Clusters**: 49 clusters, 106 individuals  
**Removed**: 57 samples (kept 1 per cluster, lowest F_MISS)  
**Result**: 1,155 → 1,098 samples

### Step 3: SNP QC for Imputation Input

**Filters Applied**:
- Removed I/D allele codes: 494 variants
- Removed duplicate positions: 396 variants
- Fixed Cyrillic homoglyphs: 2 sample IDs (м→m, Х→X)

**Result**: 654,027 → 473,081 variants

### Step 4: Michigan Imputation Server

**Settings**:
- Reference: 1000 Genomes Phase 3 v5 (GRCh38)
- Population: ALL
- Phasing: Eagle v2.4
- Mode: Quality Control + Imputation

**QC Results**:
- Samples passed: 1,098 ✓
- Strand flips: 0 ✓
- Allele switches: 0 ✓
- Palindromic SNPs: 0 ✓

**Output**: 58,886,952 imputed variants

### Step 5: Heterozygosity-Based Sample Filter

**Threshold**: |F| > 0.2  
**Identified**: 24 samples (all F < -0.2, excess heterozygosity)  
**Breakdown**:
- F < -0.3: 12 samples (extreme contamination)
- -0.3 to -0.2: 12 samples (severe excess het)

**Result**: 1,098 → 1,074 samples

### Step 6: Post-Imputation Filtering & Merging

**Process**:
1. Removed 24 het outliers from imputed VCFs
2. Extracted high-quality variants (R² ≥ 0.8, MAF ≥ 0.001)
3. Merged all 22 chromosomes
4. Converted to PLINK format

**Processing Time**: ~7 hours total
- Filtering: ~5.5 hours
- Concatenation: ~26 minutes
- Indexing: ~16 minutes
- PLINK conversion: ~70 minutes

**Final Result**: 10,846,569 high-quality variants

---

## File Manifest

### Primary Analysis Files

```
filtered_clean/
├─ UZB_imputed_HQ_clean.bed (2.8 GB)    ← PRIMARY: Binary genotypes
├─ UZB_imputed_HQ_clean.bim (291 MB)    ← PRIMARY: Variant info
├─ UZB_imputed_HQ_clean.fam (32 KB)     ← PRIMARY: Sample info
├─ UZB_imputed_HQ_clean_ALL.vcf.gz (40 GB) ← BACKUP: VCF format
├─ UZB_imputed_HQ_clean_ALL.vcf.gz.tbi  ← VCF index
└─ chr*.clean.vcf.gz (22 files)         ← Per-chromosome VCFs
```

### Reference Files

```
3_post-imputation/
├─ ConvSK_final_clean.bed/bim/fam      ← Clean typed data (1,074 samples)
├─ het_outliers_remove.txt              ← 24 removed samples
└─ het_outliers_vcf_format.txt          ← VCF format IDs

imputation_results/unz/
├─ HQ_variant_ids.txt                   ← 10M variant ID list
├─ UZB_all.HQ_imputed.R2ge0p8.MAFge0p001.tsv  ← Variant metrics + Uzbek AF
├─ UZB_chr*.AF_DS.tsv                   ← Per-chromosome frequencies
└─ chr*.dose.vcf.gz                     ← Original Michigan outputs (1,098 samples)
```

---

## Data Transformation Summary

| Stage | Samples | Variants | Change |
|-------|---------|----------|--------|
| **Raw Array** | 1,247 | 654,027 | Starting point |
| After Missingness | 1,155 | 654,027 | −92 samples |
| After Deduplication | 1,098 | 654,027 | −57 duplicates |
| After SNP QC | 1,098 | 473,081 | −180K variants |
| **After Imputation** | 1,098 | 58,886,952 | +58M variants |
| After Het Filter | 1,074 | 473,081 (typed) | −24 outliers |
| **FINAL CLEAN** | **1,074** | **10,846,569** | **High-quality set** |

**Net Changes**:
- Samples: 1,247 → 1,074 (−173, 13.9% removed)
- Variants: 654K → 10.8M (23× expansion)
- Genotyping: 98.6% → 100% (+1.4%)

---

## Ready For These Analyses

### Immediate (Phase 2)

- [x] ✅ Data QC Complete
- [ ] LD Pruning (create independent SNP set)
- [ ] PCA (internal population structure)
- [ ] ADMIXTURE (ancestry proportions K=2-10)
- [ ] ROH Analysis (runs of homozygosity)
- [ ] Population-Specific Variants (Uzbek-enriched)

### Advanced (Phase 3)

- [ ] PCA with 1000G (compare to reference populations)
- [ ] FST Calculations (genetic differentiation)
- [ ] IBD Network (cryptic relatedness)
- [ ] GWAS (if case-control phenotype available)
- [ ] Rare Variant Burden Tests

---

## Important Notes

### HWE Filtering

⚠️ **NOT YET APPLIED** - Apply based on analysis type:

**For Population Genetics (PCA, ADMIXTURE, FST):**
```bash
plink --bfile UZB_imputed_HQ_clean --hwe 1e-6 --make-bed --out UZB_for_popgen
```

**For GWAS (Case-Control):**
```bash
# Two-stage filter
plink --bfile UZB_imputed_HQ_clean --filter-controls --hwe 1e-6 --make-bed --out temp
plink --bfile temp --hwe 1e-10 --hwe-all --make-bed --out UZB_for_gwas
```

**For Frequency Analysis:**
- No HWE filter needed (already R² ≥ 0.8 quality)

### Phenotype Data

⚠️ **NOT INCLUDED** - All samples have phenotype = -9 (missing) in .fam file  
- Merge phenotype data separately before case-control analyses
- Pregnancy loss phenotype needs to be joined from clinical database

### Sample ID Format

**Format**: `FID_IID` (e.g., `1_01-01`)  
- All Cyrillic homoglyphs corrected to ASCII
- Consistent across PLINK and VCF formats
- Synchronized with `ConvSK_final_clean`

---

## Quick Start Commands

### Load Dataset

```bash
# In PLINK
plink --bfile /path/to/UZB_imputed_HQ_clean --freq --out test

# Check summary
plink --bfile /path/to/UZB_imputed_HQ_clean --missing
```

### Basic QC Checks

```bash
# Sample count
wc -l UZB_imputed_HQ_clean.fam
# Should be: 1074

# Variant count
wc -l UZB_imputed_HQ_clean.bim
# Should be: 10846569

# Check missingness (should be 0%)
plink --bfile UZB_imputed_HQ_clean --missing
# Genotyping rate should be: 1.0
```

### LD Prune for PCA

```bash
# Create LD-pruned set
plink --bfile UZB_imputed_HQ_clean \
  --maf 0.05 \
  --geno 0.01 \
  --hwe 1e-6 \
  --indep-pairwise 50 5 0.2 \
  --out UZB_ld_prune

# Extract pruned SNPs
plink --bfile UZB_imputed_HQ_clean \
  --extract UZB_ld_prune.prune.in \
  --make-bed \
  --out UZB_pruned
```

---

## Documentation Files

1. **Main Pipeline Report (HTML)**: `alsu_pipeline_report.html`
   - Comprehensive overview with visualizations
   - All QC metrics and results
   - Pre- and post-imputation statistics

2. **QC Visualizations (HTML)**: `alsu_qc_plots.html`
   - Histogram-style plots
   - Sample/variant missingness distributions
   - Heterozygosity rate and F coefficient plots
   - MAF distributions

3. **Step 6 Technical Log (Markdown)**: `step6_technical_log.md`
   - Detailed commands for post-imputation filtering
   - Processing times and file sizes
   - Complete reproducibility record

4. **Quick Reference Card (HTML)**: `final_dataset_card.html`
   - One-page summary
   - Key metrics and file locations
   - Quick start guide

5. **Complete Pipeline Summary (This File)**: `complete_pipeline_summary.md`
   - Text-based comprehensive overview
   - All steps documented
   - Command reference

---

## Reproducibility

### Software Versions

- PLINK: v1.9.0-b.7.7
- BCFtools: v1.17+
- Michigan Imputation Server: 1000G Phase 3 v5 (GRCh38)
- Eagle: v2.4 (phasing)
- Tabix: Latest

### Processing Environment

- Server: Biotech2024
- Working Directory: `/staging/ALSU-analysis/winter2025/`
- Processing Dates: November 30 - December 26, 2025
- Total Processing Time: ~15 hours

### Checksums

To verify dataset integrity:
```bash
md5sum UZB_imputed_HQ_clean.bed
md5sum UZB_imputed_HQ_clean.bim
md5sum UZB_imputed_HQ_clean.fam
```

---

## Citation & Acknowledgments

**Project**: ALSU (Uzbek Pregnancy Loss Genetics Study)  
**Data**: Controlled access via institutional IRB  
**Imputation**: Michigan Imputation Server (Das et al., Nature Genetics 2016)  
**Reference Panel**: 1000 Genomes Project Phase 3 Consortium

---

## Contact

For questions about this dataset or analysis pipeline:
- Project PI: [Contact Information]
- Data Access: Institutional IRB approval required
- Technical Questions: See documentation files listed above

---

**Last Updated**: December 26, 2025  
**Pipeline Status**: ✅ Complete - Ready for Phase 2 Population Genetics Analyses  
**Next Step**: LD Pruning + PCA

---

## Appendix: Key Thresholds Used

| QC Step | Threshold | Rationale |
|---------|-----------|-----------|
| Sample Missingness | F_MISS > 0.20 | Standard GWAS QC |
| IBD Duplication | PI_HAT ≥ 0.98 | Genetic duplicates/twins |
| Heterozygosity Outliers | \|F\| > 0.2 | Literature standard for contamination |
| Imputation Quality | R² ≥ 0.8 | Good confidence for analysis |
| MAF (Uzbek cohort) | MAF ≥ 0.001 | 0.1% frequency threshold |
| HWE (for LD pruning) | p < 1e-6 | Genotyping error detection |
| LD Pruning | r² < 0.2 | Independent SNPs for PCA |

---

**End of Complete Pipeline Summary**
