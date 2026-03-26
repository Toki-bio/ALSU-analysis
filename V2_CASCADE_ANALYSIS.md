# ALSU Pipeline Fork v2 — Cascade Analysis

## Root Cause
**Step 6.4** used `bcftools view --include 'ID=@HQ_variant_ids.txt'` to select HQ-imputed variants.
The file contained 8,634 entries with ID=`.`, which matched ALL 813,191 VCF records having ID=`.`.
Result: 10,846,569 variants instead of the intended 10,009,530 (+837,039 spurious).

## Fix
Replace ID-based filtering with position-based filtering:
```
bcftools view --targets-file HQ_targets_posrefalt.tsv  (CHROM, POS, REF, ALT)
```

## Cascade Impact by Step

| Step | Description | Affected? | Action Required | Notes |
|------|-------------|-----------|-----------------|-------|
| 1-5 | Sample QC, VCF prep, imputation | NO | None | Pre-imputation, unrelated to ID filtering |
| **6** | Post-imputation HQ filtering | **YES — ROOT** | **Re-run with --targets-file** | recalc_v2.sh Phase 1-2 |
| **7** | QC + LD pruning | **YES** | **Re-run QC on v2 data** | recalc_v2.sh Phase 4-6 |
| **8** | PCA | **YES** | Re-compute PCA on v2 LD-pruned set | validate_v2.sh does quick PCA |
| 9 | Global PCA (1000G merge) | **MAYBE** | Uses LD-pruned set → may change slightly | Low priority — re-run for consistency |
| **10** | ADMIXTURE | **YES** | Re-run K=2..10 on v2 LD-pruned set | Significant: different variant set |
| 11 | ROH | **YES** | Re-run on v2 QC set (UZB_v2_qc) | Uses full QC dataset, not LD-pruned |
| **12** | IBD | **YES** | Re-run on v2 QC set | Same as ROH |
| 13 | PBS/FST | NO | Uses pre-QC 1,199-sample dataset | Independent data path |
| 14 | Tajima's D | NO | Same as PBS — uses pre-QC data | Independent data path |
| **15** | UZB-specific SNP analysis | **YES** | Re-run on v2 QC set | Uses imputed QC variants |

## Expected Impact Sizes

### Variant counts
- **Pre-QC**: 10,846,569 → ~10,009,530 (drop 837K spurious)
- **Post-QC**: 5,383,832 → **TBD** (fewer removed by --geno because spurious variants inflated missingness)
- **LD-pruned**: **TBD** (fewer independent variants without noise)

### Sample impact
- V1 removed 12 samples at --mind 0.05
- V2 may rescue some/all of those samples (spurious variants caused artificial missingness)
- Expected: 1,062 → 1,062..1,074 samples post-QC

### PCA / ADMIXTURE
- PC loadings will change slightly
- ADMIXTURE proportions will change
- Not necessarily "wrong" in V1 — QC removed most spurious variants
- But cleaner data = more trustworthy results

## Script Files
| Script | Purpose | Runtime |
|--------|---------|---------|
| `recalc_v2.sh` | Full recalculation pipeline (Phases 0-7) | ~6-8 hours |
| `validate_v2.sh` | Post-recalculation comparison & PCA | ~30 min |

## Deployment Checklist
1. [ ] Upload `recalc_v2.sh` to DRAGEN server
2. [ ] Run with `nohup bash recalc_v2.sh > recalc_v2.log 2>&1 &`
3. [ ] Monitor: `tail -f recalc_v2.log`
4. [ ] After completion, run `validate_v2.sh`
5. [ ] Download v2 PCA eigenvec for local visualization
6. [ ] Re-run ADMIXTURE separately (K=2..10, ~12-24 hours)
7. [ ] Re-run ROH/IBD on v2_qc dataset
8. [ ] Update step documentation in HTML files
