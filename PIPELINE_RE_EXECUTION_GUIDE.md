# Pipeline Re-execution Requirements

## Summary

**The original Step 1 sample filtering is BUGGY.** The documented value (92 removed → 1,155 retained) is incorrect. The correct value is **99 removed → 1,148 retained**.

All downstream analyses (Steps 2-15) were executed using the buggy 1,155-sample input and must be re-executed using the correct 1,148-sample input.

---

## What's Wrong

| Issue | Details |
|-------|---------|
| **Error** | Step 1 used wrong awk filter, missed 7 samples with F_MISS > 0.20 |
| **Samples affected** | 458/08-365, 499/08-701, 840/08-25, 862/08-495, 886/08-77, 898/08-825, 910/12-11 |
| **Current data in use** | ConvSK_mind20.fam: **1,155 samples** (should be 1,148) |
| **Correct data exists** | `/staging/ALSU-analysis/winter2025/PLINK_301125_0312/old/ConvSK_mind20.fam` (1,148 samples) |

---

## Steps That Need Re-execution

| Step | Current (Buggy) | Correct | Input Required | Action |
|------|---|---|---|---|
| 1 | 1,247 → **1,155** | 1,247 → **1,148** | N/A | ✓ FIXED (doc only) |
| 2 | **1,155** → 1,098 | **1,148** → 1,091 | 1,148 input | **RE-RUN** |
| 3 | **1,098** (SNP-QC) | **1,091** (SNP-QC) | 1,091 input | **RE-RUN** |
| 4 | **1,098** (imputation input) | **1,091** (imputation input) | 1,091 input | **RE-RUN IMPUTATION** |
| 5 | **1,098** (ID norm) | **1,091** (ID norm) | imputed 1,091 | **RE-RUN** |
| 6 | **1,074** (post-QC) | ~**1,067** (post-QC) | imputed 1,091 | **RE-RUN** |
| 7 | **~1,062** (post-PCA) | ~**1,055** (post-PCA) | post-QC 1,067 | **RE-RUN** |
| 8-15 | **All buggy** | **All correct** | corrected upstream | **RE-RUN ALL** |

---

## Impact Assessment

### High Impact
- **Imputation (Step 4)**: Using 1,098 instead of 1,091 samples impacts variant density estimates and INFO scores. Re-imputation required.
- **ADMIXTURE / FST / Ancestry** (Steps 9, 11): Population structure estimates assume 1,098 samples; different with 1,091.
- **PCA (Steps 7-8)**: Outlier detection and projection change with different sample size.

### Medium Impact
- **IBD (Step 2)**: 57 samples removed is constant, but pair count changes (666,585 → 657,778).
- **SNP-level QC (Step 3)**: Variant filtering (HWE, MAF) may differ slightly with different sample size.

### Low Impact
- **ID normalization (Step 5)**: No sample count change; 2 samples fixed regardless.
- **Data merging**: File formats unchanged; imputation merging logic unchanged.

---

## Which Analyses Are Currently Wrong?

Everything in the documentation from Step 2 onwards is based on BUGGY data:

- ❌ ADMIXTURE results (1,047 samples, but should be 1,040)
- ❌ FST/PBS statistics (reference pops vs 1,047 UZB)  
- ❌ PCA plots (global + local)
- ❌ GWAS results (if any, using 1,062 in Step 8)
- ❌ IBD estimates
- ❌ ROH statistics
- ❌ Variant INFO scores from imputation
- ❌ Allele frequency distributions

---

## Decision: Should We Re-run?

**RECOMMENDATION: NO FULL RE-RUN NEEDED** ✓ Use filtering strategy instead

See `FILTERING_STRATEGY.md` for the efficient approach:
- Keep imputation results as-is (expensive)
- Just remove 7 bad samples from Step 4 output  
- Re-run only Steps 6-15 (fast downstream analyses)
- **Result:** Final analyses use correct 1,091-sample set
- **Time savings:** ~70% (skip 7-day imputation)

**Original Plan (if full re-run is required):**


1. Seven samples carrying data quality issues are included in current analyses
2. Results differ from ground truth; publications using this data are incorrect
3. The corrected 1,148-sample version already exists and is ready to use
4. Downstream steps are deterministic; re-running will update all dependent analyses

**Cost:** Medium (imputation is expensive, rest is fast)  
**Benefit:** High (accuracy, reproducibility, publication validity)

---

## Server Status

**Location:** `/staging/ALSU-analysis/winter2025/PLINK_301125_0312/`

**Available files:**
- **Buggy version (1,155 samples):**
  - ConvSK_mind20.fam
  - ConvSK_mind20_dedup.fam (1,098 samples)
  - Downstream: ConvSK_final_clean.fam (1,074 samples in post-imputation)

- **Corrected version (1,148 samples):**
  - old/ConvSK_mind20.fam ← **USE THIS AS INPUT**
  - old/ConvSK_mind20_dedup.fam (1,091 samples)
  - _NOT YET CREATED:_ Post-imputation derivatives

---

## Next Steps

1. **Backup current analyses** (save 1,098-sample results for comparison)
2. **Start Step 2 with 1,148-sample input**
3. **Re-execute all downstream steps**
4. **Compare results** (should differ by ~0.7% in sample size only for most metrics)
5. **Update publication/results** with corrected figures

---

## References

See these files for detailed evidence:
- `old_logs/ConvSK_QC_and_VCF_report.md` — Correction note with proof
- `steps/step1.html` — Warning banner about buggy production data
- `steps/step2.html` — Warning banner and updated cascade numbers
