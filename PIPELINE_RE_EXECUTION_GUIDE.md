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

Only Steps 6-15 need to be re-run (using filtered 1,091-sample input):

| Step | Current Input | Filtered Input | Action |
|------|---|---|---|
| 1-5 | 1,098 samples (with 7 bad) | Keep as is | **NO CHANGE** |
| 6 | 1,074 after QC | ~1,067 after filtering | **RE-RUN** |
| 7 | ~1,062 after PCA | ~1,055 after filtering | **RE-RUN** |
| 8-15 | All based on buggy input | All corrected output | **RE-RUN** |

---

## Impact Assessment

### Why Filtering Without Re-imputation Works

The 7 samples with high missingness (F_MISS > 0.20) problem:
- **Original genotypes:** Poor quality (that's why F_MISS > 0.20)
- **Imputed genotypes:** Still valid (imputation computes them from reference panel, independent of original quality)
- **Solution:** Use imputed genotypes, just exclude samples from final analyses

It's not like the imputation was "wrong" — it was computed correctly. We're just deciding to exclude the sample afterwards.

### What Gets Re-run

- **High Impact (Re-run needed):**
  - **PCA (Steps 7-8):** Different outlier cutoffs with 1,091 vs 1,098 samples
  - **ADMIXTURE/FST/Ancestry** (Steps 9, 11): Population structure changes with corrected sample set
  
- **No Impact (No re-run):**
  - **Imputation (Step 4):** ✓ SKIP — imputed genotypes are valid, just filter samples
  - **SNP-level QC (Step 3):** Already done; filtering samples doesn't change variant statistics
  - **IBD (Step 2):** Already done; filtering samples afterwards is OK

### Bottom Line

Only downstream population-level analyses (PCA, ADMIXTURE, FST) need re-running because they're sensitive to sample composition. Imputation and QC are indifferent.

---

## Which Analyses Are Currently Wrong?

**Only downstream population-level analyses** are sensitive to having those 7 bad samples:

- ❌ ADMIXTURE results (based on 1,098 instead of 1,091)
- ❌ FST/PBS statistics (reference pops vs 1,098 UZB vs 1,091)  
- ❌ PCA plots (local + global, with different outlier cutoffs)

**NOT affected** (no re-run needed):
- ✓ Imputation results (valid regardless of sample quality)
- ✓ Variant QC statistics (HWE, MAF, info scores)
- ✓ IBD results (already computed; filtering samples afterwards is OK)
- ✓ Individual-level data (genotypes, phenotypes, covariates)

---

## Decision: Should We Re-run?

**✓ RECOMMENDED: Filter bad samples, re-run only Steps 6-15**

See `FILTERING_STRATEGY.md` for detailed explanation.

**Why:**
- Imputation doesn't need re-run (valid regardless of original sample quality)
- Only population-level analyses (PCA, ADMIXTURE) are sensitive to sample composition
- Time: ~2-3 days (vs. 10+ days for full re-run)
- Result: Final analyses use correct 1,091-sample set

**Alternative** (if full re-run wanted):


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

## Next Steps (Filtering Strategy)

1. **Filter the post-imputation data**
   ```bash
   plink --bfile ConvSK_final_clean \
     --remove remove_7_highMiss.txt \
     --make-bed --out ConvSK_final_clean_filtered
   ```

2. **Re-run Steps 6-15 using filtered input**
   - Step 6: Post-imputation QC
   - Step 7: Local PCA  
   - Step 8: Global PCA (with 1000G)
   - Step 9: FST analysis
   - Step 10: Covariate validation
   - Step 11: ADMIXTURE
   - Step 15: ROH analysis

3. **Update published results** with corrected figures (1,091 samples instead of 1,098)

4. **Document in reproducibility notes** that final analyses use 1,091-sample subset after excluding high-missingness individuals

---

## Validation Checklist

- [ ] Filter operation produces 1,067 samples (1,074 - 7)
- [ ] PCA plots updated with correct sample count
- [ ] ADMIXTURE results recalculated  
- [ ] FST statistics re-computed
- [ ] All figures/tables updated
