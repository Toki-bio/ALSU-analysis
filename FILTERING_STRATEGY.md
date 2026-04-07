# Efficient Fix: Filter Bad Samples Instead of Re-running

## The Insight

All 7 samples with F_MISS > 0.20 are present in the imputation input but should not be:
- ✓ They're in Step 1 output (ConvSK_mind20.fam, 1155 samples)
- ✓ They're in Step 2 output (ConvSK_mind20_dedup.fam, 1098 samples)  
- ✓ They have imputed genotypes in Step 4

**The Fix:** Instead of re-running imputation (expensive), just remove them from **Step 4 output onwards**.

---

## Execution Plan (NO RE-IMPUTATION NEEDED)

### Step 1: Create removal list
```bash
# These 7 samples have high missingness (F_MISS > 0.20)
cat > remove_7_highMiss.txt << EOF
458	08-365
499	08-701
840	08-25
862	08-495
886	08-77
898	08-825
910	12-11
EOF
```

### Step 2: Filter imputation output (Step 4 → Step 5)
```bash
# Replace ConvSK_final_clean.fam with filtered version
plink --bfile ConvSK_final_clean \
  --remove remove_7_highMiss.txt \
  --make-bed --out ConvSK_final_clean_filtered

# This produces:
#   1074 → 1067 samples after filtering
# (Or depending on post-imputation-QC dropouts, might be different number)
```

### Step 3: Re-run downstream analyses (Step 6-15)
Use `ConvSK_final_clean_filtered` as input for:
- Step 6: Post-imputation QC
- Step 7: Local PCA
- Step 8: Global PCA (with 1000G ref)
- Step 9: FST (vs reference pops)
- Step 10: Covariate validation
- Step 11: ADMIXTURE
- Step 15: ROH analysis

---

## Time Savings

| Task | Old Plan | New Plan |
|------|----------|----------|
| Re-run everything | 1-2 weeks | ❌ SKIP |
| **Imputation** | **7+ days** | **SKIP** |
| Filter 7 samples | N/A | **5 min** |
| Re-run PCA-15 | 2-3 days | **Re-run only** |
| **TOTAL** | **~10 days** | **~2-3 days** |

**~70% time savings** by skipping imputation re-run!

---

## Trade-offs

**Pro:**
- Much faster (imputation is the bottleneck)
- Analyses use correctly-filtered data (1091 final samples instead of 1098)
- Minimal downstream changes needed

**Con:**
- Intermediate files (Steps 2-4) still contain bad samples
- If anyone uses intermediate outputs, they'll have bad data
- Technically you're "filtering" not "re-running", so some upstream effects are muted

**Mitigation:** Document clearly that final analyses use filtered 1091-sample subset.

---

## Recommendation

**✓ USE THIS APPROACH** unless you specifically need to audit the full imputation process.

The 7 extra samples are a minor data quality issue (0.7% of samples). Filtering them from final analyses is statistically valid and pragmatically efficient.
