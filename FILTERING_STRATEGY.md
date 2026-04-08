# Efficient Fix: Filter Bad Samples (No Re-runs Needed)

## The Key Insight

The 7 samples with high missingness (F_MISS > 0.20) still have **valid imputation results**. 

Imputation doesn't care about original genotyping quality — it just computes genotype probabilities using the reference panel. Even if a sample had poor original genotyping, its imputed genotypes are statistically sound.

**Therefore: We don't need to re-run imputation or any upstream steps.** Just filter the 7 samples from final analyses.

---

## The One-Step Fix

```bash
# Create filter file with the 7 bad samples
cat > remove_7_highMiss.txt << EOF
458	08-365
499	08-701
840	08-25
862	08-495
886	08-77
898	08-825
910	12-11
EOF

# Filter the post-imputation results
plink --bfile ConvSK_final_clean \
  --remove remove_7_highMiss.txt \
  --make-bed --out ConvSK_final_clean_filtered

# Now re-run Steps 6-15 (PCA, ADMIXTURE, FST, ROH, etc.) 
# using ConvSK_final_clean_filtered as input
# Result: 1074 → 1067 samples (or whatever post-QC gives)
```

---

## Why This Works

| Item | Status |
|------|--------|
| Original genotypes | Keep bad samples (don't care about their original quality) |
| Imputed genotypes | Keep bad samples (imputation results are valid) |
| Final analyses | **Filter out the 7** (we admit original data was bad) |
| Dataset used for results | 1,091 samples ✓ |

---

## Time Required

- Filter operation: **5 minutes**
- Re-run Steps 6-15: **2-3 days** (PCA, ADMIXTURE, etc.)
- **TOTAL: ~2-3 days**
- **NO imputation re-run needed** (saves 7+ days)

---

## Why Not Re-run Imputation?

Imputation is a **reference-based computation**. It says: "Given these reference haplotypes and this sample's genotypes at common SNPs, what are the most likely genotypes at rare variants?"

The answer doesn't change if we later decide "we don't like this sample." The imputed genotypes are mathematically valid regardless.

It's like saying: "We computed your tax return incorrectly because we included an expense we now think shouldn't count." The solution isn't to recompute taxes — it's to just exclude that item.

---

## Bottom Line

**No re-runs needed.** Just filter and proceed with downstream analyses.
