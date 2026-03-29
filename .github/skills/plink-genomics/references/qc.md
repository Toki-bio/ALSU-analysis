# QC Parameter Reference

## Sample QC Decision Flowchart

```
START: What type of data?
│
├── Genotyping array (raw, pre-imputation)
│   ├── Missingness (--mind): 0.20
│   ├── SNP missing rate (--geno): 0.01
│   ├── MAF (--maf): 0.01
│   └── HWE (--hwe): 1e-6
│
├── Post-imputation (Michigan/TOPMed output)
│   ├── Missingness (--mind): 0.05
│   ├── SNP missing rate (--geno): 0.05
│   ├── MAF (--maf): 0.01 (use 0.001 for rare variant studies)
│   └── HWE (--hwe): 1e-6
│
└── GWAS-ready final dataset
    ├── Missingness (--mind): 0.02
    ├── SNP missing rate (--geno): 0.02
    ├── MAF (--maf): 0.05 (common variants) or 0.01
    └── HWE (--hwe): 1e-6 (controls only if case-control design)
```

## HWE — Special Cases

**Case-control studies:** Apply HWE filter to controls only, not cases
(disease may cause HWE deviation in cases):
```bash
# Filter HWE in controls only
plink --bfile {PREFIX} \
  --hwe 1e-6 --hwe-all \
  --keep controls.txt \
  --make-bed --out {PREFIX}_hwe_controls

# Extract passing SNPs, apply to full dataset
plink --bfile {PREFIX} \
  --extract {PREFIX}_hwe_controls_passing_snps.txt \
  --make-bed --out {PREFIX}_hwe
```

**Regions to exclude from HWE:** MHC/HLA region (chr6:28-34 Mb in hg19) often shows
HWE deviation due to strong selection — consider excluding from HWE filter:
```bash
plink --bfile {PREFIX} \
  --hwe 1e-6 \
  --exclude-range 6:28000000-34000000 \
  --make-bed --out {PREFIX}_hwe
```

## Missingness Histogram (inspect before choosing threshold)

```python
import pandas as pd
import matplotlib.pyplot as plt

imiss = pd.read_csv('{PREFIX}_miss.imiss', sep='\s+')
plt.figure(figsize=(8,4))
plt.hist(imiss['F_MISS'], bins=50, color='steelblue', edgecolor='white')
plt.axvline(0.05, color='orange', linestyle='--', label='0.05 threshold')
plt.axvline(0.20, color='red', linestyle='--', label='0.20 threshold')
plt.xlabel('Per-sample missingness rate (F_MISS)')
plt.ylabel('Number of samples')
plt.title('Sample Missingness Distribution')
plt.legend()
plt.tight_layout()
plt.savefig('missingness_hist.png', dpi=150)
```

## Expected QC Metrics (ALSU cohort as reference)

| Step | Before | After | Notes |
|---|---|---|---|
| Raw array | 1,247 samples, 654,027 SNPs | — | ConvSK_raw |
| Missingness (mind=0.20) | 1,247 | 1,155 (−92, 7.4%) | ConvSK_mind20 |
| IBD dedup (PI_HAT≥0.98) | 1,155 | 1,098 (−57, 4.9%) | ConvSK_mind20_dedup |
| SNP QC (maf=0.01, hwe=1e-6, geno=0.01) | 654,027 SNPs | 472,463 SNPs (−28%) | For imputation |
| Post-imputation QC (mind=0.05, geno=0.05, maf=0.01) | 1,074 / 10.8M | 1,062 / 5.38M | UZB_imputed_HQ_qc |

These are reference values only — your cohort will differ.

## Sex Check (optional but recommended)

```bash
plink --bfile {PREFIX} --check-sex --out {PREFIX}_sexcheck
# Review: {PREFIX}_sexcheck.sexcheck
# STATUS column: OK = pass, PROBLEM = discordant
awk '$5=="PROBLEM" {print $1, $2}' {PREFIX}_sexcheck.sexcheck > sex_problems.txt
```

## Heterozygosity Outliers

Samples with extreme heterozygosity (> 3 SD from mean) may indicate contamination or inbreeding:

```python
import pandas as pd
import numpy as np

het = pd.read_csv('{PREFIX}.het', sep='\s+')
het['F_HET'] = (het['O(HOM)'] - het['E(HOM)']) / het['N(NM)']
mean_het = het['F_HET'].mean()
sd_het = het['F_HET'].std()

# Flag outliers (>3 SD)
outliers = het[np.abs(het['F_HET'] - mean_het) > 3 * sd_het]
outliers[['FID','IID']].to_csv('het_outliers.txt', sep='\t', index=False)
print(f"Flagged {len(outliers)} heterozygosity outliers")
```

```bash
plink --bfile {PREFIX} --het --out {PREFIX}_het
```
