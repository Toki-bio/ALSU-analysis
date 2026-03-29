# ADMIXTURE Analysis Reference

## Full Pipeline

```bash
# 1. LD pruning (r²=0.3 — more permissive than PCA to retain more SNPs)
plink --bfile {PREFIX} \
  --indep-pairwise 50 5 0.3 \
  --out {PREFIX}_prune_admix

plink --bfile {PREFIX} \
  --extract {PREFIX}_prune_admix.prune.in \
  --make-bed --out {PREFIX}_for_admixture

# 2. Run ADMIXTURE K=2 to K=8 with cross-validation
for K in {2..8}; do
  admixture --cv -j{THREADS} {PREFIX}_for_admixture.bed $K \
    | tee admix_log_K${K}.out
done

# 3. Extract CV errors
grep "CV error" admix_log_K*.out
```

## K Selection: Cross-Validation Error

```python
import re
import pandas as pd
import matplotlib.pyplot as plt

cv_errors = {}
for k in range(2, 9):
    with open(f'admix_log_K{k}.out') as f:
        for line in f:
            m = re.search(r'CV error \(K=(\d+)\): ([0-9.]+)', line)
            if m:
                cv_errors[int(m.group(1))] = float(m.group(2))

cv_df = pd.DataFrame(list(cv_errors.items()), columns=['K', 'CV_error'])
cv_df = cv_df.sort_values('K')

optimal_k = cv_df.loc[cv_df['CV_error'].idxmin(), 'K']
print(f"Optimal K (minimum CV error): {optimal_k}")
print(cv_df.to_string(index=False))

plt.figure(figsize=(7,4))
plt.plot(cv_df['K'], cv_df['CV_error'], 'o-', color='steelblue')
plt.axvline(optimal_k, color='red', linestyle='--', label=f'Optimal K={optimal_k}')
plt.xlabel('K (number of ancestral populations)')
plt.ylabel('Cross-validation error')
plt.title('ADMIXTURE Model Selection')
plt.legend()
plt.tight_layout()
plt.savefig('admixture_cv_error.png', dpi=150)
```

## sNMF Validation (R / LEA package)

```r
# Install: BiocManager::install("LEA")
library(LEA)

# Convert to .geno format
# Requires VCF input — convert from PLINK first:
# plink --bfile {PREFIX}_for_admixture --recode vcf --out {PREFIX}_for_snmf
geno_file <- vcf2geno("{PREFIX}_for_snmf.vcf")

# Run sNMF K=2-10, 10 replicates
obj <- snmf(geno_file,
  K = 2:10,
  repetitions = 10,
  entropy = TRUE,   # compute cross-entropy for K selection
  alpha = 10,       # L2 regularization
  iterations = 200,
  tolerance = 1e-5,
  project = "new",
  CPU = {THREADS})

# Plot cross-entropy
plot(obj, cex = 1.2, col = "lightblue", pch = 19)

# Get optimal K
ce_means <- sapply(2:10, function(k) mean(cross.entropy(obj, K=k)))
optimal_k <- which.min(ce_means) + 1
cat("Optimal K (min cross-entropy):", optimal_k, "\n")

# Extract Q-matrix at optimal K (best run)
best_run <- which.min(cross.entropy(obj, K=optimal_k))
Q_matrix <- Q(obj, K=optimal_k, run=best_run)
```

## Stacked Bar Plot (publication quality)

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Load Q matrix
K = {OPTIMAL_K}
Q = np.loadtxt(f'{PREFIX}_for_admixture.{K}.Q')
fam = pd.read_csv(f'{PREFIX}_for_admixture.fam', sep='\s+',
                   header=None, names=['FID','IID','PAT','MAT','SEX','PHE'])

# Load population labels
pop_labels = pd.read_csv('pop_labels.txt', sep='\t')  # FID IID POP SUPERPOP
data = pd.DataFrame(Q, columns=[f'Q{i+1}' for i in range(K)])
data['IID'] = fam['IID'].values
data = data.merge(pop_labels[['IID','POP','SUPERPOP']], on='IID', how='left')

# Sort by population, then by dominant component
data = data.sort_values(['SUPERPOP','POP'] + [f'Q{i+1}' for i in range(K)],
                         ascending=[True, True] + [False]*K)

# Color palette — adjust to match your K
colors = ['#2166AC','#F4A261','#2DC653','#E41A1C',
          '#984EA3','#FF7F00','#A65628','#F781BF'][:K]

fig, ax = plt.subplots(figsize=(16, 3))
bottoms = np.zeros(len(data))
for i, (col, color) in enumerate(zip([f'Q{j+1}' for j in range(K)], colors)):
    ax.bar(range(len(data)), data[col], bottom=bottoms,
           color=color, width=1.0, linewidth=0)
    bottoms += data[col].values

# Population dividers and labels
pops = data['POP'].values
pop_starts = [0] + [i+1 for i in range(len(pops)-1) if pops[i] != pops[i+1]]
pop_ends = pop_starts[1:] + [len(pops)]
for s in pop_starts[1:]:
    ax.axvline(s - 0.5, color='white', linewidth=1.0)

for start, end, pop in zip(pop_starts, pop_ends, data.iloc[pop_starts]['POP']):
    ax.text((start + end) / 2, -0.08, pop,
            ha='center', va='top', fontsize=8, rotation=45,
            transform=ax.get_xaxis_transform())

ax.set_xlim(-0.5, len(data) - 0.5)
ax.set_ylim(0, 1)
ax.set_yticks([0, 0.5, 1.0])
ax.set_xticks([])
ax.set_ylabel('Ancestry proportion')
ax.set_title(f'ADMIXTURE K={K}')
plt.tight_layout()
plt.savefig(f'admixture_K{K}_barplot.png', dpi=300, bbox_inches='tight')
plt.savefig(f'admixture_K{K}_barplot.pdf', bbox_inches='tight')
```

## ALSU Reference Results

**Global (Uzbek + 1000G, 2,095 samples × 380,376 SNPs):**
- Optimal K = 7 (CV error minimum)
- Uzbek ancestry at K=7: 70.2% Central Asian, 12.3% East Asian,
  4.5% N.European, 4.4% C.Asian (minor), 4.3% S.European, 3.9% S.Asian, 0.4% African
- Evanno ΔK peaks at K=3 (continental split) — use CV error, not ΔK

**Uzbek-only (1,074 samples × 380,376 SNPs):**
- Optimal K = 2 (both CV error and sNMF cross-entropy)
- K=2 reflects EUR↔EAS admixture cline within cohort
- Geographic, not ethnic: birthplace significant (p~10⁻²²), ethnicity non-significant

## Key Statistical Notes

- **Do NOT use Evanno ΔK for ADMIXTURE** — designed for STRUCTURE (MCMC), not ADMIXTURE (ML)
- Validate with sNMF independently — concordance r>0.99 confirms robust results
- For GWAS: use K=optimal Q-values as covariates, NOT binary ethnic labels
- Clinal populations (Central Asian, admixed): no single "true K" exists —
  K selection identifies most parsimonious description, not historical truth
