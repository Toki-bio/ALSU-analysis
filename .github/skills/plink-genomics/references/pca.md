# PCA Visualization Reference

## R Code (ggplot2) — Publication Quality

```r
library(ggplot2)
library(dplyr)

# Load PCA results
eigenvec <- read.table("{PREFIX}_pca.eigenvec", header = TRUE, comment.char = "")
eigenval <- read.table("{PREFIX}_pca.eigenval", header = FALSE)
colnames(eigenvec)[1:2] <- c("FID", "IID")

# Calculate % variance explained
pve <- data.frame(
  PC = 1:nrow(eigenval),
  pve = (eigenval$V1 / sum(eigenval$V1)) * 100
)

# Basic plot (single population, no color groups)
pca_plot <- ggplot(eigenvec, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.6, color = "darkblue", size = 1.5) +
  theme_bw() +
  labs(
    title = "PCA — {COHORT_NAME}",
    subtitle = paste0("N=", nrow(eigenvec), " samples, ", 
                      "{N_SNPS} independent SNPs"),
    x = paste0("PC1 (", round(pve$pve[1], 2), "%)"),
    y = paste0("PC2 (", round(pve$pve[2], 2), "%)")
  ) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("{PREFIX}_pca.png", pca_plot, width = 8, height = 6, dpi = 300)
ggsave("{PREFIX}_pca.pdf", pca_plot, width = 8, height = 6)
cat("PC1:", round(pve$pve[1],2), "%\nPC2:", round(pve$pve[2],2), "%\n")
```

## R Code — Multi-population (with 1000G reference)

```r
library(ggplot2)
library(RColorBrewer)

eigenvec <- read.table("{PREFIX}_pca.eigenvec", header = TRUE)
eigenval <- read.table("{PREFIX}_pca.eigenval", header = FALSE)
colnames(eigenvec)[1:2] <- c("FID", "IID")
pve <- eigenval$V1 / sum(eigenval$V1) * 100

# Load population labels
# pop_labels.txt: FID  IID  POP  SUPERPOP (tab-separated)
pop_labels <- read.table("pop_labels.txt", header = TRUE)
eigenvec <- merge(eigenvec, pop_labels, by = c("FID","IID"))

# Color palette for superpopulations
superpop_colors <- c(
  "AFR" = "#E41A1C",
  "AMR" = "#FF7F00",
  "EAS" = "#4DAF4A",
  "EUR" = "#377EB8",
  "SAS" = "#984EA3",
  "UZB" = "#F7DC6F"   # Your cohort — yellow/gold for visibility
)

# Plot — reference populations small, your cohort larger
eigenvec$point_size <- ifelse(eigenvec$SUPERPOP == "UZB", 2.0, 0.8)
eigenvec$point_alpha <- ifelse(eigenvec$SUPERPOP == "UZB", 0.8, 0.4)

p <- ggplot(eigenvec, aes(x = PC1, y = PC2, color = SUPERPOP)) +
  geom_point(aes(size = point_size, alpha = point_alpha)) +
  scale_color_manual(values = superpop_colors) +
  scale_size_identity() +
  scale_alpha_identity() +
  theme_bw() +
  labs(
    title = "Global PCA with 1000 Genomes Reference",
    x = paste0("PC1 (", round(pve[1], 2), "%)"),
    y = paste0("PC2 (", round(pve[2], 2), "%)"),
    color = "Population"
  ) +
  theme(legend.position = "right")

ggsave("{PREFIX}_global_pca.png", p, width = 10, height = 7, dpi = 300)
```

## Python Code (matplotlib/seaborn)

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Load data
eigenvec = pd.read_csv('{PREFIX}_pca.eigenvec', sep='\s+',
                        header=0, comment='#')
eigenval = pd.read_csv('{PREFIX}_pca.eigenval', header=None)[0]
eigenvec.columns = ['FID','IID'] + [f'PC{i}' for i in range(1, len(eigenvec.columns)-1)]

pve = eigenval / eigenval.sum() * 100

# Single-population plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(eigenvec['PC1'], eigenvec['PC2'],
           alpha=0.6, s=15, color='steelblue', linewidths=0)
ax.set_xlabel(f'PC1 ({pve[0]:.2f}%)', fontsize=12)
ax.set_ylabel(f'PC2 ({pve[1]:.2f}%)', fontsize=12)
ax.set_title('{COHORT_NAME} — Principal Component Analysis', fontsize=13)
plt.tight_layout()
plt.savefig('{PREFIX}_pca.png', dpi=300, bbox_inches='tight')
plt.savefig('{PREFIX}_pca.pdf', bbox_inches='tight')
print(f"PC1: {pve[0]:.2f}%, PC2: {pve[1]:.2f}%")
```

## Scree Plot (variance explained per PC)

```python
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(range(1, len(pve)+1), pve, color='steelblue', alpha=0.8)
ax.set_xlabel('Principal Component')
ax.set_ylabel('Variance Explained (%)')
ax.set_title('PCA Scree Plot')
ax.set_xticks(range(1, len(pve)+1))
plt.tight_layout()
plt.savefig('{PREFIX}_scree.png', dpi=150)
```

## Outlier Detection and Removal

Samples > 6 SD from mean on any PC are likely population outliers:

```python
def flag_pca_outliers(eigenvec, n_pcs=10, sd_threshold=6):
    pc_cols = [f'PC{i}' for i in range(1, n_pcs+1)]
    outlier_mask = pd.Series(False, index=eigenvec.index)
    for pc in pc_cols:
        mean, sd = eigenvec[pc].mean(), eigenvec[pc].std()
        outlier_mask |= (np.abs(eigenvec[pc] - mean) > sd_threshold * sd)
    outliers = eigenvec[outlier_mask][['FID','IID']]
    print(f"Flagged {len(outliers)} PCA outliers")
    return outliers

outliers = flag_pca_outliers(eigenvec)
outliers.to_csv('pca_outliers_remove.txt', sep='\t', index=False, header=False)
```

## ALSU Reference Values

- Dataset: 1,062 Uzbek samples, 88,315 independent SNPs (r²<0.05)
- PC1: 28.93% variance, PC2: 11.51% — strong internal structure
- Interpretation: geographic cline West–East within Uzbekistan
- PC1+PC2 used as GWAS covariates alongside ADMIXTURE Q-values
