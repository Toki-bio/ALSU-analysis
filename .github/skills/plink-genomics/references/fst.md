# FST Analysis Reference

## Full Genome-Wide FST Pipeline

### Prerequisites
- Uzbek/target data in hg19 (liftover if needed — see merge_1000g.md)
- Reference population VCFs (e.g., 1000G Phase 3)
- populations.txt assignment file

### populations.txt format
```
# FID  IID  GROUP  (tab-separated, no header or with header)
sample001  sample001  TARGET
HG00096    HG00096    EUR
HG00097    HG00097    EUR
```

### Per-Chromosome Pipeline (automate with loop)

```bash
for chr in {1..22}; do
  echo "=== CHR ${chr} ==="

  # [1] Extract reference population from 1000G VCF
  plink --vcf 1000G/ALL.chr${chr}.vcf.gz \
    --keep ref_samples.txt \
    --make-bed --out ref_chr${chr}

  # [2] Extract target chromosome
  plink --bfile target_hg19 --chr ${chr} --make-bed --out target_chr${chr}

  # [3] Rename to chr:pos format (enables matching across datasets)
  awk '{$2=$1":"$4; print}' OFS='\t' target_chr${chr}.bim \
    > target_chr${chr}_pos.bim
  cp target_chr${chr}.bed target_chr${chr}_pos.bed
  cp target_chr${chr}.fam target_chr${chr}_pos.fam

  awk '{$2=$1":"$4; print}' OFS='\t' ref_chr${chr}.bim \
    > ref_chr${chr}_pos.bim
  cp ref_chr${chr}.bed ref_chr${chr}_pos.bed
  cp ref_chr${chr}.fam ref_chr${chr}_pos.fam

  # [4] Find overlapping positions
  awk '{print $1":"$4}' target_chr${chr}_pos.bim | sort > target_pos.txt
  awk '{print $1":"$4}' ref_chr${chr}_pos.bim | sort > ref_pos.txt
  comm -12 target_pos.txt ref_pos.txt > positions_chr${chr}.txt
  echo "  Overlapping SNPs: $(wc -l < positions_chr${chr}.txt)"

  # [5] Extract overlapping SNPs
  plink --bfile target_chr${chr}_pos \
    --extract positions_chr${chr}.txt --make-bed --out target_chr${chr}_match
  plink --bfile ref_chr${chr}_pos \
    --extract positions_chr${chr}.txt --make-bed --out ref_chr${chr}_match

  # [6] First merge attempt
  plink --bfile target_chr${chr}_match \
    --bmerge ref_chr${chr}_match \
    --make-bed --out merged_chr${chr}_temp 2>/dev/null

  # [7] Handle strand/allele conflicts
  if [ -f merged_chr${chr}_temp-merge.missnp ]; then
    echo "  Handling merge conflicts..."
    plink --bfile target_chr${chr}_match \
      --exclude merged_chr${chr}_temp-merge.missnp \
      --make-bed --out target_chr${chr}_clean
    plink --bfile ref_chr${chr}_match \
      --exclude merged_chr${chr}_temp-merge.missnp \
      --make-bed --out ref_chr${chr}_clean
    plink --bfile target_chr${chr}_clean \
      --bmerge ref_chr${chr}_clean \
      --make-bed --out merged_chr${chr}
  else
    mv merged_chr${chr}_temp.bed merged_chr${chr}.bed
    mv merged_chr${chr}_temp.bim merged_chr${chr}.bim
    mv merged_chr${chr}_temp.fam merged_chr${chr}.fam
  fi

  snps=$(wc -l < merged_chr${chr}.bim)
  echo "  ✓ Chr${chr}: ${snps} SNPs merged"

  # Cleanup intermediates
  rm -f target_chr${chr}.* ref_chr${chr}.* \
        target_chr${chr}_pos.* ref_chr${chr}_pos.* \
        target_chr${chr}_match.* ref_chr${chr}_match.* \
        target_chr${chr}_clean.* ref_chr${chr}_clean.* \
        target_pos.txt ref_pos.txt positions_chr${chr}.txt \
        merged_chr${chr}_temp.*
done
```

### Merge All Chromosomes

```bash
# Create merge list (chr2-22)
for chr in {2..22}; do echo "merged_chr${chr}"; done > merge_list.txt

plink --bfile merged_chr1 \
  --merge-list merge_list.txt \
  --make-bed --out merged_all
```

### Calculate FST

```bash
plink --bfile merged_all \
  --fst \
  --within populations.txt \
  --out fst_results

# Output: fst_results.fst
# Columns: CHR  SNP  POS  NMISS  FST
```

## FST Manhattan Plot

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Load FST results
fst = pd.read_csv('fst_results.fst', sep='\s+')
fst['FST'] = fst['FST'].clip(lower=0)  # treat negative as 0

# Cumulative position for Manhattan layout
chromosomes = fst['CHR'].unique()
chrom_offsets = {}
offset = 0
for chrom in sorted(chromosomes):
    chrom_offsets[chrom] = offset
    offset += fst[fst['CHR'] == chrom]['POS'].max()

fst['cum_pos'] = fst.apply(
    lambda r: r['POS'] + chrom_offsets[r['CHR']], axis=1)

# Alternating colors
colors = ['#2166AC', '#4DAC26']
fst['color'] = fst['CHR'].apply(lambda c: colors[c % 2])

fig, ax = plt.subplots(figsize=(16, 5))
ax.scatter(fst['cum_pos'], fst['FST'],
           c=fst['color'], s=0.5, alpha=0.5, linewidths=0)

# Add threshold lines
ax.axhline(0.10, color='orange', linestyle='--', linewidth=0.8, label='FST=0.10')
ax.axhline(0.30, color='red', linestyle='--', linewidth=0.8, label='FST=0.30')

# Chromosome labels
xticks = [fst[fst['CHR']==c]['cum_pos'].mean() for c in sorted(chromosomes)]
ax.set_xticks(xticks)
ax.set_xticklabels([str(c) for c in sorted(chromosomes)], fontsize=7)
ax.set_xlabel('Chromosome')
ax.set_ylabel('FST')
ax.set_title('Genome-Wide FST — {POP1} vs {POP2}')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig('fst_manhattan.png', dpi=300, bbox_inches='tight')
```

## Top Differentiated Loci Annotation

```python
# Extract top FST hits for annotation
fst_clean = fst[fst['FST'] > 0].copy()
top_hits = fst_clean.nlargest(500, 'FST')
top_hits.to_csv('top_fst_hits.txt', sep='\t', index=False)

# Summary statistics
print(f"Total SNPs: {len(fst)}")
print(f"FST > 0.10: {(fst['FST'] > 0.10).sum()}")
print(f"FST > 0.30: {(fst['FST'] > 0.30).sum()}")
print(f"FST > 0.50: {(fst['FST'] > 0.50).sum()}")
print(f"Mean FST: {fst_clean['FST'].mean():.4f}")
print(f"Weighted FST: reported in PLINK log")
```

## ALSU Reference Results (Uzbek vs EUR)

- Samples: 1,199 Uzbek + 503 EUR (1000G: CEU/GBR/FIN/TSI/IBS)
- Overlapping SNPs: 376,208
- Mean FST: 0.0160
- Weighted FST: **0.0204**
- FST > 0.50: 256 loci
- Chr6 highest per-chromosome mean (0.0192) — HLA/MHC region
