# Merging with 1000 Genomes — Liftover and Integration

## When is this needed?
- FST analysis against reference populations
- Global PCA (placing your cohort in global context)
- ADMIXTURE with continental anchors
- Comparing allele frequencies

## Step 1: Check Genome Build

```bash
# Check chromosome naming in your bim file
head -5 {PREFIX}.bim
# If chr1, chr2... → UCSC style (common in hg38 VCFs)
# If 1, 2...       → Ensembl style (PLINK standard, 1000G)

# 1000G Phase 3 is hg19 (GRCh37), Ensembl-style chromosome names
```

## Step 2: Liftover hg38 → hg19 (if needed)

```bash
# Convert bim to BED format for liftOver
awk '{print "chr"$1, $4-1, $4, $2}' OFS='\t' {PREFIX}.bim \
  > positions_hg38.bed

# Run UCSC liftOver
# Download chain file: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
liftOver positions_hg38.bed hg38ToHg19.over.chain \
  positions_hg19_lifted.bed unlifted.bed

echo "Lifted: $(wc -l < positions_hg19_lifted.bed)"
echo "Failed: $(wc -l < unlifted.bed)"
# Expect: >99% success rate for standard SNP arrays

# Create update files
awk '{print $4, $3}' positions_hg19_lifted.bed > update_positions_hg19.txt
awk '{print $4}' positions_hg19_lifted.bed > snps_lifted.txt

# Apply new coordinates
plink --bfile {PREFIX} \
  --extract snps_lifted.txt \
  --make-bed --out {PREFIX}_hg19_temp

plink --bfile {PREFIX}_hg19_temp \
  --update-map update_positions_hg19.txt \
  --make-bed --out {PREFIX}_hg19
# Note: "Base-pair positions unsorted" warning is normal after liftover
```

## Step 3: Fix Chromosome Names (if needed)

```bash
# If your data has "chr1" but 1000G uses "1"
awk '{$1 = gensub(/^chr/, "", "g", $1); print}' OFS='\t' \
  {PREFIX}_hg19.bim > {PREFIX}_hg19_fix.bim
cp {PREFIX}_hg19.bed {PREFIX}_hg19_fix.bed
cp {PREFIX}_hg19.fam {PREFIX}_hg19_fix.fam
```

## Step 4: Download 1000G Phase 3 VCFs

```bash
# EUR samples for FST (CEU, GBR, FIN, TSI, IBS = 503 samples)
# Full 1000G for global PCA (2,504 samples)

# Sample list files available at:
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# Create population sample lists
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# Extract EUR samples
awk '$3=="EUR" {print $1, $1}' integrated_call_samples_v3.20130502.ALL.panel \
  > european_samples.txt

# Extract by superpopulation
for pop in AFR AMR EAS EUR SAS; do
  awk -v p=$pop '$3==p {print $1, $1}' \
    integrated_call_samples_v3.20130502.ALL.panel > ${pop}_samples.txt
done
```

## Step 5: Per-Chromosome Extraction from 1000G VCFs

```bash
# 1000G Phase 3 VCFs (hg19):
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

for chr in {1..22}; do
  plink --vcf 1000G/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
    --keep {POP}_samples.txt \
    --make-bed \
    --out 1000G_{POP}_chr${chr}
done
```

## After Merging: Create Population Labels File

```bash
# For FST/PCA: need file mapping each sample to population

# Your cohort samples
awk '{print $1, $2, "TARGET"}' {PREFIX}_hg19.fam > populations.txt

# Reference samples
awk '{print $1, $2, $3}' integrated_call_samples_v3.20130502.ALL.panel \
  | grep -f ref_sample_ids.txt >> populations.txt
```

## ALSU Reference: Liftover Results

- Input: 650,181 SNPs (hg38)
- Successfully lifted: 647,854 (99.6%)
- Failed: 2,327 (unmappable between assemblies)
- Typical failure reasons: SNPs in regions that differ between assemblies,
  centromeric/telomeric regions, multi-mapping positions
