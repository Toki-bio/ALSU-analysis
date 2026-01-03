# ConvSK pipeline report: raw -> mind20 -> IBD dedup -> VCF for imputation

Date: 2025-12-17

Host/run context: `Biotech2024`, working dir `/staging/ALSU-analysis/winter2025/PLINK_301125_0312`

Dataset prefixes used:

- Raw PLINK: `ConvSK_raw.{bed,bim,fam}`

- After sample missingness filter: `ConvSK_mind20.{bed,bim,fam}`

- After duplicate removal: `ConvSK_mind20_dedup.{bed,bim,fam}`

- After SNP-QC for VCF export (later stage): `ConvSK_mind20_dedup_snpqc.{bed,bim,fam}`


## Step 0. Inputs

Raw dataset summary from PLINK:

- Variants: **654,027**

- Samples: **1,247**

- Sex: all ambiguous in `.fam` (PLINK prints `0 males, 0 females, 1247 ambiguous`)


## Step 1. Sample missingness filter (mind20)

Goal: remove samples with **F_MISS > 0.20**.

Commands:

```bash
plink --bfile ConvSK_raw --missing --out ConvSK_raw_miss
# PLINK .imiss fixed columns: 1=FID 2=IID ... 6=F_MISS
awk 'NR>1 && $6+0 > 0.20{print $1"\t"$2}' ConvSK_raw_miss.imiss > remove_miss20.txt
plink --bfile ConvSK_raw --remove remove_miss20.txt --make-bed --out ConvSK_mind20
```

Observed counts:

- `awk 'NR>1 && $6>0.20{c++} END{print c}' ConvSK_raw_miss.imiss` -> **92** samples above 0.20

- Samples: **1247 -> 1155** after `--remove` (92 removed)


Notes / pitfall fixed:

- Earlier, `remove_miss20.txt` was accidentally written with the *wrong columns* (printing the whole line or wrong field), producing ~1000 malformed lines.

- Using the **fixed .imiss column layout** (FID, IID, ..., F_MISS) avoids header/whitespace parsing issues.


## Step 2. Recompute missingness on mind20

Purpose: get per-sample `F_MISS` values to decide which sample to keep inside each duplicate cluster.

```bash
plink --bfile ConvSK_mind20 --missing --out ConvSK_mind20_miss
```


## Step 3. SNP QC + LD pruning for IBD

Goal: create a pruned set of common, high-quality autosomal SNPs for robust IBD estimation.

Parameters used:

- Autosomes only: `--chr 1-22`

- `--maf 0.1`

- `--geno 0.01`

- `--hwe 1e-06 midp`

- LD pruning: `--indep-pairwise 50 5 0.1`


Command:

```bash
plink --bfile ConvSK_mind20 --chr 1-22 \
  --maf 0.10 --geno 0.01 --hwe 1e-6 midp \
  --indep-pairwise 50 5 0.1 --out ConvSK_mind20_ibd_prune
```

Key outputs from the run:

- After SNP filters: **145,120 variants** retained.

- Pruning removed **100,338** of 145,120 variants; kept **44,782** in `ConvSK_mind20_ibd_prune.prune.in`.


## Step 4. IBD / relatedness computation

Command:

```bash
plink --bfile ConvSK_mind20 --chr 1-22 \
  --extract ConvSK_mind20_ibd_prune.prune.in \
  --genome --out ConvSK_mind20_ibd
```

Observed:

- `--extract`: **44,782** variants used

- Samples: **1,155**

- Pairs with `PI_HAT >= 0.98`: **65**


## Step 5. Deduplication from PI_HAT>=0.98 edges

Method:

- Build a graph where nodes are samples and edges are `PI_HAT >= 0.98`.

- Connected components = duplicate clusters.

- Keep 1 sample per cluster: lowest `F_MISS` from `ConvSK_mind20_miss.imiss`.

- Output:

  - `remove_dups_pihat098.txt` (FID IID) for PLINK `--remove`

  - `dup_clusters_summary.tsv` (cluster_id, kept sample, size, members)


Critical bug fixed:

- Keys must be **FID+IID**, not IID-only. The working version uses a composite key `FID:IID` everywhere.

  (IID-only keys can silently collide and produce empty / wrong FID values.)


Final dedup stats from the successful run:

- Duplicate nodes involved: **106**

- Duplicate clusters: **49**

- Kept (1 per cluster): **49**

- Removed: **57** (`remove_dups_pihat098.txt` lines)


## Step 6. Apply deduplication

```bash
plink --bfile ConvSK_mind20 --remove remove_dups_pihat098.txt --make-bed --out ConvSK_mind20_dedup
```

Observed:

- Samples: **1155 -> 1098** (57 removed)


## Step 7. SNP-QC dataset used for VCF export

A later stage produced `ConvSK_mind20_dedup_snpqc` with:

- Variants: **473,081**

- Samples: **1,098**


When exporting to VCF (`plink --recode vcf bgz`), PLINK warned about allele codes violating the VCF spec.

Root cause: `.bim` contains indel alleles encoded as `I` / `D` (insertion/deletion shorthand), which are invalid VCF alleles.

Observed count:

- `awk '($5=="I"||$5=="D"||$6=="I"||$6=="D")' *.bim` -> **494** variants with I/D alleles.


## Step 8. Export per-chromosome VCFs with only VCF-safe variants

Fix strategy:

- Build `vcf_ok_snps.txt` with variant IDs whose alleles are VCF-safe (A/C/G/T/N, `*`, `<...>`; and allow `.` only for ALT during checks).

- Export each chromosome using `--extract vcf_ok_snps.txt`.


Example:

```bash
plink --bfile ConvSK_mind20_dedup_snpqc --chr 19 --extract vcf_ok_snps.txt --recode vcf bgz --out chr19
# tabix index each chr
for c in {1..22}; do tabix -f -p vcf chr${c}.vcf.gz; done
```


## Step 9. Build imputation input VCF (autosomes)

Commands:
```bash
# concatenate autosomes in order
bcftools concat -Oz -o impute_in.autosomes.vcf.gz chr1.vcf.gz chr2.vcf.gz chr3.vcf.gz chr4.vcf.gz chr5.vcf.gz chr6.vcf.gz chr7.vcf.gz chr8.vcf.gz chr9.vcf.gz chr10.vcf.gz chr11.vcf.gz chr12.vcf.gz chr13.vcf.gz chr14.vcf.gz chr15.vcf.gz chr16.vcf.gz chr17.vcf.gz chr18.vcf.gz chr19.vcf.gz chr20.vcf.gz chr21.vcf.gz chr22.vcf.gz
tabix -f -p vcf impute_in.autosomes.vcf.gz

# biallelic SNPs only
bcftools view -m2 -M2 -v snps -Oz -o impute_in.autosomes.clean.vcf.gz impute_in.autosomes.vcf.gz
tabix -f -p vcf impute_in.autosomes.clean.vcf.gz

# remove duplicate positions
bcftools norm -d all -Oz -o impute_in.autosomes.clean.nodup.vcf.gz impute_in.autosomes.clean.vcf.gz
tabix -f -p vcf impute_in.autosomes.clean.nodup.vcf.gz
```

Observed:
- Duplicate positions removed by `bcftools norm -d all`: **396** in the later run.

## Step 10. Fix non-ASCII sample IDs

Detected non-ASCII sample IDs (2):

- `808_03-25м` (Cyrillic `м`, UTF-8 bytes `d0 bc`)

- `1038_08-176Х-00006` (Cyrillic `Х`, UTF-8 bytes `d0 a5`)


Mapping used:

```text
808_03-25м	808_03-25m
1038_08-176Х-00006	1038_08-176X-00006
```

Apply:

```bash
bcftools reheader -s rename_samples.tsv -o impute_in.autosomes.clean.nodup.rehead.vcf.gz impute_in.autosomes.clean.nodup.vcf.gz
tabix -f -p vcf impute_in.autosomes.clean.nodup.rehead.vcf.gz
```


## Step 11. Validation checks (final file)

File:

- `impute_in.autosomes.clean.nodup.rehead.vcf.gz`


Checks performed:

- Samples: **1098** (`bcftools query -l | wc -l`)

- Variants: **472,191** (`bcftools index -n`)

- ALT='.' count: **0**

- Duplicate CHROM:POS count: **0**

- Non-ASCII sample IDs: **0**

- `bcftools stats` written to `impute_in.autosomes.stats.txt`


## Next steps (recommended)

1. Decide the target reference / build for imputation (hg19/GRCh37 vs hg38/GRCh38) and ensure your CHROM naming matches the server expectation.

2. Confirm REF/ALT are aligned to the chosen reference (strand + ref check). Typical preflight: `bcftools +fixref` / `bcftools norm --check-ref` against the reference FASTA.

3. If the imputation server requires specific variant IDs (or forbids duplicates), consider setting IDs to `CHR:POS:REF:ALT` before upload.

4. Upload `impute_in.autosomes.clean.nodup.rehead.vcf.gz` and its `.tbi` to the imputation service.
