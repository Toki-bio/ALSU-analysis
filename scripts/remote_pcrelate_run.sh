set -euo pipefail

BASE="/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/cross_gwas2026_alsu_winter"
WORK="/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/pcrelate_env"
RLIB="$WORK/Rlib"
RUN="$WORK/run"
PLINK2="/staging/conda/envs/bioinfo/bin/plink2"
export LD_LIBRARY_PATH="/staging/miniconda3/lib:/staging/conda/envs/bioinfo/lib:${LD_LIBRARY_PATH:-}"

mkdir -p "$RUN"

echo "PCRELATE_RUN"
echo "BASE	$BASE"
echo "WORK	$WORK"
echo "RLIB	$RLIB"
echo "RUN	$RUN"
echo "LD_LIBRARY_PATH	$LD_LIBRARY_PATH"
date

R_LIBS_USER="$RLIB" Rscript - <<'RS'
required <- c("gdsfmt", "SNPRelate", "GWASTools", "GENESIS", "BiocParallel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
for (pkg in required) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  version <- if (ok) as.character(utils::packageVersion(pkg)) else "NA"
  cat("PKG", pkg, ok, version, "\n", sep = "\t")
}
if (length(missing) > 0) {
  cat("STATUS\tPCRELATE_PACKAGES_MISSING\t", paste(missing, collapse = ","), "\n", sep = "")
  quit(save = "no", status = 20)
}
cat("STATUS\tPCRELATE_PACKAGES_READY\n")
RS

if [ ! -s "$RUN/merged_pruned.bed" ]; then
  "$PLINK2" \
    --bfile "$BASE/merged" \
    --extract "$BASE/cross_prune.prune.in" \
    --make-bed \
    --out "$RUN/merged_pruned" \
    --threads 16
fi

R_LIBS_USER="$RLIB" Rscript - <<'RS'
library(gdsfmt)
library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(BiocParallel)

base <- "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/cross_gwas2026_alsu_winter"
run <- "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/pcrelate_env/run"
bed <- file.path(run, "merged_pruned.bed")
bim <- file.path(run, "merged_pruned.bim")
fam <- file.path(run, "merged_pruned.fam")
gds_file <- file.path(run, "merged_pruned.gds")
pcrelate_pairs <- file.path(run, "pcrelate_kinBtwn.tsv")
pcrelate_summary <- file.path(run, "pcrelate_comparison_summary.txt")

if (!file.exists(gds_file)) {
  snpgdsBED2GDS(bed.fn = bed, bim.fn = bim, fam.fn = fam, out.gdsfn = gds_file)
}

gds <- snpgdsOpen(gds_file)
on.exit({ if (!is.null(gds)) try(snpgdsClose(gds), silent = TRUE) }, add = TRUE)

king <- snpgdsIBDKING(
  gds,
  autosome.only = FALSE,
  remove.monosnp = TRUE,
  missing.rate = 0.01,
  type = "KING-robust",
  useMatrix = TRUE,
  num.thread = 16,
  verbose = TRUE
)
kin_mat <- king$kinship
rownames(kin_mat) <- king$sample.id
colnames(kin_mat) <- king$sample.id
snps_used <- king$snp.id
snpgdsClose(gds)
gds <- NULL

geno_reader <- GdsGenotypeReader(filename = gds_file)
geno_data <- GenotypeData(geno_reader)
on.exit({ try(close(geno_reader), silent = TRUE) }, add = TRUE)

pcair_fit <- pcair(
  geno_data,
  kinobj = kin_mat,
  divobj = kin_mat,
  snp.include = snps_used,
  kin.thresh = 2^(-9/2),
  div.thresh = -2^(-9/2),
  verbose = TRUE
)

pc_count <- min(10, ncol(pcair_fit$vectors))
pcs <- pcair_fit$vectors[, seq_len(pc_count), drop = FALSE]
geno_iter <- GenotypeBlockIterator(geno_data, snpInclude = snps_used)
pcrel <- pcrelate(
  geno_iter,
  pcs = pcs,
  training.set = pcair_fit$unrels,
  ibd.probs = FALSE,
  BPPARAM = BiocParallel::SerialParam(),
  verbose = TRUE
)

kin_btwn <- pcrel$kinBtwn
write.table(kin_btwn, file = pcrelate_pairs, sep = "\t", quote = FALSE, row.names = FALSE)

id1_col <- intersect(c("ID1", "scanID1", "sample.id1", "sample1"), names(kin_btwn))[1]
id2_col <- intersect(c("ID2", "scanID2", "sample.id2", "sample2"), names(kin_btwn))[1]
kin_col <- intersect(c("kin", "kinship"), names(kin_btwn))[1]
if (any(is.na(c(id1_col, id2_col, kin_col)))) {
  stop("Could not infer PC-Relate pair columns from kinBtwn names: ", paste(names(kin_btwn), collapse = ","))
}

sample_map <- read.table(file.path(base, "sample_map.tsv"), header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
cohort_by_iid <- setNames(sample_map$cohort, sample_map$safe_iid)

pair_type <- function(a, b) {
  ca <- cohort_by_iid[[as.character(a)]]
  cb <- cohort_by_iid[[as.character(b)]]
  if (is.null(ca) || is.null(cb)) return("unknown")
  paste(sort(c(ca, cb)), collapse = "+")
}

kin_values <- kin_btwn[[kin_col]]
pair_types <- mapply(pair_type, kin_btwn[[id1_col]], kin_btwn[[id2_col]])
thresholds <- c(0.0442, 0.0884, 0.177, 0.354)

lines <- c(
  "PCRELATE_COMPARISON_SUMMARY",
  paste("KINBTWN_ROWS", nrow(kin_btwn), sep = "\t"),
  paste("PCAIR_UNRELATED", length(pcair_fit$unrels), sep = "\t"),
  paste("PCAIR_RELATED", length(pcair_fit$rels), sep = "\t"),
  paste("PCAIR_PCRELATE_SNPS", length(snps_used), sep = "\t"),
  paste("PCS_USED", pc_count, sep = "\t")
)

for (threshold in thresholds) {
  lines <- c(lines, paste("PCRELATE_GE", threshold, sum(kin_values >= threshold, na.rm = TRUE), sep = "\t"))
}

for (pt in sort(unique(pair_types))) {
  idx <- which(pair_types == pt)
  lines <- c(lines, paste("PAIR_TYPE", pt, length(idx), mean(kin_values[idx], na.rm = TRUE), max(kin_values[idx], na.rm = TRUE), sep = "\t"))
  for (threshold in thresholds) {
    lines <- c(lines, paste("PAIR_TYPE_GE", pt, threshold, sum(kin_values[idx] >= threshold, na.rm = TRUE), sep = "\t"))
  }
}

writeLines(lines, pcrelate_summary)
cat(paste(lines, collapse = "\n"), "\n")
cat("PAIR_OUT\t", pcrelate_pairs, "\n", sep = "")
cat("SUMMARY_OUT\t", pcrelate_summary, "\n", sep = "")
RS