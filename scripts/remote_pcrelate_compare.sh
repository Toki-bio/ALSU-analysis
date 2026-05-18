set -euo pipefail

BASE="/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/cross_gwas2026_alsu_winter"
RUN="/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/pcrelate_env/run"
export LD_LIBRARY_PATH="/staging/miniconda3/lib:/staging/conda/envs/bioinfo/lib:${LD_LIBRARY_PATH:-}"

echo "PCRELATE_COMPARE"
echo "BASE	$BASE"
echo "RUN	$RUN"
date

Rscript - <<'RS'
base <- "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/cross_gwas2026_alsu_winter"
run <- "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/pcrelate_env/run"

pair_file <- file.path(base, "admixture_king_pihat_pairs.tsv")
pcrelate_file <- file.path(run, "pcrelate_kinBtwn.tsv")
merged_file <- file.path(run, "pcrelate_vs_king_pihat_admix.tsv")
summary_file <- file.path(run, "pcrelate_vs_king_pihat_admix_summary.tsv")

pairs <- read.delim(pair_file, stringsAsFactors = FALSE, na.strings = c("NA", ""), check.names = FALSE)
pcrelate <- read.delim(pcrelate_file, stringsAsFactors = FALSE, na.strings = c("NA", ""), check.names = FALSE)

stopifnot(all(c("iid1", "iid2", "status", "pair_type", "pihat") %in% names(pairs)))
stopifnot(all(c("ID1", "ID2", "kin", "nsnp") %in% names(pcrelate)))

pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "__")
pairs$key <- pair_key(pairs$iid1, pairs$iid2)
pcrelate$key <- pair_key(pcrelate$ID1, pcrelate$ID2)
pcrelate_small <- pcrelate[, c("key", "kin", "nsnp")]
names(pcrelate_small) <- c("key", "pcrelate_kin", "pcrelate_nsnp")

merged <- merge(pairs, pcrelate_small, by = "key", all.x = TRUE, sort = FALSE)
merged <- merged[match(pairs$key, merged$key), ]
write.table(merged, file = merged_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

thresholds <- c(0.0442, 0.0884, 0.177, 0.354)

summarize_subset <- function(data, group_field, group_value) {
  kin <- data$pcrelate_kin
  row <- data.frame(
    group_field = group_field,
    group_value = group_value,
    n_pairs = nrow(data),
    n_pcrelate_matched = sum(!is.na(kin)),
    mean_pihat = mean(data$pihat, na.rm = TRUE),
    mean_king_kinship = mean(data$king_kinship, na.rm = TRUE),
    mean_pcrelate_kin = mean(kin, na.rm = TRUE),
    median_pcrelate_kin = median(kin, na.rm = TRUE),
    max_pcrelate_kin = max(kin, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  for (threshold in thresholds) {
    label <- paste0("pcrelate_ge_", gsub("\\.", "", sprintf("%.4f", threshold)))
    row[[label]] <- sum(kin >= threshold, na.rm = TRUE)
    row[[paste0(label, "_rate")]] <- mean(kin >= threshold, na.rm = TRUE)
  }
  row
}

summary_rows <- list(summarize_subset(merged, "all", "all"))
for (field in c("status", "pair_type", "q5_dist_bin")) {
  values <- sort(unique(merged[[field]][!is.na(merged[[field]])]))
  for (value in values) {
    summary_rows[[length(summary_rows) + 1]] <- summarize_subset(merged[merged[[field]] == value & !is.na(merged[[field]]), ], field, value)
  }
}

for (field in c("status_pair_type", "status_q5_dist_bin")) {
  parts <- strsplit(field, "_")[[1]]
  left <- parts[1]
  right <- paste(parts[-1], collapse = "_")
  values <- sort(unique(paste(merged[[left]], merged[[right]], sep = ":")))
  for (value in values) {
    split_value <- strsplit(value, ":", fixed = TRUE)[[1]]
    idx <- merged[[left]] == split_value[1] & merged[[right]] == split_value[2]
    idx[is.na(idx)] <- FALSE
    if (sum(idx) > 0) {
      summary_rows[[length(summary_rows) + 1]] <- summarize_subset(merged[idx, ], field, value)
    }
  }
}

summary <- do.call(rbind, summary_rows)
write.table(summary, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat("MERGED_OUT\t", merged_file, "\n", sep = "")
cat("SUMMARY_OUT\t", summary_file, "\n", sep = "")
cat("SUMMARY_HEAD\n")
print(summary[seq_len(min(nrow(summary), 20)), ], row.names = FALSE)
RS