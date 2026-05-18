#!/usr/bin/env bash
set -euo pipefail

WORK="/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/relateadmix_env"
PCRELATE_WORK="/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/pcrelate_env"
PCRELATE_RUN="$PCRELATE_WORK/run"
BASE="/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/cross_gwas2026_alsu_winter"
RUN="$WORK/run"
SRC="$WORK/relateAdmix"
PLINK2="/staging/conda/envs/bioinfo/bin/plink2"
ADMIXTURE="/staging/conda/envs/bioinfo/bin/admixture"
THREADS="${THREADS:-16}"

mkdir -p "$WORK" "$RUN"

echo "HOSTNAME\t$(hostname)"
echo "WORK\t$WORK"
echo "RUN\t$RUN"
echo "THREADS\t$THREADS"
echo "DATE_UTC\t$(date -u +%Y-%m-%dT%H:%M:%SZ)"

echo "SECTION\tpreflight"
command -v git >/dev/null && echo "GIT\t$(command -v git)" || { echo "ERROR\tgit not found"; exit 1; }
command -v make >/dev/null && echo "MAKE\t$(command -v make)" || { echo "ERROR\tmake not found"; exit 1; }
command -v g++ >/dev/null && echo "GXX\t$(command -v g++)" || { echo "ERROR\tg++ not found"; exit 1; }
test -x "$PLINK2" && echo "PLINK2\t$PLINK2" || { echo "ERROR\tplink2 not executable at $PLINK2"; exit 1; }
test -x "$ADMIXTURE" && echo "ADMIXTURE\t$ADMIXTURE" || { echo "ERROR\tadmixture not executable at $ADMIXTURE"; exit 1; }

echo "SECTION\tclone_build_relateadmix"
if [[ ! -d "$SRC/.git" ]]; then
  git clone https://github.com/aalbrechtsen/relateAdmix.git "$SRC"
else
  git -C "$SRC" remote -v | sed 's/^/RELATEADMIX_REMOTE\t/'
  git -C "$SRC" fetch --quiet origin || true
fi
git -C "$SRC" rev-parse HEAD | sed 's/^/RELATEADMIX_COMMIT\t/'

make -C "$SRC/src" -f CPP_Makefile clean >/dev/null 2>&1 || true
make -C "$SRC/src" -f CPP_Makefile
RELATEADMIX_BIN="$SRC/src/relateAdmix"
test -x "$RELATEADMIX_BIN" || { echo "ERROR\trelateAdmix binary was not created"; exit 1; }
echo "RELATEADMIX_BIN\t$RELATEADMIX_BIN"

echo "SECTION\tprepare_pruned_plink"
if [[ -s "$PCRELATE_RUN/merged_pruned.bed" && -s "$PCRELATE_RUN/merged_pruned.bim" && -s "$PCRELATE_RUN/merged_pruned.fam" ]]; then
  cp -f "$PCRELATE_RUN/merged_pruned.bed" "$RUN/merged_pruned.bed"
  cp -f "$PCRELATE_RUN/merged_pruned.bim" "$RUN/merged_pruned.bim"
  cp -f "$PCRELATE_RUN/merged_pruned.fam" "$RUN/merged_pruned.fam"
  echo "PRUNED_SOURCE\t$PCRELATE_RUN/merged_pruned"
else
  "$PLINK2" \
    --bfile "$BASE/merged" \
    --extract "$BASE/cross_prune.prune.in" \
    --make-bed \
    --out "$RUN/merged_pruned" \
    --threads "$THREADS"
  echo "PRUNED_SOURCE\trebuilt_from_cross_prune"
fi
wc -l "$RUN/merged_pruned.fam" | awk '{print "PRUNED_SAMPLES\t"$1}'
wc -l "$RUN/merged_pruned.bim" | awk '{print "PRUNED_SNPS\t"$1}'

echo "SECTION\tadmixture_k5_for_relateadmix"
if [[ ! -s "$RUN/merged_pruned.5.Q" || ! -s "$RUN/merged_pruned.5.P" ]]; then
  (cd "$RUN" && "$ADMIXTURE" -j"$THREADS" merged_pruned.bed 5 | tee admixture_k5_for_relateadmix.log)
else
  echo "ADMIXTURE_K5\treusing_existing_P_Q"
fi
test -s "$RUN/merged_pruned.5.Q" || { echo "ERROR\tmissing $RUN/merged_pruned.5.Q"; exit 1; }
test -s "$RUN/merged_pruned.5.P" || { echo "ERROR\tmissing $RUN/merged_pruned.5.P"; exit 1; }
wc -l "$RUN/merged_pruned.5.Q" | awk '{print "ADMIXTURE_Q_ROWS\t"$1}'
wc -l "$RUN/merged_pruned.5.P" | awk '{print "ADMIXTURE_P_ROWS\t"$1}'

echo "SECTION\trun_relateadmix"
RELATE_PREFIX="$RUN/relateadmix_merged_pruned_k5"
if [[ ! -s "$RELATE_PREFIX.k" && ! -s "$RELATE_PREFIX" && ! -s "$RELATE_PREFIX.out" ]]; then
  "$RELATEADMIX_BIN" \
    -plink "$RUN/merged_pruned" \
    -f "$RUN/merged_pruned.5.P" \
    -q "$RUN/merged_pruned.5.Q" \
    -o "$RELATE_PREFIX" \
    -P "$THREADS" \
    2>&1 | tee "$RUN/relateadmix_run.log"
else
  echo "RELATEADMIX\treusing_existing_output"
fi

RELATE_OUT=""
for candidate in "$RELATE_PREFIX.k" "$RELATE_PREFIX" "$RELATE_PREFIX.out"; do
  if [[ -s "$candidate" ]]; then
    RELATE_OUT="$candidate"
    break
  fi
done
if [[ -z "$RELATE_OUT" ]]; then
  echo "ERROR\tNo RelateAdmix output found for prefix $RELATE_PREFIX"
  ls -lh "$RUN" | sed 's/^/RUN_FILE\t/'
  exit 1
fi
echo "RELATEADMIX_OUT\t$RELATE_OUT"
wc -l "$RELATE_OUT" | awk '{print "RELATEADMIX_LINES\t"$1}'

echo "SECTION\tcompare_relateadmix_to_existing_methods"
cat > "$RUN/compare_relateadmix.R" <<'RSCRIPT'
run_dir <- Sys.getenv("RUN_DIR")
relate_out <- Sys.getenv("RELATE_OUT")
fam_path <- file.path(run_dir, "merged_pruned.fam")
existing_path <- "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/pcrelate_env/run/pcrelate_vs_king_pihat_admix.tsv"
pair_path <- "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/cross_gwas2026_alsu_winter/admixture_king_pihat_pairs.tsv"

make_key <- function(a, b) {
  first <- ifelse(a <= b, a, b)
  second <- ifelse(a <= b, b, a)
  paste(first, second, sep = "__")
}

fam <- read.table(fam_path, stringsAsFactors = FALSE)
rel <- read.table(relate_out, header = TRUE, stringsAsFactors = FALSE)
names(rel) <- tolower(names(rel))
required <- c("ind1", "ind2", "k0", "k1", "k2")
if (!all(required %in% names(rel))) {
  stop("RelateAdmix output columns not recognized: ", paste(names(rel), collapse = ","))
}

rel$iid1 <- fam$V2[rel$ind1 + 1]
rel$iid2 <- fam$V2[rel$ind2 + 1]
if (any(is.na(rel$iid1)) || any(is.na(rel$iid2))) {
  stop("RelateAdmix individual indices did not map cleanly to FAM rows; expected zero-based indices.")
}
rel$key <- make_key(rel$iid1, rel$iid2)
rel$relateadmix_pihat_equiv <- rel$k2 + 0.5 * rel$k1
rel$relateadmix_kin_equiv <- rel$relateadmix_pihat_equiv / 2

pair_screen <- read.table(pair_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
pair_screen$key <- make_key(pair_screen$iid1, pair_screen$iid2)

if (file.exists(existing_path)) {
  existing <- read.table(existing_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
} else {
  existing <- pair_screen
}
existing$key <- make_key(existing$iid1, existing$iid2)

rel_subset <- rel[, c("key", "iid1", "iid2", "k0", "k1", "k2", "relateadmix_pihat_equiv", "relateadmix_kin_equiv")]
merged <- merge(existing, rel_subset, by = "key", all.x = TRUE, suffixes = c("", "_relateadmix"))
if (!"iid1" %in% names(merged) && "iid1.x" %in% names(merged)) names(merged)[names(merged) == "iid1.x"] <- "iid1"
if (!"iid2" %in% names(merged) && "iid2.x" %in% names(merged)) names(merged)[names(merged) == "iid2.x"] <- "iid2"

merged_out <- file.path(run_dir, "relateadmix_vs_king_pihat_pcrelate.tsv")
write.table(merged, merged_out, sep = "\t", quote = FALSE, row.names = FALSE)

thresholds <- c(0.0442, 0.0884, 0.177, 0.354)
summarize <- function(df, label) {
  values <- df$relateadmix_kin_equiv
  pc_values <- if ("pcrelate_kin" %in% names(df)) df$pcrelate_kin else rep(NA_real_, nrow(df))
  c(
    group = label,
    n = nrow(df),
    mean_relateadmix_kin = sprintf("%.9f", mean(values, na.rm = TRUE)),
    median_relateadmix_kin = sprintf("%.9f", median(values, na.rm = TRUE)),
    max_relateadmix_kin = sprintf("%.9f", max(values, na.rm = TRUE)),
    relateadmix_ge_00442 = sum(values >= 0.0442, na.rm = TRUE),
    relateadmix_ge_00884 = sum(values >= 0.0884, na.rm = TRUE),
    relateadmix_ge_0177 = sum(values >= 0.177, na.rm = TRUE),
    relateadmix_ge_0354 = sum(values >= 0.354, na.rm = TRUE),
    pcrelate_ge_00884 = sum(pc_values >= 0.0884, na.rm = TRUE)
  )
}

summary_rows <- list(summarize(merged, "pihat_screen_all"))
if ("status" %in% names(merged)) {
  for (status in sort(unique(merged$status))) {
    summary_rows[[length(summary_rows) + 1]] <- summarize(merged[merged$status == status, , drop = FALSE], paste0("status:", status))
  }
}
if ("pair_type" %in% names(merged)) {
  for (pair_type in sort(unique(merged$pair_type))) {
    summary_rows[[length(summary_rows) + 1]] <- summarize(merged[merged$pair_type == pair_type, , drop = FALSE], paste0("pair_type:", pair_type))
  }
}
if ("q5_dist_bin" %in% names(merged)) {
  for (bin in sort(unique(merged$q5_dist_bin))) {
    summary_rows[[length(summary_rows) + 1]] <- summarize(merged[merged$q5_dist_bin == bin, , drop = FALSE], paste0("q5_dist_bin:", bin))
  }
}

summary <- as.data.frame(do.call(rbind, summary_rows), stringsAsFactors = FALSE)
summary_out <- file.path(run_dir, "relateadmix_vs_king_pihat_pcrelate_summary.tsv")
write.table(summary, summary_out, sep = "\t", quote = FALSE, row.names = FALSE)

all_rel_keys <- rel$key
screen_keys <- pair_screen$key
all_rel_summary <- data.frame(
  metric = c(
    "all_relateadmix_pairs",
    "all_relateadmix_ge_00884",
    "all_relateadmix_ge_00884_in_pihat_screen",
    "all_relateadmix_ge_00884_not_in_pihat_screen",
    "all_relateadmix_ge_0177",
    "all_relateadmix_ge_0354"
  ),
  value = c(
    nrow(rel),
    sum(rel$relateadmix_kin_equiv >= 0.0884, na.rm = TRUE),
    sum(rel$relateadmix_kin_equiv >= 0.0884 & all_rel_keys %in% screen_keys, na.rm = TRUE),
    sum(rel$relateadmix_kin_equiv >= 0.0884 & !(all_rel_keys %in% screen_keys), na.rm = TRUE),
    sum(rel$relateadmix_kin_equiv >= 0.177, na.rm = TRUE),
    sum(rel$relateadmix_kin_equiv >= 0.354, na.rm = TRUE)
  )
)
all_summary_out <- file.path(run_dir, "relateadmix_allpairs_summary.tsv")
write.table(all_rel_summary, all_summary_out, sep = "\t", quote = FALSE, row.names = FALSE)

top <- rel[order(-rel$relateadmix_kin_equiv), c("iid1", "iid2", "k0", "k1", "k2", "relateadmix_pihat_equiv", "relateadmix_kin_equiv")]
top_out <- file.path(run_dir, "relateadmix_top_pairs.tsv")
write.table(head(top, 100), top_out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("MERGED_OUT\t", merged_out, "\n", sep = "")
cat("SUMMARY_OUT\t", summary_out, "\n", sep = "")
cat("ALLPAIRS_SUMMARY_OUT\t", all_summary_out, "\n", sep = "")
cat("TOP_OUT\t", top_out, "\n", sep = "")
print(summary, row.names = FALSE)
print(all_rel_summary, row.names = FALSE)
cat("TOP20\n")
print(head(top, 20), row.names = FALSE)
RSCRIPT

RUN_DIR="$RUN" RELATE_OUT="$RELATE_OUT" /usr/bin/Rscript "$RUN/compare_relateadmix.R" | tee "$RUN/compare_relateadmix.log"

echo "SECTION\tdone"
echo "STATUS\tRELATEADMIX_DONE"