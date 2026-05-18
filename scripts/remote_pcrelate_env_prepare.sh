set -euo pipefail

WORK="/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/pcrelate_env"
RLIB="$WORK/Rlib"
export LD_LIBRARY_PATH="/staging/miniconda3/lib:/staging/conda/envs/bioinfo/lib:${LD_LIBRARY_PATH:-}"
mkdir -p "$RLIB"

echo "HOST $(hostname)"
echo "DATE $(date -Is)"
echo "WORK $WORK"
echo "RLIB $RLIB"
echo "LD_LIBRARY_PATH $LD_LIBRARY_PATH"
echo "RSCRIPT $(command -v Rscript || true)"
Rscript --version || true

R_LIBS_USER="$RLIB" Rscript - <<'RS'
cat("R_VERSION\t", R.version.string, "\n", sep="")
cat("LIB_PATHS\t", paste(.libPaths(), collapse=";"), "\n", sep="")

needed <- c("GENESIS", "SNPRelate", "GWASTools", "gdsfmt", "igraph", "BiocManager")
for (pkg in needed) {
  available <- requireNamespace(pkg, quietly = TRUE)
  version <- if (available) as.character(utils::packageVersion(pkg)) else "NA"
  cat("PKG\t", pkg, "\t", available, "\t", version, "\n", sep="")
}

missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(setdiff(missing, "BiocManager")) == 0) {
  cat("STATUS\tPCRELATE_READY\n")
  quit(save="no", status=0)
}

cat("STATUS\tPCRELATE_MISSING\t", paste(missing, collapse=","), "\n", sep="")
cat("INSTALL_ATTEMPT\tBEGIN\n")
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = Sys.getenv("R_LIBS_USER"), quiet = FALSE)
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph", lib = Sys.getenv("R_LIBS_USER"), quiet = FALSE)
}
BiocManager::install(c("gdsfmt", "SNPRelate", "GWASTools", "GENESIS"),
                     lib = Sys.getenv("R_LIBS_USER"),
                     ask = FALSE,
                     update = FALSE,
                     quiet = FALSE)
cat("INSTALL_ATTEMPT\tEND\n")

for (pkg in needed) {
  available <- requireNamespace(pkg, quietly = TRUE)
  version <- if (available) as.character(utils::packageVersion(pkg)) else "NA"
  cat("PKG_AFTER\t", pkg, "\t", available, "\t", version, "\n", sep="")
}

if (!requireNamespace("GENESIS", quietly = TRUE) || !requireNamespace("SNPRelate", quietly = TRUE)) {
  quit(save="no", status=20)
}
cat("STATUS\tPCRELATE_READY_AFTER_INSTALL\n")
RS