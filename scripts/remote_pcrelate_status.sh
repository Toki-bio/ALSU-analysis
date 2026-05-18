set -euo pipefail

WORK="/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/pcrelate_env"
RLIB="$WORK/Rlib"
export LD_LIBRARY_PATH="/staging/miniconda3/lib:/staging/conda/envs/bioinfo/lib:${LD_LIBRARY_PATH:-}"

echo "PCRELATE_ENV_STATUS"
date
echo "LD_LIBRARY_PATH $LD_LIBRARY_PATH"

echo "ACTIVE_R_BUILD_PROCESSES"
ps -eo pid,etime,pcpu,pmem,cmd \
  | awk '/Rscript|R CMD|gcc|g\+\+|cc1plus/ && !/awk/ {print}' \
  | head -n 40

echo "LOCK_DIRS"
find "$RLIB" -maxdepth 1 -type d -name '00LOCK*' -printf '%f\n' 2>/dev/null | sort || true

echo "PACKAGE_DIRS"
find "$RLIB" -maxdepth 1 -mindepth 1 -type d -printf '%f\n' 2>/dev/null | sort | head -n 160 || true

echo "LIBICU75_CANDIDATES"
find /usr /usr/local /opt /staging -name 'libicui18n.so.75*' -print 2>/dev/null | head -n 40 || true

echo "PACKAGE_FLAGS"
R_LIBS_USER="$RLIB" Rscript - <<'RS'
for (pkg in c("BiocManager", "gdsfmt", "SNPRelate", "GWASTools", "GENESIS", "igraph")) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  where <- if (ok) find.package(pkg) else ""
  cat(pkg, ok, where, "\n")
}
cat("R_LIBS_USER", Sys.getenv("R_LIBS_USER"), "\n")
cat("LD_LIBRARY_PATH", Sys.getenv("LD_LIBRARY_PATH"), "\n")
cat("R_VERSION", as.character(getRversion()), "\n")
RS