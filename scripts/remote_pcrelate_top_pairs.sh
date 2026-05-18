set -euo pipefail

RUN="/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/pcrelate_env/run"

Rscript - <<'RS'
path <- "/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/pcrelate_env/run/pcrelate_vs_king_pihat_admix.tsv"
x <- read.delim(path, na.strings = c("NA", ""), stringsAsFactors = FALSE)
x <- x[order(-x$pcrelate_kin), ]
cols <- c("iid1", "iid2", "pair_type", "status", "pihat", "king_kinship", "pcrelate_kin", "q5_dist", "q5_dist_bin")
print(head(x[, cols], 20), row.names = FALSE)
RS