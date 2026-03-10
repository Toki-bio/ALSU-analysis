#!/usr/bin/env Rscript
# Completion script: extract Q-matrices, plot, compare with ADMIXTURE, output JSON
# Picks up from where run_snmf.R crashed on plot()

user_lib <- file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library",
                      paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]))
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
library(LEA)

BED_DIR <- "/staging/ALSU-analysis/admixture_analysis/global_admixture"
OUT_DIR <- "/staging/ALSU-analysis/admixture_analysis/snmf_results"
K_RANGE <- 2:10
setwd(OUT_DIR)

# Load existing project
obj <- load.snmfProject("global_for_admixture.snmfProject")
ce_summary <- read.csv("cross_entropy_summary.csv")
optimal_k <- ce_summary$K[which.min(ce_summary$mean_CE)]
cat("Loaded project. Optimal K =", optimal_k, "\n")

# Step 4: Plot (fixed - don't pass xlab/ylab to avoid conflict with method)
png(file.path(OUT_DIR, "snmf_cross_entropy.png"), width = 800, height = 500, res = 120)
plot(ce_summary$K, ce_summary$mean_CE, type = "b", col = "steelblue", pch = 19, cex = 1.2,
     main = "sNMF Cross-Entropy by K",
     xlab = "Number of ancestral populations (K)",
     ylab = "Cross-entropy (mean of 10 reps)")
arrows(ce_summary$K, ce_summary$mean_CE - ce_summary$sd_CE,
       ce_summary$K, ce_summary$mean_CE + ce_summary$sd_CE,
       angle = 90, code = 3, length = 0.05, col = "steelblue")
abline(v = optimal_k, col = "red", lty = 2, lwd = 2)
text(optimal_k + 0.3, max(ce_summary$mean_CE), paste0("K=", optimal_k), col = "red", font = 2)
dev.off()
cat("Cross-entropy plot saved\n")

# Step 5: Extract Q-matrices
cat("\nExtracting Q-matrices...\n")
for (k in K_RANGE) {
    best_run <- which.min(cross.entropy(obj, K = k))
    Q_mat <- Q(obj, K = k, run = best_run)
    out_file <- file.path(OUT_DIR, sprintf("snmf_K%d.Q", k))
    write.table(Q_mat, out_file, col.names = FALSE, row.names = FALSE, sep = "\t")
    cat(sprintf("  K=%d Q-matrix saved (%d x %d)\n", k, nrow(Q_mat), ncol(Q_mat)))
}

# Step 6: Compare with ADMIXTURE Q-matrices
cat("\n=== Comparison with ADMIXTURE ===\n")
admix_dir <- file.path(BED_DIR, "admix_results")
for (k in c(3, 5, 7)) {
    admix_file <- file.path(admix_dir, sprintf("global_for_admixture.%d.Q", k))
    if (!file.exists(admix_file)) {
        admix_file <- file.path(admix_dir, sprintf("K%d.Q", k))
    }
    if (file.exists(admix_file)) {
        admix_Q <- as.matrix(read.table(admix_file))
        best_run <- which.min(cross.entropy(obj, K = k))
        snmf_Q <- Q(obj, K = k, run = best_run)
        cor_mat <- abs(cor(snmf_Q, admix_Q))
        matched <- numeric(k)
        used <- logical(k)
        for (i in 1:k) {
            remaining <- cor_mat[i, !used, drop = FALSE]
            best_j <- which(!used)[which.max(remaining)]
            matched[i] <- cor_mat[i, best_j]
            used[best_j] <- TRUE
        }
        cat(sprintf("  K=%d: mean component correlation = %.4f  (range: %.4f-%.4f)\n",
                    k, mean(matched), min(matched), max(matched)))
    } else {
        cat(sprintf("  K=%d: ADMIXTURE Q-file not found\n", k))
    }
}

# Step 7: JSON output
cat("\nGenerating JSON output...\n")
json_lines <- c("{")
json_lines <- c(json_lines, sprintf('  "optimal_k": %d,', optimal_k))
json_lines <- c(json_lines, sprintf('  "analysis_date": "%s",', format(Sys.time(), "%Y-%m-%d")))
json_lines <- c(json_lines, sprintf('  "n_samples": 2095,'))
json_lines <- c(json_lines, sprintf('  "n_snps": 172537,'))
json_lines <- c(json_lines, '  "cross_entropy": [')
for (i in seq_len(nrow(ce_summary))) {
    comma <- if (i < nrow(ce_summary)) "," else ""
    json_lines <- c(json_lines, sprintf(
        '    {"K": %d, "mean": %.6f, "sd": %.6f, "min": %.6f}%s',
        ce_summary$K[i], ce_summary$mean_CE[i], ce_summary$sd_CE[i], ce_summary$min_CE[i], comma
    ))
}
json_lines <- c(json_lines, "  ]")
json_lines <- c(json_lines, "}")
writeLines(json_lines, file.path(OUT_DIR, "snmf_results.json"))
cat("JSON saved to snmf_results.json\n")

cat("\n=== COMPLETION SCRIPT DONE ===\n")
