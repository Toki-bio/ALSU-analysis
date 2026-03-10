#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# sNMF analysis for ALSU pipeline — independent validation of K selection
# Run on server: Rscript run_snmf.R
# Requires: R package 'LEA' (Bioconductor)
#   install: BiocManager::install("LEA")
# Input:  PLINK .bed/.bim/.fam in BED_DIR
# Output: cross-entropy plot, Q-matrices, correlation with ADMIXTURE
# ──────────────────────────────────────────────────────────────────────────────

# Ensure user library is on the path (LEA installed there)
user_lib <- file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library",
                      paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]))
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

library(LEA)

# ── Configuration ─────────────────────────────────────────────────────────────
BED_DIR    <- "/staging/ALSU-analysis/admixture_analysis/global_admixture"
BED_FILE   <- file.path(BED_DIR, "global_for_admixture.bed")
OUT_DIR    <- "/staging/ALSU-analysis/admixture_analysis/snmf_results"
K_RANGE    <- 2:10
N_REPS     <- 10
ALPHA      <- 10     # L2 regularization strength
N_ITER     <- 200    # max iterations
TOLERANCE  <- 1e-5   # convergence tolerance

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUT_DIR)

cat("=== sNMF Analysis ===\n")
cat("Input:", BED_FILE, "\n")
cat("K range:", paste(K_RANGE, collapse=", "), "\n")
cat("Replicates per K:", N_REPS, "\n")
cat("Alpha (regularization):", ALPHA, "\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ── Step 1: Convert BED to geno format (BED → VCF → geno) ─────────────────────
cat("Converting BED to geno format via PLINK → VCF → vcf2geno...\n")

# BED basename (strip .bed extension)
bed_base <- sub("\\.bed$", "", BED_FILE)
vcf_file <- file.path(OUT_DIR, "global_for_admixture.vcf")

# Step 1a: PLINK BED → VCF
if (!file.exists(vcf_file)) {
    plink_cmd <- sprintf("plink --bfile %s --recode vcf --out %s/global_for_admixture --allow-extra-chr",
                         bed_base, OUT_DIR)
    cat("Running:", plink_cmd, "\n")
    ret <- system(plink_cmd)
    if (ret != 0) stop("PLINK conversion failed!")
} else {
    cat("VCF file already exists, skipping PLINK conversion\n")
}

# Step 1b: VCF → .geno  
geno_file <- vcf2geno(vcf_file)
cat("Geno file:", geno_file, "\n\n")

# ── Step 2: Run sNMF ─────────────────────────────────────────────────────────
cat("Running sNMF for K =", min(K_RANGE), "to", max(K_RANGE), "...\n")
t0 <- proc.time()

obj <- snmf(geno_file,
            K = K_RANGE,
            repetitions = N_REPS,
            entropy = TRUE,           # compute cross-entropy
            alpha = ALPHA,            # regularization
            iterations = N_ITER,
            tolerance = TOLERANCE,
            project = "new",
            CPU = 16)

elapsed <- (proc.time() - t0)[3]
cat(sprintf("sNMF completed in %.1f minutes\n\n", elapsed / 60))

# ── Step 3: Cross-entropy results ────────────────────────────────────────────
cat("=== Cross-Entropy Results ===\n")
ce_summary <- data.frame(K = integer(), mean_CE = numeric(), sd_CE = numeric(),
                          min_CE = numeric(), best_run = integer())

for (k in K_RANGE) {
    ce <- cross.entropy(obj, K = k)
    ce_summary <- rbind(ce_summary, data.frame(
        K = k,
        mean_CE = mean(ce),
        sd_CE = sd(ce),
        min_CE = min(ce),
        best_run = which.min(ce)
    ))
    cat(sprintf("K=%d:  mean CE = %.6f  (SD = %.6f)  min = %.6f  [best run: %d]\n",
                k, mean(ce), sd(ce), min(ce), which.min(ce)))
}

# Identify optimal K
optimal_k <- ce_summary$K[which.min(ce_summary$mean_CE)]
cat(sprintf("\n>>> Optimal K by cross-entropy: %d (mean CE = %.6f)\n\n", 
            optimal_k, min(ce_summary$mean_CE)))

# Save summary table
write.csv(ce_summary, file.path(OUT_DIR, "cross_entropy_summary.csv"), row.names = FALSE)

# ── Step 4: Save cross-entropy plot ──────────────────────────────────────────
png(file.path(OUT_DIR, "snmf_cross_entropy.png"), width = 800, height = 500, res = 120)
plot(obj, col = "steelblue", pch = 19, cex = 1.2,
     main = "sNMF Cross-Entropy by K",
     xlab = "Number of ancestral populations (K)",
     ylab = "Cross-entropy")
abline(v = optimal_k, col = "red", lty = 2, lwd = 2)
text(optimal_k + 0.3, max(ce_summary$mean_CE), paste0("K=", optimal_k), col = "red", font = 2)
dev.off()
cat("Cross-entropy plot saved to snmf_cross_entropy.png\n")

# ── Step 5: Extract Q-matrices for best runs ─────────────────────────────────
cat("\nExtracting Q-matrices...\n")
for (k in K_RANGE) {
    best_run <- which.min(cross.entropy(obj, K = k))
    Q_mat <- Q(obj, K = k, run = best_run)
    out_file <- file.path(OUT_DIR, sprintf("snmf_K%d.Q", k))
    write.table(Q_mat, out_file, col.names = FALSE, row.names = FALSE, sep = "\t")
    cat(sprintf("  K=%d Q-matrix saved (%d x %d)\n", k, nrow(Q_mat), ncol(Q_mat)))
}

# ── Step 6: Compare with ADMIXTURE Q-matrices (if available) ─────────────────
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
        
        # Compute max abs correlation for best column alignment
        cor_mat <- abs(cor(snmf_Q, admix_Q))
        # Greedy matching
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
        cat(sprintf("  K=%d: ADMIXTURE Q-file not found at %s\n", k, admix_file))
    }
}

# ── Step 7: Save full results as JSON for webpage ────────────────────────────
cat("\nGenerating JSON output...\n")
json_lines <- c("{")
json_lines <- c(json_lines, sprintf('  "optimal_k": %d,', optimal_k))
json_lines <- c(json_lines, sprintf('  "analysis_date": "%s",', format(Sys.time(), "%Y-%m-%d")))
json_lines <- c(json_lines, sprintf('  "alpha": %d,', ALPHA))
json_lines <- c(json_lines, sprintf('  "n_reps": %d,', N_REPS))
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
cat("JSON results saved to snmf_results.json\n")

cat(sprintf("\n=== DONE === (Total time: %.1f min)\n", (proc.time() - t0)[3] / 60))
cat(sprintf("Optimal K: %d\n", optimal_k))
cat("Output dir:", OUT_DIR, "\n")
