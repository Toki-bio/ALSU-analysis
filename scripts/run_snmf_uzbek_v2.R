#!/usr/bin/env Rscript
# sNMF analysis on Uzbek-only V2 dataset (1,047 samples × 88,722 SNPs)
# Independent validation of K selection for within-Uzbek structure

user_lib <- file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library",
                      paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]))
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
library(LEA)

BED_DIR  <- "/home/copilot/v2/admixture"
BED_FILE <- file.path(BED_DIR, "UZB_v2_admix.bed")
OUT_DIR  <- "/home/copilot/v2/snmf_uzbek"
K_RANGE  <- 2:10
N_REPS   <- 10
ALPHA    <- 10

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUT_DIR)

cat("=== sNMF Uzbek-Only V2 Analysis ===\n")
cat("Input:", BED_FILE, "\n")
cat("K range:", paste(K_RANGE, collapse=", "), "\n")
cat("Replicates per K:", N_REPS, "\n")
cat("Start time:", format(Sys.time()), "\n\n")

# Step 1: Convert BED -> VCF -> geno
bed_base <- sub("\\.bed$", "", BED_FILE)
vcf_file <- file.path(OUT_DIR, "UZB_v2_admix.vcf")

if (!file.exists(vcf_file)) {
    plink_cmd <- sprintf("plink --bfile %s --recode vcf --out %s/UZB_v2_admix --allow-extra-chr",
                         bed_base, OUT_DIR)
    cat("Running:", plink_cmd, "\n")
    ret <- system(plink_cmd)
    if (ret != 0) stop("PLINK conversion failed!")
}

geno_file <- vcf2geno(vcf_file)
cat("Geno file:", geno_file, "\n\n")

# Step 2: Run sNMF
cat("Running sNMF...\n")
t0 <- proc.time()

obj <- snmf(geno_file,
            K = K_RANGE,
            repetitions = N_REPS,
            entropy = TRUE,
            alpha = ALPHA,
            iterations = 200,
            tolerance = 1e-5,
            project = "new",
            CPU = 16)

elapsed <- (proc.time() - t0)[3]
cat(sprintf("sNMF completed in %.1f minutes\n\n", elapsed / 60))

# Step 3: Cross-entropy results
cat("=== Cross-Entropy Results ===\n")
ce_summary <- data.frame(K = integer(), mean_CE = numeric(), sd_CE = numeric(),
                          min_CE = numeric(), best_run = integer())

for (k in K_RANGE) {
    ce <- cross.entropy(obj, K = k)
    ce_summary <- rbind(ce_summary, data.frame(
        K = k, mean_CE = mean(ce), sd_CE = sd(ce),
        min_CE = min(ce), best_run = which.min(ce)
    ))
    cat(sprintf("K=%d:  mean CE = %.6f  (SD = %.6f)  min = %.6f  [best run: %d]\n",
                k, mean(ce), sd(ce), min(ce), which.min(ce)))
}

optimal_k <- ce_summary$K[which.min(ce_summary$mean_CE)]
cat(sprintf("\n>>> Optimal K by cross-entropy: %d (mean CE = %.6f)\n\n",
            optimal_k, min(ce_summary$mean_CE)))

write.csv(ce_summary, file.path(OUT_DIR, "cross_entropy_summary.csv"), row.names = FALSE)

# Step 4: Plot
png(file.path(OUT_DIR, "snmf_cross_entropy_uzbek_v2.png"), width = 800, height = 500, res = 120)
plot(ce_summary$K, ce_summary$mean_CE, type = "b", col = "steelblue", pch = 19, cex = 1.2,
     main = "sNMF Cross-Entropy by K (Uzbek-only V2, n=1047)",
     xlab = "Number of ancestral populations (K)",
     ylab = "Cross-entropy (mean of 10 reps)")
arrows(ce_summary$K, ce_summary$mean_CE - ce_summary$sd_CE,
       ce_summary$K, ce_summary$mean_CE + ce_summary$sd_CE,
       angle = 90, code = 3, length = 0.05, col = "steelblue")
abline(v = optimal_k, col = "red", lty = 2, lwd = 2)
text(optimal_k + 0.3, max(ce_summary$mean_CE), paste0("K=", optimal_k), col = "red", font = 2)
dev.off()
cat("Plot saved\n")

# Step 5: Q-matrices
cat("\nExtracting Q-matrices...\n")
for (k in K_RANGE) {
    best_run <- which.min(cross.entropy(obj, K = k))
    Q_mat <- Q(obj, K = k, run = best_run)
    out_file <- file.path(OUT_DIR, sprintf("snmf_uzbek_v2_K%d.Q", k))
    write.table(Q_mat, out_file, col.names = FALSE, row.names = FALSE, sep = "\t")
    cat(sprintf("  K=%d Q-matrix saved (%d x %d)\n", k, nrow(Q_mat), ncol(Q_mat)))
}

# Step 6: Compare with V2 ADMIXTURE Q-matrices
cat("\n=== Comparison with ADMIXTURE V2 (Uzbek-only) ===\n")
for (k in c(2, 3, 4, 5)) {
    admix_file <- file.path(BED_DIR, sprintf("UZB_v2_admix.%d.Q", k))
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
        cat(sprintf("  K=%d: ADMIXTURE Q-file not found, skipping\n", k))
    }
}

cat("\n=== Done ===\n")
cat("End time:", format(Sys.time()), "\n")
