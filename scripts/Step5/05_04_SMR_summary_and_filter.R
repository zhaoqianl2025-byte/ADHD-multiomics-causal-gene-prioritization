rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
})

# =========================================================
# Step 5.4
# Summarize SMR results and extract SMR + HEIDI passed genes
# =========================================================
project_root <- "../.."
smr_file <- file.path(project_root, "results", "Step5", "03_merge", "Step5_SMR_BrainMeta_allchr.tsv")

out_dir <- file.path(project_root, "results", "Step5", "04_summary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sig_file    <- file.path(out_dir, "Step5_SMR_BrainMeta_significant.tsv")
passed_file <- file.path(out_dir, "Step5_SMR_BrainMeta_passed.tsv")
stats_file  <- file.path(out_dir, "Step5_SMR_BrainMeta_summary_stats.tsv")
gene_list_file <- file.path(out_dir, "Step5_SMR_BrainMeta_passed_genes.txt")

cat("========================================================\n")
cat("Step 5.4: Summarize SMR results and apply HEIDI filter\n")
cat("========================================================\n")
cat("Input file : ", smr_file, "\n", sep = "")
cat("Output dir : ", out_dir, "\n\n", sep = "")

stopifnot(file.exists(smr_file))

smr <- fread(smr_file)

required_cols <- c("Gene", "p_SMR", "p_HEIDI")
missing_cols <- setdiff(required_cols, names(smr))
if (length(missing_cols) > 0) {
  stop("Missing required columns in SMR result file: ", paste(missing_cols, collapse = ", "))
}

smr[, Gene := as.character(Gene)]
smr[, p_SMR := suppressWarnings(as.numeric(p_SMR))]
smr[, p_HEIDI := suppressWarnings(as.numeric(p_HEIDI))]

smr <- smr[!is.na(Gene) & Gene != ""]
smr <- smr[!is.na(p_SMR) & p_SMR > 0 & p_SMR <= 1]

cat("Total tests: ", nrow(smr), "\n", sep = "")

smr_sig <- smr[p_SMR < 5e-6]
cat("SMR significant (p_SMR < 5e-6): ", nrow(smr_sig), "\n", sep = "")

smr_passed <- smr_sig[!is.na(p_HEIDI) & p_HEIDI > 0.01]
cat("SMR + HEIDI passed: ", nrow(smr_passed), "\n", sep = "")

setorder(smr_sig, p_SMR)
setorder(smr_passed, p_SMR)

fwrite(smr_sig, sig_file, sep = "\t")
fwrite(smr_passed, passed_file, sep = "\t")

passed_genes <- unique(smr_passed$Gene)
writeLines(passed_genes, gene_list_file)

summary_dt <- data.table(
  metric = c(
    "n_total_tests",
    "n_smr_significant",
    "n_smr_heidi_passed",
    "n_unique_genes_passed"
  ),
  value = c(
    nrow(smr),
    nrow(smr_sig),
    nrow(smr_passed),
    uniqueN(smr_passed$Gene)
  )
)

fwrite(summary_dt, stats_file, sep = "\t")

cat("\nUnique genes passed:\n")
print(passed_genes)

cat("\nOutput files:\n")
cat(sig_file, "\n")
cat(passed_file, "\n")
cat(stats_file, "\n")
cat(gene_list_file, "\n")

message("05_04_SMR_summary_and_filter.R finished successfully.")