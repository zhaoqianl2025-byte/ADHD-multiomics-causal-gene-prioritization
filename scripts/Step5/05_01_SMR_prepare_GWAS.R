rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
})

# =========================================================
# Step 5.1
# Prepare GWAS summary statistics for SMR analysis
# =========================================================
project_root <- "../.."
infile  <- file.path(project_root, "results", "Step1", "Step1_GWAS_QC_standardized.tsv.gz")
out_dir <- file.path(project_root, "results", "Step5", "01_prepare_gwas")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_smr <- file.path(out_dir, "Step5_ADHD2022_for_SMR.txt")

cat("========================================================\n")
cat("Step 5.1: Prepare GWAS summary statistics for SMR\n")
cat("========================================================\n")
cat("Input file : ", infile, "\n", sep = "")
cat("Output file: ", out_smr, "\n\n", sep = "")

stopifnot(file.exists(infile))

dt <- fread(infile)

required_cols <- c("SNP", "A1", "A2", "beta", "se", "pval", "eaf", "N_total")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

dt_smr <- dt[, .(
  SNP  = as.character(SNP),
  A1   = as.character(A1),
  A2   = as.character(A2),
  freq = as.numeric(eaf),
  b    = as.numeric(beta),
  se   = as.numeric(se),
  p    = as.numeric(pval),
  n    = as.numeric(N_total)
)]

dt_smr[, A1 := toupper(A1)]
dt_smr[, A2 := toupper(A2)]

dt_smr <- dt_smr[
  !is.na(SNP) & SNP != "" &
    !is.na(A1) & A1 %chin% c("A", "C", "G", "T") &
    !is.na(A2) & A2 %chin% c("A", "C", "G", "T") &
    A1 != A2 &
    !is.na(freq) & freq > 0 & freq < 1 &
    !is.na(b) & is.finite(b) &
    !is.na(se) & se > 0 & is.finite(se) &
    !is.na(p) & p > 0 & p <= 1 &
    !is.na(n) & n > 0
]

dt_smr <- unique(dt_smr, by = "SNP")

setorder(dt_smr, p)

fwrite(dt_smr, out_smr, sep = "\t")

summary_dt <- data.table(
  metric = c(
    "n_total_for_smr",
    "min_p",
    "max_p",
    "mean_freq",
    "mean_n"
  ),
  value = c(
    nrow(dt_smr),
    min(dt_smr$p, na.rm = TRUE),
    max(dt_smr$p, na.rm = TRUE),
    mean(dt_smr$freq, na.rm = TRUE),
    mean(dt_smr$n, na.rm = TRUE)
  )
)

summary_file <- file.path(out_dir, "Step5_SMR_prepare_GWAS_summary.tsv")
fwrite(summary_dt, summary_file, sep = "\t")

cat("GWAS rows written for SMR: ", nrow(dt_smr), "\n", sep = "")
cat("Summary file            : ", summary_file, "\n", sep = "")
cat("Done.\n")

message("05_01_SMR_prepare_GWAS.R finished successfully.")