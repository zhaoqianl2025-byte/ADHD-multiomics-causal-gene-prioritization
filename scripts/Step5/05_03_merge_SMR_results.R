rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
})

# =========================================================
# Step 5.3
# Merge chromosome-wise SMR results
# =========================================================
project_root <- "../.."
smr_dir  <- file.path(project_root, "results", "Step5", "02_smr_raw")
out_dir  <- file.path(project_root, "results", "Step5", "03_merge")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "Step5_SMR_BrainMeta_allchr.tsv")
summary_file <- file.path(out_dir, "Step5_SMR_merge_summary.tsv")

cat("========================================================\n")
cat("Step 5.3: Merge chromosome-wise SMR results\n")
cat("========================================================\n")
cat("SMR directory : ", smr_dir, "\n", sep = "")
cat("Output file   : ", out_file, "\n", sep = "")
cat("Summary file  : ", summary_file, "\n\n", sep = "")

stopifnot(dir.exists(smr_dir))

files <- list.files(smr_dir, pattern = "\\.smr$", full.names = TRUE)

if (length(files) == 0) {
  stop("No .smr files were found in: ", smr_dir)
}

cat("SMR files detected: ", length(files), "\n", sep = "")
print(basename(files))

smr_list <- lapply(files, function(f) {
  dt <- fread(f)
  dt[, source_file := basename(f)]
  
  chr_id <- sub("^.*chr([0-9]+)\\.smr$", "\\1", basename(f))
  if (!identical(chr_id, basename(f))) {
    dt[, source_chr := as.integer(chr_id)]
  } else {
    dt[, source_chr := NA_integer_]
  }
  
  dt
})

smr_all <- rbindlist(smr_list, fill = TRUE)

if (nrow(smr_all) == 0) {
  stop("Merged SMR result is empty.")
}

fwrite(smr_all, out_file, sep = "\t")

summary_dt <- data.table(
  metric = c(
    "n_smr_files",
    "n_total_records",
    "n_unique_source_chr"
  ),
  value = c(
    length(files),
    nrow(smr_all),
    uniqueN(smr_all$source_chr)
  )
)

fwrite(summary_dt, summary_file, sep = "\t")

cat("\nTotal SMR records merged: ", nrow(smr_all), "\n", sep = "")
cat("Unique chromosomes      : ", uniqueN(smr_all$source_chr), "\n", sep = "")
cat("Done.\n")

message("05_03_merge_SMR_results.R finished successfully.")