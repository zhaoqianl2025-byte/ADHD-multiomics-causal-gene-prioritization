rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
})

# =========================================================
# Step 1. GWAS QC and standardization
# Purpose:
#   1) Read raw ADHD GWAS summary statistics
#   2) Standardize key fields
#   3) Perform basic QC filtering
#   4) Export standardized GWAS table
#   5) Extract genome-wide significant variants
#   6) Prepare PLINK clumping input
# =========================================================

# -------------------------------
# 0) Project paths
# -------------------------------
project_root <- ".."
data_dir     <- file.path(project_root, "data", "Step1")
result_dir   <- file.path(project_root, "results", "Step1")

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

gwas_in <- file.path(data_dir, "ADHD2022_iPSYCH_deCODE_PGC.meta.tsv")
stopifnot(file.exists(gwas_in))

# -------------------------------
# 1) Helper functions
# -------------------------------
pick_first_match <- function(candidates, header) {
  hit <- intersect(candidates, header)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

safe_integer <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_min <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

# -------------------------------
# 2) Read header and detect columns
# -------------------------------
header <- names(fread(gwas_in, nrows = 0))

chr_col <- pick_first_match(c("CHR", "chr", "Chromosome", "chrom", "chromosome"), header)
bp_col  <- pick_first_match(c("BP", "POS", "Position", "base_pair_location", "bp", "pos"), header)
snp_col <- pick_first_match(c("SNP", "RSID", "rsid", "MarkerName", "ID", "variant", "variant_id"), header)

if (is.na(chr_col) || is.na(bp_col) || is.na(snp_col)) {
  stop(
    "Failed to detect required columns.\n",
    "Detected: CHR = ", chr_col, ", BP = ", bp_col, ", SNP = ", snp_col, "\n",
    "Available header: ", paste(header, collapse = ", ")
  )
}

other_candidates <- c(
  "A1", "A2", "OR", "SE", "P", "Nca", "Nco",
  "FRQ_U_186843", "FRQ_A_38691", "INFO"
)
other_exist <- intersect(other_candidates, header)

sel_cols <- unique(c(chr_col, bp_col, snp_col, other_exist))
dt <- fread(gwas_in, select = sel_cols, showProgress = TRUE)

setnames(dt, c(chr_col, bp_col, snp_col), c("CHR", "BP", "SNP"))

# -------------------------------
# 3) Type conversion and standardization
# -------------------------------
dt[, CHR := safe_integer(CHR)]
dt[, BP  := safe_integer(BP)]

if ("OR" %in% names(dt)) dt[, OR := safe_numeric(OR)]
if ("SE" %in% names(dt)) dt[, SE := safe_numeric(SE)]
if ("P"  %in% names(dt)) dt[, P  := safe_numeric(P)]

if ("A1" %in% names(dt)) dt[, A1 := toupper(as.character(A1))]
if ("A2" %in% names(dt)) dt[, A2 := toupper(as.character(A2))]

required_cols <- c("CHR", "BP", "SNP", "OR", "SE", "P", "A1", "A2")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ", paste(missing_cols, collapse = ", "), "\n",
    "Available columns: ", paste(names(dt), collapse = ", ")
  )
}

dt[, beta := log(OR)]
dt[, se   := SE]
dt[, pval := P]

# Prefer control allele frequency
if ("FRQ_U_186843" %in% names(dt)) {
  dt[, eaf := safe_numeric(FRQ_U_186843)]
} else if ("FRQ_A_38691" %in% names(dt)) {
  dt[, eaf := safe_numeric(FRQ_A_38691)]
} else {
  dt[, eaf := NA_real_]
}

# Sample size
if (all(c("Nca", "Nco") %in% names(dt))) {
  dt[, Nca := safe_numeric(Nca)]
  dt[, Nco := safe_numeric(Nco)]
  dt[, N_total := Nca + Nco]
  dt[, s := Nca / (Nca + Nco)]
} else {
  dt[, N_total := NA_real_]
  dt[, s := NA_real_]
}

if ("INFO" %in% names(dt)) {
  dt[, INFO := safe_numeric(INFO)]
} else {
  dt[, INFO := NA_real_]
}

# -------------------------------
# 4) Quality control filtering
# -------------------------------
n_raw <- nrow(dt)

dt <- dt[!is.na(SNP) & SNP != ""]
dt <- dt[!is.na(CHR) & !is.na(BP)]
dt <- dt[CHR %in% 1:22]
dt <- dt[BP > 0]
dt <- dt[!is.na(OR) & OR > 0]
dt <- dt[!is.na(se) & se > 0]
dt <- dt[!is.na(pval) & pval > 0 & pval <= 1]
dt <- dt[!is.na(A1) & !is.na(A2)]
dt <- dt[A1 %chin% c("A", "C", "G", "T")]
dt <- dt[A2 %chin% c("A", "C", "G", "T")]
dt <- dt[A1 != A2]
dt <- dt[is.finite(beta) & is.finite(se) & is.finite(pval)]

# Remove duplicated SNPs, keep row with smallest p-value;
# if tied, keep the first after ordering
setorder(dt, SNP, pval, CHR, BP)
dt <- dt[, .SD[1], by = SNP]

n_qc <- nrow(dt)

# -------------------------------
# 5) Standardized GWAS output
# -------------------------------
gwas_std <- dt[, .(
  CHR, BP, SNP, A1, A2,
  beta, se, pval,
  eaf, N_total, INFO
)]

std_file <- file.path(result_dir, "Step1_GWAS_QC_standardized.tsv.gz")
fwrite(gwas_std, std_file, sep = "\t")

# -------------------------------
# 6) Genome-wide significant variants
# -------------------------------
sig <- gwas_std[pval < 5e-8]
sig_file <- file.path(result_dir, "Step1_GWAS_QC_genomewide_significant.tsv")
fwrite(sig, sig_file, sep = "\t")

# -------------------------------
# 7) PLINK clumping input
# -------------------------------
clump_in <- gwas_std[, .(SNP, P = pval, CHR, BP)]
clump_file <- file.path(result_dir, "Step1_GWAS_QC_for_PLINK_clumping.txt")
fwrite(clump_in, clump_file, sep = "\t")

# -------------------------------
# 8) QC summary
# -------------------------------
qc <- data.table(
  metric = c(
    "n_raw",
    "n_after_qc",
    "n_genomewide_significant",
    "min_p",
    "max_p",
    "mean_eaf",
    "chr_min",
    "chr_max"
  ),
  value = c(
    n_raw,
    n_qc,
    nrow(sig),
    safe_min(gwas_std$pval),
    safe_max(gwas_std$pval),
    safe_mean(gwas_std$eaf),
    safe_min(gwas_std$CHR),
    safe_max(gwas_std$CHR)
  )
)

qc_file <- file.path(result_dir, "Step1_GWAS_QC_summary.tsv")
fwrite(qc, qc_file, sep = "\t")

# -------------------------------
# 9) Console messages
# -------------------------------
cat("==== Step 1 completed ====\n")
cat("Input file:      ", gwas_in, "\n")
cat("Raw variants:    ", n_raw, "\n")
cat("After QC:        ", n_qc, "\n")
cat("Significant SNPs:", nrow(sig), "\n")
cat("Standardized:    ", std_file, "\n")
cat("Significant:     ", sig_file, "\n")
cat("Clumping input:  ", clump_file, "\n")
cat("QC summary:      ", qc_file, "\n")