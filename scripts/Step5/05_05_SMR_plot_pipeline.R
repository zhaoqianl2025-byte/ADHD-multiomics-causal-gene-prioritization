rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
})

# =========================================================
# Step 5.5
# SMR candidate gene selection and visualization pipeline
# =========================================================
project_root <- "../.."
infile  <- file.path(project_root, "results", "Step5", "03_merge", "Step5_SMR_BrainMeta_allchr.tsv")
glist_file <- file.path(project_root, "data", "Step5", "glist-hg19.txt")

out_dir  <- file.path(project_root, "results", "Step5", "05_plots")
gene_dir <- file.path(out_dir, "genes")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(gene_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 1) PARAMETERS
# =========================================================
p_smr_cutoff   <- 5e-6
p_heidi_cutoff <- 0.01

cat("========================================================\n")
cat("Step 5.5: SMR candidate gene pipeline\n")
cat("========================================================\n")

stopifnot(file.exists(infile))

# =========================================================
# 2) READ INPUT DATA
# =========================================================
dt <- fread(infile)

required_cols <- c("Gene", "p_SMR", "p_HEIDI", "b_SMR", "ProbeChr", "topSNP")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

# =========================================================
# 3) FILTER SMR RESULTS
# =========================================================
dt[, pass := p_SMR < p_smr_cutoff & p_HEIDI > p_heidi_cutoff]
cand <- dt[pass == TRUE]

if (nrow(cand) == 0) stop("No SMR genes passed filters.")

setorder(cand, p_SMR)

# Keep the most significant record for each gene
cand <- cand[, .SD[1], by = Gene]
setorder(cand, p_SMR)

cat("Candidate genes:", nrow(cand), "\n")

# =========================================================
# 4) EXPORT RESULT TABLES
# =========================================================
fwrite(cand,
       file.path(out_dir, "Step5_candidate_genes.tsv"),
       sep = "\t")

write.xlsx(as.data.frame(cand),
           file = file.path(out_dir, "Step5_candidate_genes.xlsx"),
           overwrite = TRUE)

# Table formatted for manuscript reporting
paper_tab <- cand[, .(
  Gene,
  Lead_SNP = topSNP,
  Chr = ProbeChr,
  beta_SMR = b_SMR,
  p_SMR,
  p_HEIDI
)]

fwrite(paper_tab,
       file.path(out_dir, "Step5_table_for_paper.tsv"),
       sep = "\t")

# =========================================================
# 5) GENERATE GENE-SPECIFIC OUTPUT DIRECTORIES
# =========================================================
for (i in seq_len(nrow(cand))) {
  
  gene <- cand$Gene[i]
  gene_dir_i <- file.path(gene_dir, paste0(sprintf("%02d_", i), gene))
  dir.create(gene_dir_i, showWarnings = FALSE)
  
  fwrite(cand[i],
         file.path(gene_dir_i, paste0(gene, "_summary.tsv")),
         sep = "\t")
}

# =========================================================
# 6) SUMMARY VISUALIZATION
# =========================================================
plot_dt <- copy(cand)
plot_dt[, log10_p := -log10(p_SMR)]
plot_dt[, Gene := factor(Gene, levels = rev(Gene))]

p <- ggplot(plot_dt, aes(log10_p, Gene)) +
  geom_point(aes(size = abs(b_SMR))) +
  theme_bw() +
  labs(
    x = "-log10(p_SMR)",
    y = NULL,
    title = "SMR candidate genes"
  )

ggsave(file.path(out_dir, "Step5_summary_plot.png"), p, width = 8, height = 5)

# =========================================================
# 7) SUMMARY STATISTICS
# =========================================================
stats <- data.table(
  n_total = nrow(dt),
  n_pass  = nrow(cand)
)

fwrite(stats,
       file.path(out_dir, "Step5_stats.tsv"),
       sep = "\t")

cat("Done.\n")