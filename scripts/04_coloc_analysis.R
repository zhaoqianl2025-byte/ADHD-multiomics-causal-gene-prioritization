rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# =========================================================
# Step 4. GTEx v10 brain colocalization analysis and heatmap
#
# Workflow
#   1) Load standardized ADHD GWAS summary statistics
#   2) Load lead SNPs from Step2
#   3) Run Bayesian colocalization against GTEx v10 brain eQTLs
#   4) Export full coloc results and gene-level best results
#   5) Generate heatmap for PP.H4 >= 0.8
# =========================================================

# =========================================================
# 0) Project paths
# =========================================================
project_root <- ".."
gwas_file     <- file.path(project_root, "results", "Step1", "Step1_GWAS_QC_standardized.tsv.gz")
lead_file     <- file.path(project_root, "results", "Step2", "Step2_leadSNP_annotation.tsv")
gtex_brain_dir <- file.path(project_root, "data", "Step3", "GTEx_v10")
out_dir       <- file.path(project_root, "results", "Step4")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========================================================\n")
cat("Step 4: GTEx brain colocalization and heatmap\n")
cat("========================================================\n")
cat("GWAS file  : ", gwas_file, "\n", sep = "")
cat("Lead file  : ", lead_file, "\n", sep = "")
cat("GTEx dir   : ", gtex_brain_dir, "\n", sep = "")
cat("Output dir : ", out_dir, "\n\n", sep = "")

stopifnot(file.exists(gwas_file))
stopifnot(file.exists(lead_file))
stopifnot(dir.exists(gtex_brain_dir))

# =========================================================
# 1) Parameters
# =========================================================
WINDOW <- 5e5
CASE_PROP <- 38691 / (38691 + 186843)
HEATMAP_THR <- 0.8

# =========================================================
# 2) Load ADHD GWAS
# Using Step1 standardized GWAS
# Required columns: SNP, CHR, BP, beta, se, eaf, N_total
# =========================================================
gwas_raw <- fread(gwas_file)

need_cols <- c("SNP", "CHR", "BP", "beta", "se", "eaf", "N_total")
miss <- setdiff(need_cols, names(gwas_raw))
if (length(miss) > 0) {
  stop("Missing required GWAS columns: ", paste(miss, collapse = ", "))
}

gwas <- gwas_raw[, .(
  SNP  = as.character(SNP),
  CHR  = as.integer(CHR),
  BP   = as.integer(BP),
  beta = as.numeric(beta),
  se   = as.numeric(se),
  eaf  = as.numeric(eaf),
  N    = as.numeric(N_total)
)]

gwas <- gwas[
  !is.na(SNP) & SNP != "" &
    !is.na(CHR) & CHR %in% 1:22 &
    !is.na(BP) & BP > 0 &
    !is.na(beta) &
    !is.na(se) & se > 0 &
    !is.na(eaf) & eaf > 0 & eaf < 1 &
    !is.na(N) & N > 0
]

setkey(gwas, CHR, BP)

cat("GWAS SNPs retained: ", nrow(gwas), "\n\n", sep = "")

# =========================================================
# 3) Load lead SNPs from Step2
# Using Step2 lead SNP annotation table
# Required columns: rsid, chr, pos
# =========================================================
leads_raw <- fread(lead_file)

if (!all(c("rsid", "chr", "pos") %in% names(leads_raw))) {
  stop("Lead SNP file must contain rsid, chr, and pos columns.")
}

leads <- leads_raw[, .(
  SNP = as.character(rsid),
  CHR = as.integer(chr),
  BP  = as.integer(pos)
)]

leads <- unique(leads[!is.na(SNP) & SNP != "" & !is.na(CHR) & !is.na(BP)])

if (nrow(leads) == 0) {
  stop("No valid lead SNPs were loaded.")
}

cat("Lead SNPs loaded: ", nrow(leads), "\n\n", sep = "")

# =========================================================
# 4) Locate GTEx brain files
# =========================================================
gtex_files <- list.files(
  gtex_brain_dir,
  pattern = "^Brain_.*\\.v10\\.eQTLs\\.txt\\.gz$|^Brain_.*\\.txt\\.gz$|^Brain_.*\\.gz$",
  full.names = TRUE
)

if (length(gtex_files) == 0) {
  stop("No GTEx brain files were found in data/Step3/GTEx_v10")
}

cat("GTEx brain tissues detected: ", length(gtex_files), "\n", sep = "")
print(head(basename(gtex_files), 10))
cat("\n")

# =========================================================
# 5) Helper functions
# =========================================================
normalize_eqtl_for_coloc <- function(eqtl_dt) {
  nms <- names(eqtl_dt)
  
  gene_col <- intersect(c("external_gene_name", "gene_name", "Gene", "gene", "gene_id", "ensembl_gene_id"), nms)
  snp_col  <- intersect(c("rsids2", "rsid", "SNP", "snp", "variant_id"), nms)
  chr_col  <- intersect(c("chr", "CHR", "chrom"), nms)
  pos_col  <- intersect(c("pos", "position", "BP"), nms)
  beta_col <- intersect(c("beta", "Beta", "slope"), nms)
  se_col   <- intersect(c("se", "SE", "slope_se"), nms)
  eaf_col  <- intersect(c("eaf", "EAF", "af", "AF"), nms)
  n_col    <- intersect(c("n", "N", "ma_samples"), nms)
  p_col    <- intersect(c("p", "pval_nominal", "pval", "P", "p_nominal"), nms)
  
  if (length(gene_col) > 0 && gene_col[1] != "Gene") setnames(eqtl_dt, gene_col[1], "Gene")
  if (length(snp_col)  > 0 && snp_col[1]  != "SNP")  setnames(eqtl_dt, snp_col[1],  "SNP")
  if (length(chr_col)  > 0 && chr_col[1]  != "CHR")  setnames(eqtl_dt, chr_col[1],  "CHR")
  if (length(pos_col)  > 0 && pos_col[1]  != "BP")   setnames(eqtl_dt, pos_col[1],  "BP")
  if (length(beta_col) > 0 && beta_col[1] != "beta") setnames(eqtl_dt, beta_col[1], "beta")
  if (length(se_col)   > 0 && se_col[1]   != "se")   setnames(eqtl_dt, se_col[1],   "se")
  if (length(eaf_col)  > 0 && eaf_col[1]  != "eaf")  setnames(eqtl_dt, eaf_col[1],  "eaf")
  if (length(n_col)    > 0 && n_col[1]    != "n")    setnames(eqtl_dt, n_col[1],    "n")
  if (length(p_col)    > 0 && p_col[1]    != "p")    setnames(eqtl_dt, p_col[1],    "p")
  
  required <- c("Gene", "SNP", "CHR", "BP", "beta", "se", "eaf", "n")
  missing_required <- setdiff(required, names(eqtl_dt))
  if (length(missing_required) > 0) {
    stop("GTEx file missing required columns: ", paste(missing_required, collapse = ", "))
  }
  
  if (!"p" %in% names(eqtl_dt)) {
    eqtl_dt[, p := 2 * pnorm(-abs(as.numeric(beta) / as.numeric(se)))]
  }
  
  eqtl_dt[, Gene := as.character(Gene)]
  eqtl_dt[, SNP  := as.character(SNP)]
  eqtl_dt[, CHR  := suppressWarnings(as.integer(CHR))]
  eqtl_dt[, BP   := suppressWarnings(as.integer(BP))]
  eqtl_dt[, beta := suppressWarnings(as.numeric(beta))]
  eqtl_dt[, se   := suppressWarnings(as.numeric(se))]
  eqtl_dt[, eaf  := suppressWarnings(as.numeric(eaf))]
  eqtl_dt[, n    := suppressWarnings(as.numeric(n))]
  eqtl_dt[, p    := suppressWarnings(as.numeric(p))]
  
  eqtl_dt[
    !is.na(Gene) & Gene != "" &
      !is.na(SNP) & SNP != "" &
      !is.na(CHR) & CHR %in% 1:22 &
      !is.na(BP) & BP > 0 &
      !is.na(beta) &
      !is.na(se) & se > 0 &
      !is.na(eaf) & eaf > 0 & eaf < 1 &
      !is.na(n) & n > 0
  ]
}

run_coloc <- function(gwas_dt, eqtl_dt, case_prop) {
  gwas_dt <- copy(gwas_dt)
  eqtl_dt <- copy(eqtl_dt)
  
  gwas_dt <- gwas_dt[order(se)][!duplicated(SNP)]
  eqtl_dt <- eqtl_dt[order(p, se)][!duplicated(SNP)]
  
  m <- merge(gwas_dt, eqtl_dt, by = "SNP", suffixes = c(".gwas", ".eqtl"))
  
  if (nrow(m) < 50) return(NULL)
  if (anyDuplicated(m$SNP)) return(NULL)
  
  d1 <- list(
    snp = m$SNP,
    beta = m$beta.gwas,
    varbeta = (m$se.gwas)^2,
    MAF = pmin(m$eaf.gwas, 1 - m$eaf.gwas),
    N = m$N.gwas[1],
    type = "cc",
    s = case_prop
  )
  
  d2 <- list(
    snp = m$SNP,
    beta = m$beta.eqtl,
    varbeta = (m$se.eqtl)^2,
    MAF = pmin(m$eaf.eqtl, 1 - m$eaf.eqtl),
    N = m$n[1],
    type = "quant"
  )
  
  coloc.abf(d1, d2)$summary
}

# =========================================================
# 6) Main coloc loop
# =========================================================
res_list <- list()
idx <- 0L

n_locus_skip_small_gwas <- 0L
n_tissue_skip_small_eqtl <- 0L
n_gene_skip_small_eqtl <- 0L
n_coloc_null <- 0L

for (i in seq_len(nrow(leads))) {
  chr  <- leads$CHR[i]
  pos  <- leads$BP[i]
  lead <- leads$SNP[i]
  
  lo <- pos - WINDOW
  hi <- pos + WINDOW
  
  gwas_locus <- gwas[CHR == chr & BP >= lo & BP <= hi]
  
  if (nrow(gwas_locus) < 200) {
    n_locus_skip_small_gwas <- n_locus_skip_small_gwas + 1L
    next
  }
  
  if (i %% 5 == 0 || i == 1) {
    cat("Processing locus ", i, "/", nrow(leads), ": ", lead, "\n", sep = "")
  }
  
  for (f in gtex_files) {
    tissue <- basename(f)
    tissue <- sub("\\.gz$", "", tissue)
    tissue <- sub("\\.txt$", "", tissue)
    tissue <- sub("\\.v10\\.eQTLs$", "", tissue)
    
    eqtl_raw <- fread(f, showProgress = FALSE)
    eqtl <- normalize_eqtl_for_coloc(eqtl_raw)
    eqtl <- eqtl[CHR == chr & BP >= lo & BP <= hi]
    
    if (nrow(eqtl) < 200) {
      n_tissue_skip_small_eqtl <- n_tissue_skip_small_eqtl + 1L
      next
    }
    
    genes <- unique(eqtl$Gene)
    
    for (g in genes) {
      eqtl_g <- eqtl[Gene == g]
      
      if (nrow(eqtl_g) < 50) {
        n_gene_skip_small_eqtl <- n_gene_skip_small_eqtl + 1L
        next
      }
      
      sm <- tryCatch(
        run_coloc(gwas_locus, eqtl_g, case_prop = CASE_PROP),
        error = function(e) NULL
      )
      
      if (is.null(sm)) {
        n_coloc_null <- n_coloc_null + 1L
        next
      }
      
      idx <- idx + 1L
      
      res_list[[idx]] <- data.table(
        Gene = g,
        tissue = tissue,
        locus_leadSNP = lead,
        PP_H0 = unname(sm["PP.H0.abf"]),
        PP_H1 = unname(sm["PP.H1.abf"]),
        PP_H2 = unname(sm["PP.H2.abf"]),
        PP_H3 = unname(sm["PP.H3.abf"]),
        PP_H4 = unname(sm["PP.H4.abf"])
      )
    }
  }
}

res <- rbindlist(res_list, fill = TRUE)

if (nrow(res) == 0) {
  debug_empty <- data.table(
    n_lead_snps = nrow(leads),
    n_gtex_brain_tissues = length(gtex_files),
    n_locus_skip_small_gwas = n_locus_skip_small_gwas,
    n_tissue_skip_small_eqtl = n_tissue_skip_small_eqtl,
    n_gene_skip_small_eqtl = n_gene_skip_small_eqtl,
    n_coloc_null = n_coloc_null,
    n_coloc_records = 0
  )
  fwrite(debug_empty, file.path(out_dir, "Step4_GTExBrain_coloc_debug_stats.tsv"), sep = "\t")
  stop("No valid coloc results were generated. Please inspect Step4_GTExBrain_coloc_debug_stats.tsv")
}

# =========================================================
# 7) Classification and best-gene table
# =========================================================
res[, coloc_tier := fifelse(
  PP_H4 >= 0.8, "STRONG",
  fifelse(PP_H4 >= 0.5, "MID", "WEAK")
)]

setorder(res, -PP_H4)

gene_best <- res[, .SD[1], by = Gene]
setorder(gene_best, -PP_H4)

# =========================================================
# 8) Output coloc tables
# =========================================================
coloc_all_file    <- file.path(out_dir, "Step4_GTExBrain_coloc_all.tsv")
coloc_best_file   <- file.path(out_dir, "Step4_GTExBrain_coloc_gene_best.tsv")
coloc_counts_file <- file.path(out_dir, "Step4_GTExBrain_coloc_counts.tsv")
coloc_debug_file  <- file.path(out_dir, "Step4_GTExBrain_coloc_debug_stats.tsv")

fwrite(res, coloc_all_file, sep = "\t")
fwrite(gene_best, coloc_best_file, sep = "\t")

tier_counts <- res[, .N, by = coloc_tier][order(match(coloc_tier, c("STRONG", "MID", "WEAK")))]
gene_tier_counts <- gene_best[, .N, by = coloc_tier][order(match(coloc_tier, c("STRONG", "MID", "WEAK")))]

counts_dt <- rbind(
  cbind(table = "coloc_all", tier_counts),
  cbind(table = "coloc_gene_best", gene_tier_counts),
  fill = TRUE
)
fwrite(counts_dt, coloc_counts_file, sep = "\t")

debug_dt <- data.table(
  n_lead_snps = nrow(leads),
  n_gtex_brain_tissues = length(gtex_files),
  n_locus_skip_small_gwas = n_locus_skip_small_gwas,
  n_tissue_skip_small_eqtl = n_tissue_skip_small_eqtl,
  n_gene_skip_small_eqtl = n_gene_skip_small_eqtl,
  n_coloc_null = n_coloc_null,
  n_coloc_records = nrow(res),
  n_gene_best = nrow(gene_best),
  n_strong = res[coloc_tier == "STRONG", .N],
  n_mid = res[coloc_tier == "MID", .N],
  n_weak = res[coloc_tier == "WEAK", .N],
  n_strong_gene = gene_best[coloc_tier == "STRONG", .N],
  n_mid_gene = gene_best[coloc_tier == "MID", .N],
  n_weak_gene = gene_best[coloc_tier == "WEAK", .N]
)
fwrite(debug_dt, coloc_debug_file, sep = "\t")

# =========================================================
# 9) Heatmap preparation
# =========================================================
need_cols_heatmap <- c("Gene", "tissue", "locus_leadSNP", "PP_H4", "coloc_tier")
miss_heatmap <- setdiff(need_cols_heatmap, names(res))
if (length(miss_heatmap) > 0) {
  stop("Missing heatmap columns: ", paste(miss_heatmap, collapse = ", "))
}

tissue_map <- c(
  "Brain_Cerebellar_Hemisphere"           = "Cerebellar hemisphere",
  "Brain_Hypothalamus"                    = "Hypothalamus",
  "Brain_Cerebellum"                      = "Cerebellum",
  "Brain_Spinal_cord_cervical_c-1"        = "Spinal cord (cervical C1)",
  "Brain_Frontal_Cortex_BA9"              = "Frontal cortex (BA9)",
  "Brain_Caudate_basal_ganglia"           = "Caudate (basal ganglia)",
  "Brain_Cortex"                          = "Cortex",
  "Brain_Putamen_basal_ganglia"           = "Putamen (basal ganglia)",
  "Brain_Hippocampus"                     = "Hippocampus",
  "Brain_Nucleus_accumbens_basal_ganglia" = "Nucleus accumbens (basal ganglia)",
  "Brain_Anterior_cingulate_cortex_BA24"  = "Anterior cingulate cortex (BA24)",
  "Brain_Amygdala"                        = "Amygdala",
  "Brain_Substantia_nigra"                = "Substantia nigra"
)

df <- as.data.table(res) %>%
  mutate(
    tissue_pretty = ifelse(tissue %in% names(tissue_map), tissue_map[tissue], tissue),
    Gene = as.character(Gene),
    PP_H4 = as.numeric(PP_H4),
    coloc_tier = as.character(coloc_tier)
  )

df_f <- df %>%
  filter(!is.na(PP_H4), PP_H4 >= HEATMAP_THR)

if (nrow(df_f) == 0) {
  stop("No records after filtering PP.H4 >= ", HEATMAP_THR)
}

df_f2 <- df_f %>%
  group_by(Gene, tissue_pretty) %>%
  arrange(desc(PP_H4), .by_group = TRUE) %>%
  summarise(
    PP_H4 = first(PP_H4),
    coloc_tier = first(coloc_tier),
    locus_leadSNP = first(locus_leadSNP),
    .groups = "drop"
  )

heatmap_gene_tissue_file <- file.path(out_dir, "Step4_GTExBrain_coloc_heatmap_input_gene_tissue.tsv")
heatmap_matrix_file      <- file.path(out_dir, "Step4_GTExBrain_coloc_heatmap_matrix.tsv")
fwrite(as.data.table(df_f2), heatmap_gene_tissue_file, sep = "\t")

mat_wide <- df_f2 %>%
  select(Gene, tissue_pretty, PP_H4) %>%
  pivot_wider(
    names_from  = tissue_pretty,
    values_from = PP_H4,
    values_fill = list(PP_H4 = 0)
  )

mat_df <- as.data.frame(mat_wide)
rownames(mat_df) <- mat_df$Gene
mat_df$Gene <- NULL
mat <- as.matrix(mat_df)

fwrite(
  data.table(Gene = rownames(mat), as.data.frame(mat)),
  heatmap_matrix_file,
  sep = "\t"
)

row_tier <- rep("STRONG", nrow(mat))
row_tier <- factor(row_tier, levels = "STRONG")

wanted_tissue_order <- c(
  "Cerebellar hemisphere",
  "Hypothalamus",
  "Cerebellum",
  "Spinal cord (cervical C1)",
  "Frontal cortex (BA9)",
  "Caudate (basal ganglia)",
  "Cortex",
  "Putamen (basal ganglia)",
  "Hippocampus",
  "Nucleus accumbens (basal ganglia)",
  "Anterior cingulate cortex (BA24)",
  "Amygdala",
  "Substantia nigra"
)

keep_cols <- intersect(wanted_tissue_order, colnames(mat))
if (length(keep_cols) > 0) {
  mat <- mat[, keep_cols, drop = FALSE]
}

# =========================================================
# 10) Heatmap drawing
# =========================================================
col_fun <- colorRamp2(
  c(0, HEATMAP_THR, 1),
  c("#FFFFFF", "#F6C08B", "#D55E00")
)

row_anno <- rowAnnotation(
  Tier = row_tier,
  col = list(Tier = c("STRONG" = "#D55E00")),
  annotation_name_side = "top",
  width = unit(6, "mm")
)

ht <- Heatmap(
  mat,
  name = "PP.H4",
  col = col_fun,
  rect_gp = gpar(col = "grey85", lwd = 0.6),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 12, fontface = "plain"),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_title = paste0("GTEx brain colocalization (PP.H4 >= ", HEATMAP_THR, ")"),
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)

png_file <- file.path(out_dir, "Step4_GTExBrain_coloc_heatmap_PPH4_ge_0.8.png")
pdf_file <- file.path(out_dir, "Step4_GTExBrain_coloc_heatmap_PPH4_ge_0.8.pdf")

png(png_file, width = 10, height = 7.5, units = "in", res = 300, type = "cairo")
draw(
  row_anno + ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(4, 10, 4, 6), "mm")
)
dev.off()

pdf(pdf_file, width = 10, height = 7.5, useDingbats = FALSE)
draw(
  row_anno + ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(4, 10, 4, 6), "mm")
)
dev.off()

# =========================================================
# 11) Manifest
# =========================================================
manifest <- data.table(
  file = c(
    "Step4_GTExBrain_coloc_all.tsv",
    "Step4_GTExBrain_coloc_gene_best.tsv",
    "Step4_GTExBrain_coloc_counts.tsv",
    "Step4_GTExBrain_coloc_debug_stats.tsv",
    "Step4_GTExBrain_coloc_heatmap_input_gene_tissue.tsv",
    "Step4_GTExBrain_coloc_heatmap_matrix.tsv",
    "Step4_GTExBrain_coloc_heatmap_PPH4_ge_0.8.png",
    "Step4_GTExBrain_coloc_heatmap_PPH4_ge_0.8.pdf"
  ),
  full_path = c(
    coloc_all_file,
    coloc_best_file,
    coloc_counts_file,
    coloc_debug_file,
    heatmap_gene_tissue_file,
    heatmap_matrix_file,
    png_file,
    pdf_file
  )
)
fwrite(manifest, file.path(out_dir, "Step4_GTExBrain_coloc_manifest.tsv"), sep = "\t")

# =========================================================
# 12) Final report
# =========================================================
cat("\n========================================================\n")
cat("GTEx Brain coloc completed\n")
cat("========================================================\n")
cat("Total coloc records : ", nrow(res), "\n", sep = "")
cat("Gene-best records   : ", nrow(gene_best), "\n", sep = "")
cat("Heatmap genes       : ", nrow(mat), "\n", sep = "")
cat("Heatmap tissues     : ", ncol(mat), "\n\n", sep = "")

cat("coloc_all tier counts:\n")
print(tier_counts)

cat("\ncoloc_gene_best tier counts:\n")
print(gene_tier_counts)

cat("\nDebug stats:\n")
print(debug_dt)

cat("\nSaved files:\n")
cat(coloc_all_file, "\n")
cat(coloc_best_file, "\n")
cat(png_file, "\n")
cat(pdf_file, "\n")

message("04_coloc_analysis.R finished successfully.")