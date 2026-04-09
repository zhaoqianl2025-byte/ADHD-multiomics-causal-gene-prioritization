rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggrastr)
})

# =========================================================
# Step 2. LD clumping, LD proxy extraction, SNP panel,
# and Manhattan plot
# =========================================================

# -------------------------------
# 0) Project paths
# -------------------------------
project_root <- ".."
data_dir     <- file.path(project_root, "data", "Step2")
result_dir   <- file.path(project_root, "results", "Step2")
tool_dir     <- file.path(project_root, "tools", "plink")

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

# Input from Step 1
gwas_std_file    <- file.path(project_root, "results", "Step1", "Step1_GWAS_QC_standardized.tsv.gz")
clump_input_file <- file.path(project_root, "results", "Step1", "Step1_GWAS_QC_for_PLINK_clumping.txt")

# PLINK executable
plink_exe <- file.path(tool_dir, "plink.exe")

# LD reference prefix (do not include .bed/.bim/.fam suffix)
bfile_prefix <- file.path(data_dir, "LD_reference", "data_maf0.01_rs_ref")

cat("==== Step 2 Start ====\n")
cat("Standardized GWAS : ", gwas_std_file, "\n")
cat("Clump input       : ", clump_input_file, "\n")
cat("PLINK executable  : ", plink_exe, "\n")
cat("LD reference      : ", bfile_prefix, "\n")
cat("Output directory  : ", result_dir, "\n\n")

# -------------------------------
# 1) Parameters
# -------------------------------
CLUMP_P1 <- 5e-8
CLUMP_P2 <- 1e-2
CLUMP_R2 <- 0.10
CLUMP_KB <- 1000

LD_R2_TH     <- 0.80
LD_WINDOW_KB <- 1000

# -------------------------------
# 2) Basic checks
# -------------------------------
stopifnot(file.exists(gwas_std_file))
stopifnot(file.exists(clump_input_file))
stopifnot(file.exists(plink_exe))
stopifnot(file.exists(paste0(bfile_prefix, ".bed")))
stopifnot(file.exists(paste0(bfile_prefix, ".bim")))
stopifnot(file.exists(paste0(bfile_prefix, ".fam")))

# Normalize key paths for robust external calls
plink_exe_abs <- normalizePath(plink_exe, winslash = "/", mustWork = TRUE)
clump_input_file_abs <- normalizePath(clump_input_file, winslash = "/", mustWork = TRUE)

bfile_dir_abs <- normalizePath(file.path(data_dir, "LD_reference"), winslash = "/", mustWork = TRUE)
bfile_prefix_abs <- file.path(bfile_dir_abs, "data_maf0.01_rs_ref")

clump_out_prefix <- normalizePath(
  file.path(result_dir, "Step2_LD_clumped"),
  winslash = "/",
  mustWork = FALSE
)

# -------------------------------
# 3) Read standardized GWAS
# -------------------------------
gwas <- fread(gwas_std_file)

required_cols <- c("CHR", "BP", "SNP", "A1", "A2", "beta", "se", "pval", "eaf", "N_total", "INFO")
missing_cols <- setdiff(required_cols, names(gwas))
if (length(missing_cols) > 0) {
  stop("Standardized GWAS is missing required columns: ", paste(missing_cols, collapse = ", "))
}

gwas[, CHR := suppressWarnings(as.integer(as.character(CHR)))]
gwas[, BP  := suppressWarnings(as.numeric(as.character(BP)))]

gwas <- gwas[
  !is.na(SNP) & SNP != "" &
    !is.na(CHR) & CHR %in% 1:22 &
    !is.na(BP)  & BP > 0 &
    !is.na(pval) & pval > 0 & pval <= 1
]

setkey(gwas, SNP)

cat("Variants loaded: ", nrow(gwas), "\n")
cat("Genome-wide significant variants: ", gwas[pval < 5e-8, .N], "\n\n")

# -------------------------------
# 4) Run PLINK clumping
# -------------------------------
args_clump <- c(
  "--bfile", bfile_prefix_abs,
  "--clump", clump_input_file_abs,
  "--clump-p1", as.character(CLUMP_P1),
  "--clump-p2", as.character(CLUMP_P2),
  "--clump-r2", as.character(CLUMP_R2),
  "--clump-kb", as.character(CLUMP_KB),
  "--out", clump_out_prefix
)

cat("Running PLINK clumping...\n")
cat(plink_exe_abs, paste(args_clump, collapse = " "), "\n\n")

status_clump <- system2(plink_exe_abs, args = args_clump)

if (!identical(status_clump, 0L) && !identical(status_clump, 0)) {
  stop("PLINK clumping failed with exit status: ", status_clump)
}

clump_file <- paste0(clump_out_prefix, ".clumped")
if (!file.exists(clump_file) || file.size(clump_file) == 0) {
  stop("PLINK clumping output is missing or empty: ", clump_file)
}

clumped <- fread(clump_file, fill = TRUE)
if (!("SNP" %in% names(clumped))) {
  stop("The .clumped file does not contain the SNP column.")
}

lead_snps <- unique(clumped$SNP)
cat("Lead SNPs after clumping: ", length(lead_snps), "\n\n")

# -------------------------------
# 5) Export lead SNP list
# -------------------------------
lead_snps_file <- file.path(result_dir, "Step2_lead_snps.txt")
writeLines(lead_snps, lead_snps_file)

# -------------------------------
# 6) Lead SNP annotation table
# -------------------------------
lead_info <- gwas[J(lead_snps), nomatch = 0]
setorder(lead_info, pval)

lead_info <- lead_info[, .(
  rsid = SNP,
  chr = CHR,
  pos = BP,
  effect_allele = A1,
  other_allele  = A2,
  beta,
  se,
  p = pval,
  eaf,
  N_total,
  INFO
)]

lead_info_file <- file.path(result_dir, "Step2_leadSNP_annotation.tsv")
fwrite(lead_info, lead_info_file, sep = "\t")

# -------------------------------
# 7) Extract LD proxies
# -------------------------------
proxy_list <- vector("list", length(lead_snps))
names(proxy_list) <- lead_snps

for (s in lead_snps) {
  out_prefix <- file.path(result_dir, paste0("tmp_LD_", s))
  out_prefix_abs <- normalizePath(out_prefix, winslash = "/", mustWork = FALSE)
  
  args_ld <- c(
    "--bfile", bfile_prefix_abs,
    "--r2", "gz",
    "--ld-snp", s,
    "--ld-window-kb", as.character(LD_WINDOW_KB),
    "--ld-window", "999999",
    "--ld-window-r2", as.character(LD_R2_TH),
    "--out", out_prefix_abs
  )
  
  status_ld <- suppressWarnings(
    system2(plink_exe_abs, args = args_ld, stdout = TRUE, stderr = TRUE)
  )
  
  ld_file <- paste0(out_prefix_abs, ".ld.gz")
  if (!file.exists(ld_file) || file.size(ld_file) == 0) next
  
  ld <- fread(ld_file)
  
  snpA_col <- intersect(c("SNP_A", "SNP1"), names(ld))
  snpB_col <- intersect(c("SNP_B", "SNP2"), names(ld))
  r2_col   <- intersect(c("R2", "r2"), names(ld))
  
  if (length(snpA_col) == 0 || length(snpB_col) == 0 || length(r2_col) == 0) next
  
  dd <- ld[, .(
    anchor = as.character(get(snpA_col[1])),
    proxy  = as.character(get(snpB_col[1])),
    r2     = as.numeric(get(r2_col[1]))
  )]
  
  dd <- dd[
    anchor == s &
      !is.na(proxy) & proxy != "" &
      !is.na(r2) & r2 >= LD_R2_TH
  ]
  
  if (nrow(dd) > 0) {
    proxy_list[[s]] <- dd
  }
}

proxy_map <- rbindlist(proxy_list, fill = TRUE)
if (nrow(proxy_map) == 0) {
  proxy_map <- data.table(anchor = character(), proxy = character(), r2 = numeric())
}

proxy_map_file <- file.path(result_dir, "Step2_LD_proxy_map_r2ge0.8.tsv")
fwrite(proxy_map, proxy_map_file, sep = "\t")

cat("Proxy pairs: ", nrow(proxy_map), "\n")
cat("Unique proxy SNPs: ", uniqueN(proxy_map$proxy), "\n\n")

# -------------------------------
# 8) Construct lead + proxy SNP panel
# -------------------------------
panel_snps <- unique(c(lead_snps, proxy_map$proxy))

panel <- data.table(SNP = panel_snps)
panel[, SNP_type := fifelse(SNP %in% lead_snps, "LEAD", "PROXY")]

panel <- merge(
  panel,
  proxy_map[, .(SNP = proxy, anchor, r2)],
  by = "SNP",
  all.x = TRUE
)

panel <- merge(
  panel,
  gwas[, .(SNP, CHR, BP, A1, A2, beta, se, pval, eaf, N_total, INFO)],
  by = "SNP",
  all.x = TRUE
)

panel[, type_rank := fifelse(SNP_type == "LEAD", 1L, 2L)]
setorder(panel, type_rank, pval)
panel[, type_rank := NULL]

panel_file <- file.path(result_dir, "Step2_lead_proxy_panel.tsv")
fwrite(panel, panel_file, sep = "\t")

# -------------------------------
# 9) Export rsid lists
# -------------------------------
lead_rsid_file  <- file.path(result_dir, "Step2_lead_snps_rsid.txt")
panel_rsid_file <- file.path(result_dir, "Step2_panel_snps_rsid.txt")

writeLines(unique(lead_snps), lead_rsid_file)
writeLines(unique(panel$SNP), panel_rsid_file)

# -------------------------------
# 10) Summary
# -------------------------------
summary_dt <- data.table(
  item = c(
    "n_total_snps_in_gwas",
    "n_sig_P5e-8_in_gwas",
    "n_lead_snps_after_clump",
    "n_proxy_pairs_r2ge0.8",
    "n_unique_proxy_snps",
    "n_total_panel_snps"
  ),
  value = c(
    nrow(gwas),
    gwas[pval < 5e-8, .N],
    length(lead_snps),
    nrow(proxy_map),
    uniqueN(proxy_map$proxy),
    nrow(panel)
  )
)

summary_file <- file.path(result_dir, "Step2_LD_clumping_summary.tsv")
fwrite(summary_dt, summary_file, sep = "\t")

# =========================================================
# 11) Manhattan plot
# =========================================================
man_dt <- gwas[, .(CHR, BP, SNP, pval)]
man_dt <- man_dt[!is.na(CHR) & CHR %in% 1:22 & !is.na(BP) & BP > 0]
man_dt[, `:=`(
  BP   = as.numeric(BP),
  logP = -log10(pval)
)]
setorder(man_dt, CHR, BP)

chr_len <- man_dt[, .(chr_len = max(BP, na.rm = TRUE)), by = CHR][order(CHR)]

full_chr <- data.table(CHR = 1:22)
chr_len <- merge(full_chr, chr_len, by = "CHR", all.x = TRUE)

med_len <- median(chr_len$chr_len, na.rm = TRUE)
if (!is.finite(med_len) || is.na(med_len)) med_len <- 1e8
chr_len[is.na(chr_len), chr_len := med_len]

chr_len[, chr_len := as.numeric(chr_len)]
chr_len[, chr_start := shift(cumsum(chr_len), fill = 0)]
chr_len[, chr_center := chr_start + chr_len / 2]
chr_len[, `:=`(
  chr_start  = as.numeric(chr_start),
  chr_center = as.numeric(chr_center)
)]

man_dt <- merge(
  man_dt,
  chr_len[, .(CHR, chr_start, chr_center)],
  by = "CHR",
  all.x = TRUE
)

man_dt[, BP_cum := as.numeric(BP) + as.numeric(chr_start)]
man_dt[, chr_parity := factor(CHR %% 2)]

x_max <- max(chr_len$chr_start + chr_len$chr_len, na.rm = TRUE)

lead <- fread(lead_info_file)

if ("rsid" %in% names(lead)) setnames(lead, "rsid", "SNP")
if ("chr"  %in% names(lead)) setnames(lead, "chr",  "CHR")
if ("pos"  %in% names(lead)) setnames(lead, "pos",  "BP")

stopifnot(all(c("SNP", "CHR", "BP") %in% names(lead)))

lead[, CHR := suppressWarnings(as.integer(as.character(CHR)))]
lead[, BP  := suppressWarnings(as.numeric(as.character(BP)))]
lead <- lead[!is.na(SNP) & SNP != ""]
lead <- lead[!is.na(CHR) & CHR %in% 1:22]
lead <- lead[!is.na(BP) & BP > 0]

lead_plot <- merge(
  lead,
  man_dt[, .(SNP, BP_cum, logP, pval)],
  by = "SNP",
  all.x = TRUE
)

need_fallback <- lead_plot[is.na(BP_cum), .(SNP, CHR, BP)]
if (nrow(need_fallback) > 0) {
  setkey(man_dt, CHR, BP)
  setkey(need_fallback, CHR, BP)
  fb <- man_dt[need_fallback, .(SNP = i.SNP, CHR, BP, BP_cum, logP, pval)]
  lead_plot <- rbind(
    lead_plot[!is.na(BP_cum)],
    fb[!is.na(BP_cum)],
    fill = TRUE
  )
}

setkey(lead_plot, SNP)
lead_plot <- unique(lead_plot, by = "SNP")
lead_plot <- lead_plot[!is.na(BP_cum)]

max_labels <- 40
if ("p" %in% names(lead_plot)) {
  lead_plot[, p := suppressWarnings(as.numeric(as.character(p)))]
  setorder(lead_plot, p)
} else {
  setorder(lead_plot, pval)
}
lead_plot <- lead_plot[1:min(max_labels, .N)]

cat("Lead SNPs matched for plotting: ", nrow(lead_plot), "\n")

p_manhattan <- ggplot(man_dt, aes(x = BP_cum, y = logP, color = chr_parity)) +
  ggrastr::geom_point_rast(size = 0.35, raster.dpi = 600) +
  geom_hline(
    yintercept = -log10(5e-8),
    color = "red",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_x_continuous(
    breaks = chr_len$chr_center,
    labels = chr_len$CHR,
    limits = c(0, x_max),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "Manhattan Plot (ADHD2022)",
    x = "Chromosome",
    y = "-log10(P)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title   = element_text(size = 16),
    plot.margin  = margin(8, 24, 8, 10)
  ) +
  coord_cartesian(clip = "off") +
  geom_point(
    data = lead_plot,
    aes(x = BP_cum, y = logP),
    inherit.aes = FALSE,
    color = "black",
    size = 1.15
  ) +
  ggrepel::geom_text_repel(
    data = lead_plot,
    aes(x = BP_cum, y = logP, label = SNP),
    inherit.aes = FALSE,
    size = 3.2,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf
  )

ggsave(
  filename = file.path(result_dir, "Step2_Manhattan_plot.pdf"),
  plot = p_manhattan,
  device = cairo_pdf,
  width = 14,
  height = 5,
  units = "in"
)

ggsave(
  filename = file.path(result_dir, "Step2_Manhattan_plot.png"),
  plot = p_manhattan,
  width = 14,
  height = 5,
  units = "in",
  dpi = 300
)

# -------------------------------
# 12) Finish
# -------------------------------
cat("\n==== Step 2 Done ====\n")
cat("Lead SNP list       : ", lead_snps_file, "\n")
cat("Lead annotation     : ", lead_info_file, "\n")
cat("Proxy map           : ", proxy_map_file, "\n")
cat("Lead+proxy panel    : ", panel_file, "\n")
cat("Lead rsid list      : ", lead_rsid_file, "\n")
cat("Panel rsid list     : ", panel_rsid_file, "\n")
cat("Summary             : ", summary_file, "\n")
cat("Manhattan PDF       : ", file.path(result_dir, "Step2_Manhattan_plot.pdf"), "\n")
cat("Manhattan PNG       : ", file.path(result_dir, "Step2_Manhattan_plot.png"), "\n")

message("Step 2 finished successfully.")