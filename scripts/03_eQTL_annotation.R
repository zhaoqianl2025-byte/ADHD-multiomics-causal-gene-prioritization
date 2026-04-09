rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
})

# =========================================================
# Step 3: Brain eQTL annotation
# =========================================================

# =========================================================
# 0) Project paths
# =========================================================
project_root <- ".."
gwas_file    <- file.path(project_root, "results", "Step1", "Step1_GWAS_QC_standardized.tsv.gz")
panel_file   <- file.path(project_root, "results", "Step2", "Step2_lead_proxy_panel.tsv")
gtex_root    <- file.path(project_root, "data", "Step3", "GTEx_v10")
out_root     <- file.path(project_root, "results", "Step3")

dir_match  <- file.path(out_root, "matches")
dir_panels <- file.path(out_root, "panels")
dir_tables <- file.path(dir_panels, "tables")
dir_debug  <- file.path(dir_panels, "debug")
dir_plots  <- file.path(dir_panels, "plots")

dir.create(out_root,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_match,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_panels, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_debug,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots,  showWarnings = FALSE, recursive = TRUE)

cat("Running Step 3: Brain eQTL annotation...
")

stopifnot(file.exists(gwas_file))
stopifnot(file.exists(panel_file))
stopifnot(dir.exists(gtex_root))

# =========================================================
# 1) Parameters
# =========================================================
TH_STRONG   <- 5e-8
TH_MODERATE <- 1e-5
TH_WEAK     <- 1e-3

# =========================================================
# 2) Functions
# =========================================================
classify_gtex_evidence <- function(p) {
  fifelse(
    is.na(p), "NoEvidence",
    fifelse(
      p < TH_STRONG, "Strong",
      fifelse(
        p < TH_MODERATE, "Moderate",
        fifelse(p < TH_WEAK, "Weak", "NoEvidence")
      )
    )
  )
}

clean_gtex_tissue_label <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.v10\\.eQTLs$", "", x)
  x <- gsub("\\.txt$", "", x)
  x <- gsub("^Brain_", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("Anterior cingulate cortex BA24", "ACC (BA24)", x, fixed = TRUE)
  x <- gsub("Frontal Cortex BA9", "Frontal Cortex (BA9)", x, fixed = TRUE)
  x <- gsub("Caudate basal ganglia", "Caudate", x, fixed = TRUE)
  x <- gsub("Putamen basal ganglia", "Putamen", x, fixed = TRUE)
  x <- gsub("Nucleus accumbens basal ganglia", "Nucleus accumbens", x, fixed = TRUE)
  x <- gsub("Spinal cord cervical c-1", "Spinal cord (C1)", x, fixed = TRUE)
  x <- gsub("Substantia nigra", "Substantia nigra", x, fixed = TRUE)
  x
}

normalize_gtex_columns <- function(dt) {
  nms <- names(dt)

  snp_col  <- intersect(c("rsids2", "rsid", "SNP", "snp", "variant_id"), nms)
  gene_col <- intersect(c("external_gene_name", "gene_name", "Gene", "gene", "gene_id", "ensembl_gene_id"), nms)
  p_col    <- intersect(c("p", "pval_nominal", "pval", "P", "p_nominal"), nms)
  chr_col  <- intersect(c("chr", "CHR", "chrom"), nms)
  pos_col  <- intersect(c("pos", "position", "BP"), nms)
  a1_col   <- intersect(c("a1", "A1", "alt", "effect_allele"), nms)
  a2_col   <- intersect(c("a2", "A2", "ref", "other_allele"), nms)
  beta_col <- intersect(c("beta", "Beta", "slope"), nms)
  se_col   <- intersect(c("se", "SE", "slope_se"), nms)
  n_col    <- intersect(c("n", "N", "ma_samples"), nms)
  af_col   <- intersect(c("eaf", "EAF", "af", "AF"), nms)

  if (length(snp_col)  > 0 && snp_col[1]  != "SNP")  setnames(dt, snp_col[1],  "SNP")
  if (length(gene_col) > 0 && gene_col[1] != "Gene") setnames(dt, gene_col[1], "Gene")
  if (length(p_col)    > 0 && p_col[1]    != "P")    setnames(dt, p_col[1],    "P")
  if (length(chr_col)  > 0 && chr_col[1]  != "CHR")  setnames(dt, chr_col[1],  "CHR")
  if (length(pos_col)  > 0 && pos_col[1]  != "BP")   setnames(dt, pos_col[1],  "BP")
  if (length(a1_col)   > 0 && a1_col[1]   != "A1")   setnames(dt, a1_col[1],   "A1")
  if (length(a2_col)   > 0 && a2_col[1]   != "A2")   setnames(dt, a2_col[1],   "A2")
  if (length(beta_col) > 0 && beta_col[1] != "Beta") setnames(dt, beta_col[1], "Beta")
  if (length(se_col)   > 0 && se_col[1]   != "SE")   setnames(dt, se_col[1],   "SE")
  if (length(n_col)    > 0 && n_col[1]    != "N")    setnames(dt, n_col[1],    "N")
  if (length(af_col)   > 0 && af_col[1]   != "EAF")  setnames(dt, af_col[1],   "EAF")

  if (!"SNP"  %in% names(dt)) dt[, SNP  := NA_character_]
  if (!"Gene" %in% names(dt)) dt[, Gene := NA_character_]
  if (!"P"    %in% names(dt)) dt[, P    := NA_real_]
  if (!"CHR"  %in% names(dt)) dt[, CHR  := NA_integer_]
  if (!"BP"   %in% names(dt)) dt[, BP   := NA_integer_]
  if (!"A1"   %in% names(dt)) dt[, A1   := NA_character_]
  if (!"A2"   %in% names(dt)) dt[, A2   := NA_character_]
  if (!"EAF"  %in% names(dt)) dt[, EAF  := NA_real_]
  if (!"Beta" %in% names(dt)) dt[, Beta := NA_real_]
  if (!"SE"   %in% names(dt)) dt[, SE   := NA_real_]
  if (!"N"    %in% names(dt)) dt[, N    := NA_real_]

  dt[, SNP  := as.character(SNP)]
  dt[, Gene := as.character(Gene)]
  dt[is.na(Gene), Gene := ""]
  dt[, P    := suppressWarnings(as.numeric(P))]
  dt[, CHR  := suppressWarnings(as.integer(CHR))]
  dt[, BP   := suppressWarnings(as.integer(BP))]
  dt[, A1   := as.character(A1)]
  dt[, A2   := as.character(A2)]
  dt[, EAF  := suppressWarnings(as.numeric(EAF))]
  dt[, Beta := suppressWarnings(as.numeric(Beta))]
  dt[, SE   := suppressWarnings(as.numeric(SE))]
  dt[, N    := suppressWarnings(as.numeric(N))]

  idx_variant <- grepl("^(chr)?[0-9XYM]+_[0-9]+_", dt$SNP)
  if (any(idx_variant, na.rm = TRUE)) {
    tmp <- tstrsplit(dt$SNP[idx_variant], "_", fixed = TRUE)
    if (length(tmp) >= 2) {
      chr_parsed <- gsub("^chr", "", tmp[[1]], ignore.case = TRUE)
      bp_parsed  <- suppressWarnings(as.integer(tmp[[2]]))
      dt[idx_variant & is.na(CHR), CHR := suppressWarnings(as.integer(chr_parsed))]
      dt[idx_variant & is.na(BP),  BP  := bp_parsed]
    }
  }

  dt
}

read_gtex_hits_one_file <- function(f, panel_snps, panel_chr_bp) {
  message("Reading: ", basename(f))

  dt <- tryCatch(
    fread(f, showProgress = FALSE),
    error = function(e) {
      con <- gzfile(f, open = "rt")
      on.exit(close(con), add = TRUE)
      txt <- readLines(con, warn = FALSE)
      fread(text = paste(txt, collapse = "\n"), showProgress = FALSE)
    }
  )

  dt <- normalize_gtex_columns(dt)
  dt <- dt[!is.na(P)]
  if (nrow(dt) == 0) return(NULL)

  dt_rsid <- dt[!is.na(SNP) & SNP != "" & SNP %chin% panel_snps]
  dt_pos <- merge(dt[!is.na(CHR) & !is.na(BP)], panel_chr_bp, by = c("CHR", "BP"))

  dt_use <- unique(
    rbindlist(list(dt_rsid, dt_pos), fill = TRUE),
    by = c("SNP", "Gene", "CHR", "BP", "P")
  )
  if (nrow(dt_use) == 0) return(NULL)

  tissue_name <- basename(f)
  tissue_name <- sub("\\.gz$", "", tissue_name)
  tissue_name <- sub("\\.txt$", "", tissue_name)
  tissue_name <- sub("\\.tsv$", "", tissue_name)
  tissue_name <- sub("\\.signif.*$", "", tissue_name)

  dt_use[, Tissue := tissue_name]
  dt_use[, Tissue_label := clean_gtex_tissue_label(Tissue)]
  dt_use[, Gene_nonempty := Gene != ""]
  dt_use[, Gene_priority := fifelse(Gene_nonempty, 1L, 0L)]
  dt_use[]
}

# =========================================================
# 2) Load Step2 SNP panel
# =========================================================
panel <- fread(panel_file)
stopifnot("SNP" %in% names(panel))

snp_type_col <- intersect(c("SNP_type", "snp_type", "TYPE", "type"), names(panel))
if ("SNP_type" %in% names(panel)) {
  panel[, SNP_type := as.character(SNP_type)]
} else if (length(snp_type_col) > 0) {
  setnames(panel, snp_type_col[1], "SNP_type")
  panel[, SNP_type := as.character(SNP_type)]
} else {
  panel[, SNP_type := "UNKNOWN"]
}
panel[is.na(SNP_type) | SNP_type == "", SNP_type := "UNKNOWN"]
panel <- unique(panel[!is.na(SNP) & SNP != ""])

panel_map <- unique(panel[, .(
  SNP,
  SNP_type,
  CHR = if ("CHR" %in% names(panel)) as.integer(CHR) else NA_integer_,
  BP  = if ("BP"  %in% names(panel)) as.integer(BP) else NA_integer_,
  A1  = if ("A1"  %in% names(panel)) as.character(A1) else NA_character_,
  A2  = if ("A2"  %in% names(panel)) as.character(A2) else NA_character_,
  EAF = if ("eaf" %in% names(panel)) as.numeric(eaf) else if ("EAF" %in% names(panel)) as.numeric(EAF) else NA_real_,
  beta = if ("beta" %in% names(panel)) as.numeric(beta) else NA_real_,
  se   = if ("se"   %in% names(panel)) as.numeric(se) else NA_real_,
  pval = if ("pval" %in% names(panel)) as.numeric(pval) else NA_real_,
  INFO = if ("INFO" %in% names(panel)) as.numeric(INFO) else NA_real_,
  anchor = if ("anchor" %in% names(panel)) as.character(anchor) else NA_character_,
  r2 = if ("r2" %in% names(panel)) as.numeric(r2) else NA_real_,
  N_total = if ("N_total" %in% names(panel)) as.numeric(N_total) else NA_real_
)])

panel_snps   <- unique(panel_map$SNP)
panel_chr_bp <- unique(panel_map[!is.na(CHR) & !is.na(BP), .(CHR, BP)])

# =========================================================
# 3) Load Step1 GWAS map
# =========================================================
gwas_map <- fread(gwas_file, select = c("SNP", "CHR", "BP", "A1", "A2"))
gwas_map[, SNP := as.character(SNP)]
gwas_map[, CHR_gwas := as.integer(CHR)]
gwas_map[, BP_gwas  := as.integer(BP)]
gwas_map[, A1_gwas  := as.character(A1)]
gwas_map[, A2_gwas  := as.character(A2)]
gwas_map[, c("CHR", "BP", "A1", "A2") := NULL]
gwas_map <- unique(gwas_map[!is.na(SNP) & SNP != ""])

# =========================================================
# 4) Locate GTEx brain files
# =========================================================
all_gtex_files <- list.files(
  gtex_root,
  pattern = "\\.txt\\.gz$|\\.gz$",
  full.names = TRUE,
  recursive = FALSE
)

gtex_files <- all_gtex_files[
  grepl("^Brain_.*\\.txt\\.gz$|^Brain_.*\\.gz$", basename(all_gtex_files))
]
if (length(gtex_files) == 0) stop("No brain GTEx .gz files were found. Please check the GTEx directory.")

cat("Total GTEx .gz files found : ", length(all_gtex_files), "\n", sep = "")
cat("Brain GTEx .gz files used  : ", length(gtex_files), "\n", sep = "")
print(head(basename(gtex_files), 10))
cat("\n")

# =========================================================
# 5) Match panel SNPs to GTEx brain eQTL files
# =========================================================
all_hits_list <- lapply(
  gtex_files,
  read_gtex_hits_one_file,
  panel_snps = panel_snps,
  panel_chr_bp = panel_chr_bp
)
all_hits_list <- Filter(Negate(is.null), all_hits_list)
if (length(all_hits_list) == 0) stop("No GTEx brain hits were found after SNP and CHR+BP matching.")

hits <- rbindlist(all_hits_list, fill = TRUE)

panel_map_snp <- unique(panel_map[, .(
  SNP,
  SNP_type_panel = SNP_type,
  CHR_panel = CHR,
  BP_panel  = BP,
  A1_panel  = A1,
  A2_panel  = A2,
  EAF_panel = EAF
)])

panel_map_pos <- unique(panel_map[, .(
  CHR,
  BP,
  SNP_pos = SNP,
  SNP_type_pos = SNP_type,
  A1_pos = A1,
  A2_pos = A2,
  EAF_pos = EAF
)])

hits <- merge(hits, panel_map_snp, by = "SNP", all.x = TRUE)
hits <- merge(hits, panel_map_pos, by = c("CHR", "BP"), all.x = TRUE)

if (!"SNP_type" %in% names(hits)) hits[, SNP_type := NA_character_]
hits[(is.na(SNP_type) | SNP_type == "") & !is.na(SNP_type_panel) & SNP_type_panel != "", SNP_type := SNP_type_panel]
hits[(is.na(SNP_type) | SNP_type == "") & !is.na(SNP_type_pos)   & SNP_type_pos   != "", SNP_type := SNP_type_pos]
hits[is.na(SNP_type) | SNP_type == "", SNP_type := "UNKNOWN"]

hits[(is.na(A1) | A1 == "") & !is.na(A1_panel) & A1_panel != "", A1 := A1_panel]
hits[(is.na(A2) | A2 == "") & !is.na(A2_panel) & A2_panel != "", A2 := A2_panel]
hits[is.na(EAF) & !is.na(EAF_panel), EAF := EAF_panel]
hits[(is.na(A1) | A1 == "") & !is.na(A1_pos) & A1_pos != "", A1 := A1_pos]
hits[(is.na(A2) | A2 == "") & !is.na(A2_pos) & A2_pos != "", A2 := A2_pos]
hits[is.na(EAF) & !is.na(EAF_pos), EAF := EAF_pos]
hits[is.na(CHR) & !is.na(CHR_panel), CHR := CHR_panel]
hits[is.na(BP)  & !is.na(BP_panel),  BP  := BP_panel]
hits[, logP := -log10(P)]

ord_cols <- c("SNP", "SNP_type", "Tissue", "Tissue_label", "Gene", "P",
              "CHR", "BP", "A1", "A2", "EAF", "Beta", "SE", "N", "logP")
ord_cols <- intersect(ord_cols, names(hits))
setcolorder(hits, ord_cols)

# =========================================================
# 6) Save match results
# =========================================================
all_hits_file <- file.path(dir_match, "Step3_GTExBrain_hits_ALL.tsv.gz")
fwrite(hits, all_hits_file, sep = "\t")

setorder(hits, SNP, Tissue, -Gene_priority, P, Gene)
top1_hits <- hits[, .SD[1], by = .(SNP, Tissue)]
if (!"SNP_type" %in% names(top1_hits)) top1_hits[, SNP_type := "UNKNOWN"]
top1_hits[is.na(SNP_type) | SNP_type == "", SNP_type := "UNKNOWN"]

top1_hits_file <- file.path(dir_match, "Step3_GTExBrain_hits_TOP1_perSNP_Tissue.tsv.gz")
fwrite(top1_hits, top1_hits_file, sep = "\t")

summary_by_tissue <- unique(top1_hits[, .(SNP, Tissue, Tissue_label, SNP_type)])[
  , .(n_unique_SNP = uniqueN(SNP)), by = .(Tissue, Tissue_label, SNP_type)
]
summary_by_tissue_file <- file.path(dir_match, "Step3_GTExBrain_hits_summary_byTissue.tsv")
fwrite(summary_by_tissue, summary_by_tissue_file, sep = "\t")

summary_by_gene <- hits[Gene != "",
                        .(
                          n_hits = .N,
                          n_unique_SNP = uniqueN(SNP),
                          top_p = min(P, na.rm = TRUE),
                          top_tissue = Tissue[which.min(P)],
                          top_tissue_label = Tissue_label[which.min(P)],
                          top_snp = SNP[which.min(P)]
                        ),
                        by = Gene
][order(top_p)]
summary_by_gene_file <- file.path(dir_match, "Step3_GTExBrain_hits_summary_byGene.tsv")
fwrite(summary_by_gene, summary_by_gene_file, sep = "\t")

debug_missing_gene_all  <- hits[Gene == ""]
debug_missing_gene_top1 <- top1_hits[Gene == ""]
debug_missing_snptype   <- hits[is.na(SNP_type) | SNP_type == ""]

fwrite(debug_missing_gene_all,  file.path(dir_match, "Step3_GTExBrain_hits_ALL_missingGene.tsv.gz"), sep = "\t")
fwrite(debug_missing_gene_top1, file.path(dir_match, "Step3_GTExBrain_hits_TOP1_missingGene.tsv.gz"), sep = "\t")
fwrite(debug_missing_snptype,   file.path(dir_match, "Step3_GTExBrain_hits_ALL_missingSNPtype.tsv.gz"), sep = "\t")

debug_stats_match <- data.table(
  n_panel_snps = uniqueN(panel$SNP),
  n_all_hits = nrow(hits),
  n_top1_rows = nrow(top1_hits),
  n_unique_hit_snps = uniqueN(hits$SNP),
  n_unique_hit_genes = uniqueN(hits[Gene != "", Gene]),
  n_all_hits_missing_gene = nrow(debug_missing_gene_all),
  n_top1_missing_gene = nrow(debug_missing_gene_top1),
  n_all_hits_missing_snptype = nrow(debug_missing_snptype)
)
fwrite(debug_stats_match, file.path(dir_match, "Step3_GTExBrain_hits_debug_stats.tsv"), sep = "\t")

# =========================================================
# 7) Plot GTEx brain hit distribution
# =========================================================
summary_tissue_all <- unique(hits[SNP_type %in% c("LEAD", "PROXY"), .(SNP, Tissue, Tissue_label, SNP_type)])[
  , .(n_unique_SNP = uniqueN(SNP)), by = .(Tissue, Tissue_label, SNP_type)
]

summary_tissue_plot <- dcast(
  summary_tissue_all,
  Tissue + Tissue_label ~ SNP_type,
  value.var = "n_unique_SNP",
  fill = 0
)
if (!"LEAD" %in% names(summary_tissue_plot))  summary_tissue_plot[, LEAD := 0L]
if (!"PROXY" %in% names(summary_tissue_plot)) summary_tissue_plot[, PROXY := 0L]

summary_tissue_plot[, Total := LEAD + PROXY]
setorder(summary_tissue_plot, Total)

summary_tissue_long <- melt(
  summary_tissue_plot,
  id.vars = c("Tissue", "Tissue_label", "Total"),
  measure.vars = c("LEAD", "PROXY"),
  variable.name = "Type",
  value.name = "Count"
)
summary_tissue_long[, Tissue_label := factor(Tissue_label, levels = rev(summary_tissue_plot$Tissue_label))]
x_max <- max(summary_tissue_plot$Total) * 1.12

p_match <- ggplot(summary_tissue_long, aes(x = Tissue_label, y = Count, fill = Type)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("LEAD" = "#F8766D", "PROXY" = "#00BFC4"), breaks = c("LEAD", "PROXY")) +
  geom_text(
    data = summary_tissue_plot,
    aes(x = Tissue_label, y = Total, label = Total),
    inherit.aes = FALSE,
    hjust = -0.15,
    size = 3.6
  ) +
  scale_y_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "GTEx brain eQTL signals across tissues", x = NULL, y = "Unique SNPs") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    legend.title = element_blank(),
    legend.position = c(0.88, 0.48),
    legend.text = element_text(size = 11),
    plot.margin = margin(10, 40, 10, 10)
  )

match_plot_png <- file.path(dir_match, "Step3_GTExBrain_hits_distribution.png")
match_plot_pdf <- file.path(dir_match, "Step3_GTExBrain_hits_distribution.pdf")
ggsave(match_plot_png, p_match, width = 8.8, height = 6.4, dpi = 600, bg = "white")
ggsave(match_plot_pdf, p_match, device = cairo_pdf, width = 8.8, height = 6.4, bg = "white")

manifest_match <- data.table(
  file = c(
    "Step3_GTExBrain_hits_ALL.tsv.gz",
    "Step3_GTExBrain_hits_TOP1_perSNP_Tissue.tsv.gz",
    "Step3_GTExBrain_hits_summary_byTissue.tsv",
    "Step3_GTExBrain_hits_summary_byGene.tsv",
    "Step3_GTExBrain_hits_ALL_missingGene.tsv.gz",
    "Step3_GTExBrain_hits_TOP1_missingGene.tsv.gz",
    "Step3_GTExBrain_hits_ALL_missingSNPtype.tsv.gz",
    "Step3_GTExBrain_hits_debug_stats.tsv",
    "Step3_GTExBrain_hits_distribution.png",
    "Step3_GTExBrain_hits_distribution.pdf"
  ),
  full_path = c(
    all_hits_file,
    top1_hits_file,
    summary_by_tissue_file,
    summary_by_gene_file,
    file.path(dir_match, "Step3_GTExBrain_hits_ALL_missingGene.tsv.gz"),
    file.path(dir_match, "Step3_GTExBrain_hits_TOP1_missingGene.tsv.gz"),
    file.path(dir_match, "Step3_GTExBrain_hits_ALL_missingSNPtype.tsv.gz"),
    file.path(dir_match, "Step3_GTExBrain_hits_debug_stats.tsv"),
    match_plot_png,
    match_plot_pdf
  )
)
fwrite(manifest_match, file.path(dir_match, "Step3_GTExBrain_hits_manifest.tsv"), sep = "\t")

cat("Matching stage completed.\n")
cat("ALL hits         : ", nrow(hits), "\n", sep = "")
cat("TOP1 rows        : ", nrow(top1_hits), "\n", sep = "")
cat("Unique hit SNPs  : ", uniqueN(hits$SNP), "\n", sep = "")
cat("Unique hit genes : ", uniqueN(hits[Gene != '', Gene]), "\n\n", sep = "")

# =========================================================
# 8) Build brain-only GTEx evidence panels
# =========================================================
hits_panel <- copy(hits)
hits_panel[, Gene := as.character(Gene)]
hits_panel[is.na(Gene), Gene := ""]
hits_panel[, Gene_nonempty := Gene != ""]
hits_panel[, Gene_priority := fifelse(Gene_nonempty, 1L, 0L)]
hits_panel[, Tissue := as.character(Tissue)]
hits_panel[, Tissue_label := clean_gtex_tissue_label(Tissue)]

hits_panel <- hits_panel[SNP %chin% panel$SNP]
if (nrow(hits_panel) == 0) stop("No hits remain after intersecting with the Step2 panel.")
hits_panel <- hits_panel[grepl("^Brain_", Tissue)]
if (nrow(hits_panel) == 0) stop("No brain-only hits remain after tissue filtering.")
stopifnot(all(grepl("^Brain_", hits_panel$Tissue)))

hits_panel[, GTEx_evidence_class := classify_gtex_evidence(P)]
hits_panel[, logP := -log10(P)]

setorder(hits_panel, SNP, -Gene_priority, P, Tissue, Gene)
best_pick <- hits_panel[, .SD[1], by = .(SNP)]
best_pick[, evidence_source := "GTEx_brain"]

best_pick_missing_gene <- best_pick[is.na(Gene) | Gene == ""]
best_pick_with_gene    <- best_pick[!is.na(Gene) & Gene != ""]

keep_cols <- intersect(
  c("SNP","SNP_type","CHR","BP","A1","A2","eaf","EAF","beta","se","pval","INFO",
    "anchor","leadSNP","anchor_SNP","index_snp","r2","r2_with_anchor","proxy_of","N_total","type_rank"),
  names(panel)
)

panel2 <- unique(panel[, ..keep_cols])
if (!"SNP_type" %in% names(panel2)) panel2[, SNP_type := "UNKNOWN"]
panel2[is.na(SNP_type) | SNP_type == "", SNP_type := "UNKNOWN"]
panel2 <- merge(panel2, gwas_map, by = "SNP", all.x = TRUE)

snp_panel <- merge(
  panel2,
  best_pick,
  by = "SNP",
  all.x = TRUE,
  suffixes = c(".panel", ".gtex")
)

snp_panel[, has_gtex := !is.na(P)]
snp_panel[, GTEx_evidence_class := classify_gtex_evidence(P)]

CHR_panel <- if ("CHR.panel" %in% names(snp_panel)) as.integer(snp_panel[["CHR.panel"]]) else rep(NA_integer_, nrow(snp_panel))
BP_panel  <- if ("BP.panel"  %in% names(snp_panel)) as.integer(snp_panel[["BP.panel"]])  else rep(NA_integer_, nrow(snp_panel))
CHR_hit   <- if ("CHR.gtex" %in% names(snp_panel)) as.integer(snp_panel[["CHR.gtex"]]) else rep(NA_integer_, nrow(snp_panel))
BP_hit    <- if ("BP.gtex"  %in% names(snp_panel)) as.integer(snp_panel[["BP.gtex"]])  else rep(NA_integer_, nrow(snp_panel))
CHR_gwas_vec <- if ("CHR_gwas" %in% names(snp_panel)) as.integer(snp_panel$CHR_gwas) else rep(NA_integer_, nrow(snp_panel))
BP_gwas_vec  <- if ("BP_gwas"  %in% names(snp_panel)) as.integer(snp_panel$BP_gwas)  else rep(NA_integer_, nrow(snp_panel))

CHR_final <- fifelse(!is.na(CHR_panel), CHR_panel,
                     fifelse(!is.na(CHR_gwas_vec), CHR_gwas_vec, CHR_hit))
BP_final  <- fifelse(!is.na(BP_panel), BP_panel,
                     fifelse(!is.na(BP_gwas_vec), BP_gwas_vec, BP_hit))

A1_panel <- if ("A1.panel" %in% names(snp_panel)) as.character(snp_panel[["A1.panel"]]) else rep(NA_character_, nrow(snp_panel))
A2_panel <- if ("A2.panel" %in% names(snp_panel)) as.character(snp_panel[["A2.panel"]]) else rep(NA_character_, nrow(snp_panel))
A1_hit   <- if ("A1.gtex"  %in% names(snp_panel)) as.character(snp_panel[["A1.gtex"]])  else rep(NA_character_, nrow(snp_panel))
A2_hit   <- if ("A2.gtex"  %in% names(snp_panel)) as.character(snp_panel[["A2.gtex"]])  else rep(NA_character_, nrow(snp_panel))

A1_panel[A1_panel == ""] <- NA_character_
A2_panel[A2_panel == ""] <- NA_character_
A1_hit[A1_hit == ""] <- NA_character_
A2_hit[A2_hit == ""] <- NA_character_

A1_gwas_vec <- if ("A1_gwas" %in% names(snp_panel)) as.character(snp_panel$A1_gwas) else rep(NA_character_, nrow(snp_panel))
A2_gwas_vec <- if ("A2_gwas" %in% names(snp_panel)) as.character(snp_panel$A2_gwas) else rep(NA_character_, nrow(snp_panel))
A1_gwas_vec[A1_gwas_vec == ""] <- NA_character_
A2_gwas_vec[A2_gwas_vec == ""] <- NA_character_

A1_final <- fifelse(!is.na(A1_panel), A1_panel,
                    fifelse(!is.na(A1_gwas_vec), A1_gwas_vec, A1_hit))
A2_final <- fifelse(!is.na(A2_panel), A2_panel,
                    fifelse(!is.na(A2_gwas_vec), A2_gwas_vec, A2_hit))

EAF_panel_small <- if ("eaf" %in% names(snp_panel)) as.numeric(snp_panel$eaf) else rep(NA_real_, nrow(snp_panel))
EAF_panel_big   <- if ("EAF.panel" %in% names(snp_panel)) as.numeric(snp_panel[["EAF.panel"]]) else
  if ("EAF" %in% names(snp_panel)) as.numeric(snp_panel$EAF) else rep(NA_real_, nrow(snp_panel))
EAF_hit <- if ("EAF.gtex" %in% names(snp_panel)) as.numeric(snp_panel[["EAF.gtex"]]) else rep(NA_real_, nrow(snp_panel))
EAF_final <- fifelse(!is.na(EAF_panel_small), EAF_panel_small,
                     fifelse(!is.na(EAF_panel_big), EAF_panel_big, EAF_hit))
INFO_final <- if ("INFO" %in% names(snp_panel)) snp_panel$INFO else rep(NA, nrow(snp_panel))

if (!"SNP_type" %in% names(snp_panel)) {
  if ("SNP_type.panel" %in% names(snp_panel)) {
    snp_panel[, SNP_type := SNP_type.panel]
  } else if ("SNP_type.gtex" %in% names(snp_panel)) {
    snp_panel[, SNP_type := SNP_type.gtex]
  } else {
    snp_panel[, SNP_type := "UNKNOWN"]
  }
}
snp_panel[is.na(SNP_type) | SNP_type == "", SNP_type := "UNKNOWN"]

snp_panel_out <- snp_panel[, .(
  SNP,
  SNP_type,
  CHR = CHR_final,
  BP  = BP_final,
  A1  = A1_final,
  A2  = A2_final,
  EAF = EAF_final,
  Tissue = Tissue,
  Tissue_label = Tissue_label,
  Gene_GTEx = Gene,
  P_GTEx = P,
  GTEx_evidence_class,
  evidence_source,
  has_gtex,
  INFO = INFO_final
)]

hits2 <- merge(
  hits_panel,
  panel2[, .(SNP, SNP_type)],
  by = "SNP",
  all.x = TRUE,
  suffixes = c(".hit", ".panel")
)

if ("SNP_type.hit" %in% names(hits2) && "SNP_type.panel" %in% names(hits2)) {
  hits2[, SNP_type := fifelse(
    !is.na(SNP_type.hit) & SNP_type.hit != "",
    SNP_type.hit,
    SNP_type.panel
  )]
} else if ("SNP_type.hit" %in% names(hits2)) {
  hits2[, SNP_type := SNP_type.hit]
} else if ("SNP_type.panel" %in% names(hits2)) {
  hits2[, SNP_type := SNP_type.panel]
} else if (!"SNP_type" %in% names(hits2)) {
  hits2[, SNP_type := "UNKNOWN"]
}
hits2[is.na(SNP_type) | SNP_type == "", SNP_type := "UNKNOWN"]
setorder(hits2, Gene, P)

gene_support <- hits2[!is.na(Gene) & Gene != "",
                      .(
                        n_hits_support = .N,
                        n_unique_snp_support = uniqueN(SNP),
                        n_index = sum(SNP_type %in% c("INDEX","LEAD"), na.rm = TRUE),
                        n_proxy = sum(SNP_type == "PROXY", na.rm = TRUE),
                        top_p = min(P, na.rm = TRUE),
                        GTEx_evidence_class = classify_gtex_evidence(min(P, na.rm = TRUE)),
                        evidence_source = "GTEx_brain",
                        top_tissue = Tissue[which.min(P)],
                        top_tissue_label = Tissue_label[which.min(P)],
                        top_snp = SNP[which.min(P)]
                      ),
                      by = .(Gene)]
setorder(gene_support, top_p)

best_snp_loc <- unique(snp_panel_out[, .(SNP, CHR, BP, A1, A2, EAF)])
gene_panel <- merge(gene_support, best_snp_loc, by.x = "top_snp", by.y = "SNP", all.x = TRUE)

gene_panel_out <- gene_panel[, .(
  Gene,
  GTEx_evidence_class,
  top_p,
  evidence_source,
  top_tissue,
  top_tissue_label,
  top_snp,
  n_hits_support,
  n_unique_snp_support,
  n_index,
  n_proxy,
  CHR,
  BP,
  A1,
  A2,
  EAF
)]

# =========================================================
# 9) Save evidence panels
# =========================================================
count_gene <- gene_panel_out[, .N, by = .(GTEx_evidence_class)][order(GTEx_evidence_class)]
count_snp  <- snp_panel_out[,  .N, by = .(GTEx_evidence_class)][order(GTEx_evidence_class)]

summary_panel <- rbind(
  cbind(table = "GenePanel", count_gene),
  cbind(table = "SNPPanel",  count_snp),
  fill = TRUE
)

debug_missing_gene_snp    <- snp_panel_out[is.na(Gene_GTEx) | Gene_GTEx == ""]
debug_nonmissing_gene_snp <- snp_panel_out[!is.na(Gene_GTEx) & Gene_GTEx != ""]
debug_missing_A1 <- snp_panel_out[is.na(A1) | A1 == ""]
debug_missing_A2 <- snp_panel_out[is.na(A2) | A2 == ""]

debug_stats_panel <- data.table(
  n_panel_snps = nrow(panel2),
  n_hits_after_filter = nrow(hits_panel),
  n_best_pick = nrow(best_pick),
  n_best_pick_with_gene = nrow(best_pick_with_gene),
  n_best_pick_missing_gene = nrow(best_pick_missing_gene),
  n_snp_panel_total = nrow(snp_panel_out),
  n_snp_panel_with_gene = nrow(debug_nonmissing_gene_snp),
  n_snp_panel_missing_gene = nrow(debug_missing_gene_snp),
  prop_snp_panel_missing_gene = round(nrow(debug_missing_gene_snp) / max(1, nrow(snp_panel_out)), 4),
  n_snp_panel_A1_missing = nrow(debug_missing_A1),
  n_snp_panel_A2_missing = nrow(debug_missing_A2),
  prop_snp_panel_A1_missing = round(nrow(debug_missing_A1) / max(1, nrow(snp_panel_out)), 4),
  prop_snp_panel_A2_missing = round(nrow(debug_missing_A2) / max(1, nrow(snp_panel_out)), 4)
)

snp_panel_file   <- file.path(dir_tables, "Step3_GTExBrain_SNP_evidence_panel.tsv")
gene_panel_file  <- file.path(dir_tables, "Step3_GTExBrain_Gene_evidence_panel.tsv")
count_panel_file <- file.path(dir_tables, "Step3_GTExBrain_evidence_counts.tsv")

fwrite(snp_panel_out,  snp_panel_file,  sep = "\t")
fwrite(gene_panel_out, gene_panel_file, sep = "\t")
fwrite(summary_panel,  count_panel_file, sep = "\t")

fwrite(best_pick_missing_gene, file.path(dir_debug, "Step3_GTExBrain_best_pick_missing_gene.tsv"), sep = "\t")
fwrite(debug_missing_gene_snp, file.path(dir_debug, "Step3_GTExBrain_SNPPanel_missing_gene.tsv"), sep = "\t")
fwrite(debug_missing_A1,       file.path(dir_debug, "Step3_GTExBrain_SNPPanel_missing_A1.tsv"), sep = "\t")
fwrite(debug_missing_A2,       file.path(dir_debug, "Step3_GTExBrain_SNPPanel_missing_A2.tsv"), sep = "\t")
fwrite(debug_stats_panel,      file.path(dir_debug, "Step3_GTExBrain_panel_debug_stats.tsv"), sep = "\t")

# =========================================================
# 10) Save panel plots
# =========================================================
bp <- best_pick[!is.na(P)]
best_tissue_sum <- bp[, .(
  n_best = .N,
  n_unique_gene = uniqueN(Gene[Gene != ""])
), by = .(Tissue, Tissue_label)][order(-n_best)]

best_tissue_sum[, Tissue_label := factor(Tissue_label, levels = rev(best_tissue_sum$Tissue_label))]

p_best_tissue <- ggplot(best_tissue_sum, aes(x = Tissue_label, y = n_best)) +
  geom_col(fill = "#666666") +
  coord_flip() +
  labs(title = "Top GTEx brain tissue by SNP-level eQTL evidence", x = "Brain tissue", y = "Count of SNPs") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 13)
  )

ggsave(
  file.path(dir_plots, "Step3_GTExBrain_best_record_tissue_distribution.png"),
  p_best_tissue,
  width = 9,
  height = max(4, 0.32 * nrow(best_tissue_sum)),
  dpi = 300,
  bg = "white"
)

topN <- min(30, nrow(gene_panel_out))
gene_panel_top <- copy(gene_panel_out[1:topN])
gene_panel_top[, GTEx_evidence_class := factor(
  GTEx_evidence_class,
  levels = c("Strong", "Moderate", "Weak", "NoEvidence")
)]

p_top_genes <- ggplot(
  gene_panel_top,
  aes(x = reorder(Gene, -log10(top_p)), y = -log10(top_p), fill = GTEx_evidence_class)
) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Strong" = "#F8766D", "Moderate" = "#00BA38", "Weak" = "#619CFF", "NoEvidence" = "grey80"),
    drop = FALSE
  ) +
  labs(
    title = paste0("Top genes by best GTEx brain eQTL P value (Top ", topN, ")"),
    x = "Gene", y = "-log10(top P)", fill = "GTEx evidence"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

ggsave(
  file.path(dir_plots, "Step3_GTExBrain_top_genes_by_bestP.png"),
  p_top_genes,
  width = 8.5,
  height = max(4, 0.25 * topN),
  dpi = 300,
  bg = "white"
)

manifest_panels <- data.table(
  category = c(rep("tables", 3), rep("debug", 5), rep("plots", 2)),
  file = c(
    "Step3_GTExBrain_SNP_evidence_panel.tsv",
    "Step3_GTExBrain_Gene_evidence_panel.tsv",
    "Step3_GTExBrain_evidence_counts.tsv",
    "Step3_GTExBrain_best_pick_missing_gene.tsv",
    "Step3_GTExBrain_SNPPanel_missing_gene.tsv",
    "Step3_GTExBrain_SNPPanel_missing_A1.tsv",
    "Step3_GTExBrain_SNPPanel_missing_A2.tsv",
    "Step3_GTExBrain_panel_debug_stats.tsv",
    "Step3_GTExBrain_best_record_tissue_distribution.png",
    "Step3_GTExBrain_top_genes_by_bestP.png"
  ),
  full_path = c(
    snp_panel_file,
    gene_panel_file,
    count_panel_file,
    file.path(dir_debug, "Step3_GTExBrain_best_pick_missing_gene.tsv"),
    file.path(dir_debug, "Step3_GTExBrain_SNPPanel_missing_gene.tsv"),
    file.path(dir_debug, "Step3_GTExBrain_SNPPanel_missing_A1.tsv"),
    file.path(dir_debug, "Step3_GTExBrain_SNPPanel_missing_A2.tsv"),
    file.path(dir_debug, "Step3_GTExBrain_panel_debug_stats.tsv"),
    file.path(dir_plots, "Step3_GTExBrain_best_record_tissue_distribution.png"),
    file.path(dir_plots, "Step3_GTExBrain_top_genes_by_bestP.png")
  )
)
fwrite(manifest_panels, file.path(dir_panels, "Step3_GTExBrain_panel_manifest.tsv"), sep = "\t")

message("03_eQTL_annotation.R finished successfully.")
