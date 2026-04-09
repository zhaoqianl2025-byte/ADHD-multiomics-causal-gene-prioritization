rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# =========================================================
# Step 5.7
# Draw publication-style SMR regional association plot
# Current version: Panel A only for LSM6
# =========================================================


options(stringsAsFactors = FALSE)

# =========================================================
# 0) ROOT PATHS
# =========================================================
project_root <- "../.."

gwas_file <- file.path(
  project_root, "results", "Step5", "01_prepare_gwas", "Step5_ADHD2022_for_SMR.txt"
)

gene_list <- file.path(
  project_root, "data", "Step5", "glist-hg19.txt"
)

# =========================================================
# 1) SINGLE-GENE META
# Only modify this block when switching to another gene
# =========================================================
gene_symbol  <- "LSM6"
probe_id     <- "ENSG00000164167.12"
study1_label <- "BrainMeta study"
folder       <- "03_LSM6"

gene_dir <- file.path(
  project_root, "results", "Step5", "05_plots", "genes", folder
)

plot_dir <- file.path(gene_dir, "plot")

publication_dir <- file.path(
  project_root, "results", "Step5", "05_plots", "publication_plots", folder
)
if (!dir.exists(publication_dir)) dir.create(publication_dir, recursive = TRUE)

out_prefix <- file.path(
  publication_dir,
  paste0("Step5_", gene_symbol, "_regional_plot")
)

# =========================================================
# 2) PARAMETERS
# =========================================================
smr_sig_p_line <- 5e-6
max_probe_labels <- 1

show_mode <- "zoom"
zoom_center_mode <- "probe"
zoom_half_width_mb <- 0.18

# =========================================================
# 3) HELPER FUNCTIONS
# =========================================================
pretty_cap <- function(x, min_cap = 5) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(min_cap)
  m <- max(x, na.rm = TRUE)
  if (!is.finite(m)) return(min_cap)
  if (m <= min_cap) return(min_cap)
  base <- 10^(floor(log10(m)))
  ceiling(m / base) * base
}

format_p_sci <- function(p, digits = 2) {
  if (is.na(p) || !is.finite(p) || p <= 0) return("NA")
  format(p, scientific = TRUE, digits = digits)
}

read_block <- function(lines, tag) {
  idx <- which(startsWith(trimws(lines), paste0("$", tag)))[1]
  if (is.na(idx)) return(character())
  next_idx <- which(seq_along(lines) > idx & startsWith(trimws(lines), "$"))
  end <- if (length(next_idx)) min(next_idx) - 1 else length(lines)
  if (idx + 1 > end) return(character())
  out <- lines[(idx + 1):end]
  out[nzchar(trimws(out))]
}

find_plot_file <- function(gene_dir, plot_dir, gene_symbol, probe_id) {
  cand <- c(
    file.path(plot_dir, paste0("plot_", gene_symbol, ".", probe_id, ".txt")),
    file.path(gene_dir, paste0("plot_", gene_symbol, ".", probe_id, ".txt")),
    file.path(plot_dir, paste0("plot_", gene_symbol, ".txt")),
    file.path(gene_dir, paste0("plot_", gene_symbol, ".txt"))
  )
  
  hit <- cand[file.exists(cand)]
  if (length(hit) > 0) return(hit[1])
  
  if (dir.exists(plot_dir)) {
    patt <- paste0("^plot_", gsub("([.\\-+()\\[\\]{}^$|?*\\\\])", "\\\\\\1", gene_symbol), ".*\\.txt$")
    fs <- list.files(plot_dir, pattern = patt, full.names = TRUE)
    if (length(fs) > 0) return(fs[1])
  }
  
  if (dir.exists(gene_dir)) {
    patt <- paste0("^plot_", gsub("([.\\-+()\\[\\]{}^$|?*\\\\])", "\\\\\\1", gene_symbol), ".*\\.txt$")
    fs <- list.files(gene_dir, pattern = patt, full.names = TRUE)
    if (length(fs) > 0) return(fs[1])
  }
  
  return(NA_character_)
}

parse_smr_plot <- function(infile) {
  lines <- readLines(infile, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(trimws(lines))]
  
  header_line <- lines[1]
  target_probe <- sub("^\\$probe\\s+\\d+\\s+(\\S+).*", "\\1", header_line)
  
  probe_lines <- read_block(lines, "probe")
  snp_lines   <- read_block(lines, "SNP")
  eqtl_lines  <- read_block(lines, "eQTL")
  
  if (length(probe_lines) == 0) stop("Probe block not found: ", infile)
  if (length(snp_lines)   == 0) stop("SNP block not found: ", infile)
  if (length(eqtl_lines)  == 0) stop("eQTL block not found: ", infile)
  
  probe_dt <- fread(text = paste(probe_lines, collapse = "\n"), header = FALSE, fill = TRUE)
  while (ncol(probe_dt) < 9) probe_dt[, paste0("V", ncol(probe_dt) + 1) := NA]
  probe_dt <- probe_dt[, 1:9]
  setnames(probe_dt, c("probeID", "chr", "probe_bp", "gene", "gene_start", "gene_end", "strand", "p_SMR", "p_HEIDI"))
  
  probe_dt[, chr := suppressWarnings(as.integer(chr))]
  probe_dt[, probe_bp := suppressWarnings(as.numeric(probe_bp))]
  probe_dt[, gene_start := suppressWarnings(as.numeric(gene_start))]
  probe_dt[, gene_end   := suppressWarnings(as.numeric(gene_end))]
  probe_dt[, p_SMR      := suppressWarnings(as.numeric(p_SMR))]
  probe_dt[, p_HEIDI    := suppressWarnings(as.numeric(p_HEIDI))]
  probe_dt[, probeID    := as.character(probeID)]
  probe_dt[, gene       := as.character(gene)]
  probe_dt[, strand     := as.character(strand)]
  
  probe_dt[, has_smr := !is.na(p_SMR) & is.finite(p_SMR) & p_SMR > 0]
  probe_dt[, logp_smr := fifelse(has_smr, -log10(p_SMR), NA_real_)]
  
  snpinfo_dt <- fread(text = paste(snp_lines, collapse = "\n"), header = FALSE, fill = TRUE)
  if (ncol(snpinfo_dt) < 5) stop("SNP block has fewer than 5 columns.")
  setnames(snpinfo_dt, names(snpinfo_dt)[1:5], c("snp", "chr", "bp", "A1", "A2"))
  snpinfo_dt[, snp := as.character(snp)]
  snpinfo_dt[, chr := suppressWarnings(as.integer(chr))]
  snpinfo_dt[, bp  := suppressWarnings(as.numeric(bp))]
  setorder(snpinfo_dt, snp, chr, bp)
  snpinfo_dt <- snpinfo_dt[, .SD[1], by = snp]
  setkey(snpinfo_dt, snp)
  
  eqtl_dt <- fread(text = paste(eqtl_lines, collapse = "\n"), header = FALSE, fill = TRUE)
  if (ncol(eqtl_dt) < 3) stop("eQTL block has fewer than 3 columns.")
  
  if (ncol(eqtl_dt) == 3) {
    setnames(eqtl_dt, names(eqtl_dt)[1:3], c("snp", "b_eqtl", "se_eqtl"))
    eqtl_dt[, b_eqtl := as.numeric(b_eqtl)]
    eqtl_dt[, se_eqtl := as.numeric(se_eqtl)]
    eqtl_dt[, p_eqtl := 2 * pnorm(-abs(b_eqtl / se_eqtl))]
  } else {
    setnames(eqtl_dt, names(eqtl_dt)[1:4], c("snp", "b_eqtl", "se_eqtl", "p_eqtl"))
    eqtl_dt[, b_eqtl := as.numeric(b_eqtl)]
    eqtl_dt[, se_eqtl := as.numeric(se_eqtl)]
    eqtl_dt[, p_eqtl_raw := as.character(p_eqtl)]
    eqtl_dt[, p_eqtl := fifelse(
      is.na(p_eqtl_raw) | p_eqtl_raw %in% c("NA", "NaN", ""),
      2 * pnorm(-abs(b_eqtl / se_eqtl)),
      as.numeric(p_eqtl_raw)
    )]
    eqtl_dt[, p_eqtl_raw := NULL]
  }
  
  eqtl_dt[, snp := as.character(snp)]
  setkey(eqtl_dt, snp)
  eqtl_dt <- snpinfo_dt[eqtl_dt]
  eqtl_dt <- eqtl_dt[!is.na(bp) & !is.na(p_eqtl) & p_eqtl > 0]
  eqtl_dt[, logp_eqtl := -log10(p_eqtl)]
  
  target_row <- probe_dt[probeID == target_probe]
  if (nrow(target_row) == 0) {
    target_row <- probe_dt[has_smr == TRUE][which.min(p_SMR)]
  }
  if (nrow(target_row) == 0) {
    stop("No target probe found and no valid SMR values available.")
  }
  
  list(
    probe_dt     = probe_dt,
    target_probe = target_probe,
    target_row   = target_row[1],
    snpinfo_dt   = snpinfo_dt,
    eqtl_dt      = eqtl_dt
  )
}

build_gwas_dt <- function(gwas_file, snpinfo_dt) {
  gwas_all <- fread(gwas_file)
  setnames(gwas_all, tolower(names(gwas_all)))
  stopifnot(all(c("snp", "b", "se", "p") %in% names(gwas_all)))
  
  gwas_dt <- gwas_all[snp %in% snpinfo_dt$snp, .(
    snp,
    b_gwas = as.numeric(b),
    se_gwas = as.numeric(se),
    p_gwas = as.numeric(p)
  )]
  
  gwas_dt[, snp := as.character(snp)]
  setkey(gwas_dt, snp)
  setkey(snpinfo_dt, snp)
  gwas_dt <- snpinfo_dt[gwas_dt]
  gwas_dt <- gwas_dt[!is.na(bp) & !is.na(p_gwas) & p_gwas > 0]
  gwas_dt[, logp_gwas := -log10(p_gwas)]
  gwas_dt
}

build_gene_track <- function(gene_list, region_chr, xlim_bp, lead_bp) {
  glist <- fread(gene_list, header = FALSE, fill = TRUE)
  stopifnot(ncol(glist) >= 4)
  setnames(glist, 1:4, c("chr", "start", "end", "gene"))
  
  glist[, chr   := suppressWarnings(as.integer(chr))]
  glist[, start := suppressWarnings(as.numeric(start))]
  glist[, end   := suppressWarnings(as.numeric(end))]
  glist[, gene  := as.character(gene)]
  
  gtrack <- glist[
    !is.na(chr) & !is.na(start) & !is.na(end) &
      chr == region_chr &
      end >= xlim_bp[1] &
      start <= xlim_bp[2]
  ]
  
  gtrack <- unique(gtrack, by = c("gene", "start", "end"))
  
  if (nrow(gtrack) > 0) {
    gtrack[, mid := (start + end) / 2]
    gtrack[, dist_to_lead := abs(mid - lead_bp)]
    setorder(gtrack, dist_to_lead)
    gtrack <- gtrack[1:min(8, .N)]
    setorder(gtrack, start)
    gtrack[, lane := (seq_len(.N) - 1) %% 4]
  }
  
  gtrack
}

choose_probe_labels <- function(probe_dt, target_probe, max_n = 1) {
  dt <- copy(probe_dt)
  dt[, label_show := fifelse(!is.na(gene) & gene != "", gene, probeID)]
  dt <- dt[has_smr == TRUE & !is.na(probe_bp)]
  if (nrow(dt) == 0) return(dt)
  
  dt[, x_mb := probe_bp / 1e6]
  out <- dt[probeID == target_probe]
  if (nrow(out) == 0) {
    setorder(dt, p_SMR)
    out <- dt[1]
  }
  out
}

save_plot <- function(p, filename, w, h) {
  ggsave(filename = paste0(filename, ".png"), plot = p, width = w, height = h, dpi = 350)
  ggsave(filename = paste0(filename, ".pdf"), plot = p, width = w, height = h)
}

# =========================================================
# 4) MAIN WORKFLOW: PANEL A FOR LSM6 ONLY
# =========================================================
if (!file.exists(gwas_file)) stop("GWAS file does not exist: ", gwas_file)
if (!file.exists(gene_list)) stop("Gene list file does not exist: ", gene_list)
if (!dir.exists(gene_dir)) stop("Gene directory does not exist: ", gene_dir)

in_plot1 <- find_plot_file(
  gene_dir    = gene_dir,
  plot_dir    = plot_dir,
  gene_symbol = gene_symbol,
  probe_id    = probe_id
)

if (is.na(in_plot1)) stop("Input plot file does not exist for: ", gene_symbol)

message("Using plot file: ", in_plot1)

obj1 <- parse_smr_plot(in_plot1)
gwas_dt_full <- build_gwas_dt(gwas_file, obj1$snpinfo_dt)

if (nrow(gwas_dt_full) == 0) stop("No matched GWAS SNPs found for: ", gene_symbol)

probe_label1 <- obj1$target_probe
region_chr <- gwas_dt_full$chr[1]
smr_obs_p <- obj1$target_row$p_SMR[1]

lead_snp <- obj1$eqtl_dt[which.min(p_eqtl)]$snp[1]
lead_bp  <- obj1$eqtl_dt[snp == lead_snp]$bp[1]
lead_mb  <- lead_bp / 1e6

probe_bp_center <- obj1$target_row$probe_bp[1]
probe_mb_center <- probe_bp_center / 1e6

if (show_mode == "zoom") {
  center_mb <- if (zoom_center_mode == "probe") probe_mb_center else lead_mb
  xlim_mb <- c(center_mb - zoom_half_width_mb, center_mb + zoom_half_width_mb)
  xlim_bp <- xlim_mb * 1e6
} else {
  region_min <- min(gwas_dt_full$bp, na.rm = TRUE)
  region_max <- max(gwas_dt_full$bp, na.rm = TRUE)
  pad <- (region_max - region_min) * 0.02
  xlim_bp <- c(region_min - pad, region_max + pad)
  xlim_mb <- xlim_bp / 1e6
}

gwas_dt_full[, x_mb := bp / 1e6]
obj1$eqtl_dt[, x_mb := bp / 1e6]
obj1$probe_dt[, x_mb := probe_bp / 1e6]

gwas_dt <- copy(gwas_dt_full)[x_mb >= xlim_mb[1] & x_mb <= xlim_mb[2]]
obj1_eqtl_zoom  <- copy(obj1$eqtl_dt)[x_mb >= xlim_mb[1] & x_mb <= xlim_mb[2]]
obj1_probe_zoom <- copy(obj1$probe_dt)[x_mb >= xlim_mb[1] & x_mb <= xlim_mb[2]]

if (nrow(gwas_dt) == 0) stop("No GWAS points in zoom window for: ", gene_symbol)
if (nrow(obj1_eqtl_zoom) == 0) stop("No eQTL points in zoom window for: ", gene_symbol)

gtrack <- build_gene_track(gene_list, region_chr, xlim_bp, lead_bp)
if (nrow(gtrack) > 0) {
  gtrack[, xstart_mb := start / 1e6]
  gtrack[, xend_mb   := end / 1e6]
}

probe_lab1 <- choose_probe_labels(obj1_probe_zoom, probe_label1, max_n = max_probe_labels)
probe_smr1   <- obj1_probe_zoom[has_smr == TRUE]
probe_nosmr1 <- obj1_probe_zoom[has_smr == FALSE]
target_probe_dt1 <- obj1_probe_zoom[probeID == probe_label1]
top_eqtl1 <- obj1_eqtl_zoom[which.min(p_eqtl)]

# =========================================================
# 5) PANEL A
# =========================================================
ymax_gwas <- pretty_cap(c(gwas_dt$logp_gwas, probe_smr1$logp_smr), min_cap = 8)
y_probe_text <- ymax_gwas * 0.90

p_a1 <- ggplot() +
  geom_point(
    data = gwas_dt,
    aes(x = x_mb, y = logp_gwas),
    size = 1.25, color = "#9a7a2f"
  ) +
  geom_point(
    data = probe_smr1,
    aes(x = x_mb, y = logp_smr),
    shape = 23, size = 2.2, fill = "white",
    color = "#d33682", stroke = 0.6
  ) +
  {
    if (nrow(probe_nosmr1) > 0) {
      geom_point(
        data = probe_nosmr1,
        aes(x = x_mb, y = 0.35),
        shape = 24, size = 2.4, fill = "white",
        color = "black", stroke = 0.6
      )
    }
  } +
  geom_point(
    data = target_probe_dt1[target_probe_dt1$has_smr == TRUE],
    aes(x = x_mb, y = logp_smr),
    shape = 23, size = 3.2, fill = "white",
    color = "#d33682", stroke = 0.95
  ) +
  geom_hline(
    yintercept = -log10(smr_sig_p_line),
    linetype = "dashed", color = "#d33682", linewidth = 0.5
  ) +
  geom_point(
    data = top_eqtl1,
    aes(x = x_mb, y = 0.35),
    inherit.aes = FALSE,
    shape = 24, size = 3.1, fill = "white", color = "black", stroke = 0.7
  ) +
  geom_vline(
    data = probe_lab1,
    aes(xintercept = x_mb),
    inherit.aes = FALSE,
    linetype = "dotted", color = "grey55", linewidth = 0.35
  ) +
  geom_text(
    data = probe_lab1,
    aes(x = x_mb, y = y_probe_text, label = label_show),
    inherit.aes = FALSE,
    angle = 45,
    hjust = 0,
    vjust = 0,
    size = 3.4,
    color = "#d33682",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = xlim_mb[1] + 0.06 * diff(xlim_mb),
    y = ymax_gwas * 0.78,
    label = gene_symbol,
    fontface = "italic", hjust = 0, size = 6.0
  ) +
  annotate(
    "text",
    x = xlim_mb[2],
    y = -log10(smr_sig_p_line) + 0.15,
    label = paste0("P[SMR] = ", format_p_sci(smr_obs_p, digits = 2)),
    hjust = 1, color = "#d33682", size = 4
  ) +
  coord_cartesian(xlim = xlim_mb, ylim = c(0, ymax_gwas * 1.02)) +
  theme_classic(base_size = 12) +
  labs(y = expression(-log[10](P[GWAS]~or~P[SMR])), x = NULL) +
  theme(
    axis.title.x = element_blank(),
    plot.margin = margin(5, 15, 2, 10)
  )

cap_eqtl1 <- quantile(obj1_eqtl_zoom$logp_eqtl, 0.99, na.rm = TRUE)
cap_eqtl1 <- max(12, min(35, cap_eqtl1))

p_a2 <- ggplot(obj1_eqtl_zoom, aes(x = x_mb, y = logp_eqtl)) +
  geom_point(shape = 4, size = 1.2, stroke = 0.45, color = "#c2185b") +
  annotate(
    "text",
    x = xlim_mb[1] + 0.02 * diff(xlim_mb), y = cap_eqtl1 * 0.82,
    label = paste0(study1_label, "\n", probe_label1, " (", gene_symbol, ")"),
    hjust = 0, size = 4, fontface = "bold"
  ) +
  coord_cartesian(xlim = xlim_mb, ylim = c(0, cap_eqtl1)) +
  theme_classic(base_size = 12) +
  labs(y = expression(-log[10](P[eQTL])), x = NULL) +
  theme(
    axis.title.x = element_blank(),
    plot.margin = margin(2, 15, 2, 10)
  )

if (nrow(gtrack) > 0) {
  p_agene <- ggplot(gtrack) +
    geom_segment(
      aes(x = xstart_mb, xend = xend_mb, y = lane, yend = lane),
      linewidth = 0.45, color = "grey35"
    ) +
    geom_point(aes(x = xstart_mb, y = lane), size = 0.8, color = "grey35") +
    geom_point(aes(x = xend_mb, y = lane), size = 0.8, color = "grey35") +
    geom_text(
      aes(x = (xstart_mb + xend_mb) / 2, y = lane + 0.16, label = gene),
      size = 2.5, fontface = "italic", color = "grey20",
      check_overlap = TRUE
    ) +
    coord_cartesian(xlim = xlim_mb, ylim = c(-0.4, 3.5), clip = "off") +
    theme_void(base_size = 11) +
    labs(x = paste0("Chromosome ", region_chr, " position (Mb)")) +
    theme(
      plot.margin = margin(0, 15, 6, 10),
      axis.title.x = element_text(size = 11)
    )
} else {
  p_agene <- ggplot() +
    theme_void() +
    coord_cartesian(xlim = xlim_mb)
}

panel_a <- p_a1 / p_a2 / p_agene + plot_layout(heights = c(1.25, 1, 0.52))

# =========================================================
# 6) SAVE PANEL A ONLY
# =========================================================
save_plot(
  panel_a,
  paste0(out_prefix, "_panelA"),
  10,
  6.8
)

cat("Saved Panel A for ", gene_symbol, ":\n", sep = "")
cat(paste0(out_prefix, "_panelA.png"), "\n")
cat(paste0(out_prefix, "_panelA.pdf"), "\n")
cat("xlim_mb = ", paste(round(xlim_mb, 4), collapse = " ~ "), "\n", sep = "")

cat("\nDone for LSM6 Panel A only.\n")