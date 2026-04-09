rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

# =========================================================
# Step 6: Multi-evidence integration
# =========================================================

# =========================================================
# 0) PATHS
# =========================================================
project_root <- ".."

step3_gene <- file.path(project_root, "results", "Step3", "panels", "tables", "Step3_GTExBrain_Gene_evidence_panel.tsv")
step3_snp  <- file.path(project_root, "results", "Step3", "panels", "tables", "Step3_GTExBrain_SNP_evidence_panel.tsv")
coloc_file <- file.path(project_root, "results", "Step4", "Step4_GTExBrain_coloc_gene_best.tsv")
adhd_gwas  <- file.path(project_root, "data", "Step1", "ADHD2022_iPSYCH_deCODE_PGC.meta.tsv")
smr_allchr_txt <- file.path(project_root, "results", "Step5", "03_merge", "Step5_SMR_BrainMeta_allchr.tsv")
smr_chr_glob   <- file.path(project_root, "results", "Step5", "02_smr_raw", "Step5_SMR_chr*.smr")

out_root <- file.path(project_root, "results", "Step6")
dir_tables  <- file.path(out_root, "tables")
dir_figures <- file.path(out_root, "figures")
dir_debug   <- file.path(out_root, "debug")
dir_summary <- file.path(out_root, "summary")

dir.create(out_root,    showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tables,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_figures, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_debug,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_summary, showWarnings = FALSE, recursive = TRUE)

cat("Running Step 6 integration...
")

# =========================================================
# 1) PARAMETERS
# =========================================================
WINDOW <- 5e5
WINDOW_BIN <- 5e5
SMR_P_THR   <- 5e-6
HEIDI_P_THR <- 0.01
VALID_BASES <- c("A", "C", "G", "T")

# =========================================================
# 2) UTILITIES
# =========================================================
rename_if_exists <- function(dt, old, new) {
  if (old %in% names(dt) && !(new %in% names(dt))) {
    setnames(dt, old, new)
  }
}

drop_if_exists <- function(dt, cols) {
  cols <- intersect(cols, names(dt))
  if (length(cols) > 0) dt[, (cols) := NULL]
  invisible(dt)
}

pick1 <- function(cands, nms) {
  hit <- intersect(cands, nms)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

rank_coloc_tier <- function(x) {
  x <- as.character(x)
  fifelse(x == "STRONG", 3L,
          fifelse(x == "MID", 2L,
                  fifelse(x == "WEAK", 1L, 0L)))
}

brain_tissue_map <- c(
  "Brain_Amygdala" = "Amygdala",
  "Brain_Anterior_cingulate_cortex_BA24" = "ACC (BA24)",
  "Brain_Caudate_basal_ganglia" = "Caudate",
  "Brain_Cerebellar_Hemisphere" = "Cerebellar Hemisphere",
  "Brain_Cerebellum" = "Cerebellum",
  "Brain_Cortex" = "Cortex",
  "Brain_Frontal_Cortex_BA9" = "Frontal Cortex (BA9)",
  "Brain_Hippocampus" = "Hippocampus",
  "Brain_Hypothalamus" = "Hypothalamus",
  "Brain_Nucleus_accumbens_basal_ganglia" = "Nucleus accumbens",
  "Brain_Putamen_basal_ganglia" = "Putamen",
  "Brain_Spinal_cord_cervical_c-1" = "Spinal cord (C1)",
  "Brain_Substantia_nigra" = "Substantia nigra"
)

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
  x
}

abbr_tissue <- function(x) {
  x <- as.character(x)
  x2 <- gsub("\\.v10\\.eQTLs$", "", x)
  x2 <- gsub("\\.txt$", "", x2)
  out <- unname(brain_tissue_map[x2])
  out[is.na(out)] <- clean_gtex_tissue_label(x[is.na(out)])
  out
}

map_coloc_level <- function(x) {
  x <- as.character(x)
  fifelse(x == "STRONG", "Strong",
          fifelse(x == "MID", "Moderate",
                  fifelse(x == "WEAK", "Weak", "NoEvidence")))
}

is_strict_snp <- function(a1, a2) {
  !is.na(a1) & !is.na(a2) &
    nchar(a1) == 1 & nchar(a2) == 1 &
    a1 %in% VALID_BASES & a2 %in% VALID_BASES
}

fmt_p <- function(x) {
  ifelse(is.na(x), NA_character_, format(x, scientific = TRUE, digits = 3))
}

# =========================================================
# 3) CHECK INPUT FILES
# =========================================================
need_files <- c(step3_gene, step3_snp, coloc_file, adhd_gwas)
miss_files <- need_files[!file.exists(need_files)]
if (length(miss_files) > 0) stop("Missing input files:\n", paste(miss_files, collapse = "\n"))

# =========================================================
# 4) LOAD STEP 3 PANELS
# =========================================================
g_gene <- fread(step3_gene)
g_snp  <- fread(step3_snp)

old_gene_fields <- c("Final_Level_GTEx", "bestP", "bestLevel_GTEx", "bestTissue", "bestSNP")
old_snp_fields  <- c("Final_Level_GTEx", "level_GTEx", "GTEx_Tissue", "GTEx_Gene", "GTEx_P")

if (any(old_gene_fields %in% names(g_gene))) {
  stop("Old GTEx naming detected in Step3 gene panel. Please regenerate Step3 outputs using the updated Step3 script.")
}
if (any(old_snp_fields %in% names(g_snp))) {
  stop("Old GTEx naming detected in Step3 SNP panel. Please regenerate Step3 outputs using the updated Step3 script.")
}

stopifnot(all(c("Gene","GTEx_evidence_class","top_p","top_tissue","top_snp") %in% names(g_gene)))
stopifnot(all(c("SNP","Gene_GTEx","P_GTEx","Tissue","GTEx_evidence_class") %in% names(g_snp)))

g_gene[, bestP_GTEx := as.numeric(top_p)]
g_gene[, bestSNP := as.character(top_snp)]
g_gene[, bestTissue := as.character(top_tissue)]
g_gene[, bestLevel_GTEx := as.character(GTEx_evidence_class)]

if (!("evidence_source" %in% names(g_gene))) g_gene[, evidence_source := "GTEx_brain"]
if (!("evidence_source" %in% names(g_snp)))  g_snp[, evidence_source := "GTEx_brain"]
if (!("Tissue_label" %in% names(g_snp))) g_snp[, Tissue_label := abbr_tissue(Tissue)]

# =========================================================
# 5) LOAD COLOC
# =========================================================
coloc <- fread(coloc_file)
rename_if_exists(coloc, "id", "Gene")
rename_if_exists(coloc, "PP.H4.abf", "PP_H4")

stopifnot("Gene" %in% names(coloc))
stopifnot("PP_H4" %in% names(coloc))
coloc[, PP_H4 := as.numeric(PP_H4)]

if (!("coloc_tier" %in% names(coloc))) {
  coloc[, coloc_tier := fifelse(PP_H4 >= 0.80, "STRONG",
                                fifelse(PP_H4 >= 0.50, "MID", "WEAK"))]
}

tissue_col <- pick1(c("tissue", "bestTissue", "Tissue", "GTEx_Tissue"), names(coloc))
lead_col   <- pick1(c("locus_leadSNP", "leadSNP", "lead_snp", "lead"), names(coloc))

if (!is.na(tissue_col) && tissue_col != "coloc_bestTissue") setnames(coloc, tissue_col, "coloc_bestTissue")
if (!is.na(lead_col) && lead_col != "coloc_leadSNP") setnames(coloc, lead_col, "coloc_leadSNP")

setorderv(coloc, cols = c("Gene", "PP_H4"), order = c(1, -1))
coloc_gene <- coloc[, .SD[1], by = Gene][, .(
  Gene,
  PP_H4,
  coloc_tier,
  coloc_bestTissue = if ("coloc_bestTissue" %in% names(coloc)) as.character(coloc_bestTissue) else NA_character_,
  coloc_leadSNP    = if ("coloc_leadSNP" %in% names(coloc)) as.character(coloc_leadSNP) else NA_character_
)]

# =========================================================
# 6) LOAD SMR
# =========================================================
smr <- NULL
if (file.exists(smr_allchr_txt)) {
  smr <- fread(smr_allchr_txt)
} else {
  smr_files <- Sys.glob(smr_chr_glob)
  if (length(smr_files) == 0) {
    stop("No SMR inputs found: neither Step5_SMR_BrainMeta_allchr.tsv nor Step5_SMR_chr*.smr")
  }
  smr <- rbindlist(lapply(smr_files, fread), fill = TRUE)
}

rename_if_exists(smr, "probeID", "Gene")
rename_if_exists(smr, "gene", "Gene")
rename_if_exists(smr, "topSNP", "SMR_topSNP")
rename_if_exists(smr, "topSNP_chr", "SMR_topSNP_chr")
rename_if_exists(smr, "topSNP_bp", "SMR_topSNP_bp")
rename_if_exists(smr, "p_heidi", "p_HEIDI")

need_smr_cols <- c("Gene", "p_SMR", "p_HEIDI")
if (!all(need_smr_cols %in% names(smr))) {
  stop("SMR table missing required columns: ", paste(setdiff(need_smr_cols, names(smr)), collapse = ", "))
}

smr[, p_SMR   := as.numeric(p_SMR)]
smr[, p_HEIDI := as.numeric(p_HEIDI)]
smr[, SMR_sig           := !is.na(p_SMR)   & (p_SMR < SMR_P_THR)]
smr[, HEIDI_pass        := !is.na(p_HEIDI) & (p_HEIDI > HEIDI_P_THR)]
smr[, SMR_HEIDI_support := SMR_sig & HEIDI_pass]

smr_best <- smr[order(-as.integer(SMR_HEIDI_support), p_SMR)][, .SD[1], by = .(Gene)]
smr_best <- smr_best[, .(
  Gene,
  SMR_topSNP,
  SMR_topSNP_chr,
  SMR_topSNP_bp,
  p_SMR,
  p_HEIDI,
  SMR_sig,
  HEIDI_pass,
  SMR_HEIDI_support
)]

# =========================================================
# 7) LOAD GWAS MAP
# =========================================================
gwas_map <- fread(
  adhd_gwas,
  select = c("SNP", "CHR", "BP", "A1", "A2", "FRQ_A_38691", "FRQ_U_186843", "INFO", "OR", "SE", "P")
)
setnames(gwas_map, c("FRQ_A_38691", "FRQ_U_186843"), c("FRQ_A_gwas", "FRQ_U_gwas"))

gwas_map[, SNP := as.character(SNP)]
gwas_map[, CHR := as.integer(CHR)]
gwas_map[, BP  := as.integer(BP)]
gwas_map[, A1  := as.character(A1)]
gwas_map[, A2  := as.character(A2)]
gwas_map[, FRQ_A_gwas := suppressWarnings(as.numeric(FRQ_A_gwas))]
gwas_map[, FRQ_U_gwas := suppressWarnings(as.numeric(FRQ_U_gwas))]
gwas_map[, OR  := suppressWarnings(as.numeric(OR))]
gwas_map[, SE  := suppressWarnings(as.numeric(SE))]
gwas_map[, P   := suppressWarnings(as.numeric(P))]
gwas_map[, INFO := suppressWarnings(as.numeric(INFO))]
gwas_map <- unique(gwas_map[!is.na(SNP) & SNP != ""])

# =========================================================
# 8) GENE-LEVEL INTEGRATION
# =========================================================
gene_panel <- copy(g_gene)
drop_if_exists(gene_panel, c("PP_H4", "coloc_tier", "coloc_bestTissue", "coloc_leadSNP"))
gene_panel <- merge(gene_panel, coloc_gene, by = "Gene", all.x = TRUE)
gene_panel <- merge(gene_panel, smr_best,   by = "Gene", all.x = TRUE)

gene_panel[, GTEx_L1      := !is.na(bestLevel_GTEx) & bestLevel_GTEx == "Strong"]
gene_panel[, GTEx_L2      := !is.na(bestLevel_GTEx) & bestLevel_GTEx == "Moderate"]
gene_panel[, coloc_STRONG := !is.na(coloc_tier) & coloc_tier == "STRONG"]
gene_panel[, coloc_MID    := !is.na(coloc_tier) & coloc_tier == "MID"]
gene_panel[, SMR_support  := !is.na(SMR_HEIDI_support) & SMR_HEIDI_support == TRUE]

gene_panel[, gene_strong_score := rowSums(cbind(
  as.integer(GTEx_L1),
  as.integer(coloc_STRONG),
  as.integer(SMR_support)
), na.rm = TRUE)]

gene_panel[, Gene_Evidence_Level := "Unclassified"]
gene_panel[gene_strong_score >= 2, Gene_Evidence_Level := "HighConfidence"]
gene_panel[Gene_Evidence_Level == "Unclassified" & gene_strong_score == 1, Gene_Evidence_Level := "ModerateConfidence"]
gene_panel[Gene_Evidence_Level == "Unclassified" & (coloc_MID | GTEx_L2), Gene_Evidence_Level := "Suggestive"]

gene_panel[, Gene_Evidence_Reason := NA_character_]
gene_panel[Gene_Evidence_Level %in% c("HighConfidence", "ModerateConfidence"),
           Gene_Evidence_Reason := {
             v <- character(0)
             if (GTEx_L1)      v <- c(v, "GTEx_strong")
             if (coloc_STRONG) v <- c(v, "coloc_STRONG")
             if (SMR_support)  v <- c(v, "SMR_support")
             if (length(v) == 0) NA_character_ else paste(v, collapse = "+")
           }, by = Gene]

gene_panel[Gene_Evidence_Level == "Suggestive",
           Gene_Evidence_Reason := {
             v <- character(0)
             if (coloc_MID) v <- c(v, "coloc_MID")
             if (GTEx_L2)   v <- c(v, "GTEx_moderate")
             if (length(v) == 0) NA_character_ else paste(v, collapse = "+")
           }, by = Gene]

gene_panel[, GTEx_Evidence_Level := bestLevel_GTEx]
gene_panel[, GTEx_Best_Evidence_Level := bestLevel_GTEx]
gene_panel[, Coloc_Evidence_Level := map_coloc_level(coloc_tier)]
gene_panel[, GTEx_Best_Tissue_Abbrev := abbr_tissue(bestTissue)]
gene_panel[, Coloc_Best_Tissue_Abbrev := abbr_tissue(coloc_bestTissue)]

# =========================================================
# 9) SNP-LEVEL INTEGRATION
# =========================================================
snp_panel <- copy(g_snp)
snp_panel[, Gene := fifelse(!is.na(Gene_GTEx) & Gene_GTEx != "", as.character(Gene_GTEx), NA_character_)]
if ("CHR" %in% names(snp_panel)) snp_panel[, CHR := as.integer(CHR)]
if ("BP"  %in% names(snp_panel)) snp_panel[, BP  := as.integer(BP)]

snp_panel <- merge(
  snp_panel,
  gwas_map[, .(
    SNP,
    CHR_gwas = CHR,
    BP_gwas = BP,
    A1_gwas = A1,
    A2_gwas = A2,
    FRQ_A_gwas,
    FRQ_U_gwas,
    INFO_gwas = INFO,
    OR_gwas   = OR,
    SE_gwas   = SE,
    P_gwas    = P
  )],
  by = "SNP",
  all.x = TRUE
)

if (!("A1" %in% names(snp_panel))) snp_panel[, A1 := NA_character_]
if (!("A2" %in% names(snp_panel))) snp_panel[, A2 := NA_character_]
if (!("CHR" %in% names(snp_panel))) snp_panel[, CHR := NA_integer_]
if (!("BP"  %in% names(snp_panel))) snp_panel[, BP  := NA_integer_]

snp_panel[is.na(A1) | A1 == "", A1 := A1_gwas]
snp_panel[is.na(A2) | A2 == "", A2 := A2_gwas]
snp_panel[is.na(CHR), CHR := CHR_gwas]
snp_panel[is.na(BP),  BP  := BP_gwas]

snp_panel[, EA  := A1]
snp_panel[, NEA := A2]
if (!("EAF" %in% names(snp_panel))) snp_panel[, EAF := FRQ_A_gwas]

snp_panel[, SMR_topSNP_any  := FALSE]
snp_panel[, SMR_topSNP_pass := FALSE]
if ("SMR_topSNP" %in% names(smr_best)) {
  snp_panel[SNP %in% smr_best$SMR_topSNP, SMR_topSNP_any := TRUE]
  snp_panel[SNP %in% smr_best[SMR_HEIDI_support == TRUE, SMR_topSNP], SMR_topSNP_pass := TRUE]
}

snp_panel[, coloc_locus_hit  := FALSE]
snp_panel[, coloc_locus_tier := NA_character_]

if ("coloc_leadSNP" %in% names(coloc_gene) && all(c("CHR", "BP") %in% names(snp_panel))) {
  lead_map <- unique(
    coloc_gene[
      !is.na(coloc_leadSNP) & coloc_leadSNP != "" & coloc_tier %in% c("STRONG", "MID"),
      .(coloc_leadSNP, coloc_tier)
    ]
  )
  lead_map <- merge(
    lead_map,
    gwas_map[, .(SNP, CHR, BP)],
    by.x = "coloc_leadSNP",
    by.y = "SNP",
    all.x = TRUE
  )
  lead_map <- lead_map[!is.na(CHR) & !is.na(BP)]

  if (nrow(lead_map) > 0) {
    for (ii in seq_len(nrow(lead_map))) {
      chr0  <- lead_map$CHR[ii]
      bp0   <- lead_map$BP[ii]
      tier0 <- as.character(lead_map$coloc_tier[ii])

      lo <- bp0 - WINDOW
      hi <- bp0 + WINDOW
      hits_idx <- which(snp_panel$CHR == chr0 & snp_panel$BP >= lo & snp_panel$BP <= hi)

      if (length(hits_idx) > 0) {
        snp_panel[hits_idx, coloc_locus_hit := TRUE]
        snp_panel[hits_idx, coloc_locus_tier := {
          cur <- coloc_locus_tier
          ifelse(rank_coloc_tier(tier0) > rank_coloc_tier(cur), tier0, cur)
        }]
      }
    }
  }
}

snp_panel <- merge(
  snp_panel,
  gene_panel[, .(Gene, Gene_Evidence_Level, gene_strong_score, Gene_Evidence_Reason)],
  by = "Gene",
  all.x = TRUE
)

snp_panel[, snp_regstrong   := !is.na(GTEx_evidence_class) & GTEx_evidence_class == "Strong"]
snp_panel[, snp_regmoderate := !is.na(GTEx_evidence_class) & GTEx_evidence_class == "Moderate"]
snp_panel[, snp_coloc       := !is.na(coloc_locus_hit) & coloc_locus_hit == TRUE]
snp_panel[, snp_coloc_STRONG := !is.na(coloc_locus_tier) & coloc_locus_tier == "STRONG"]
snp_panel[, snp_coloc_MID    := !is.na(coloc_locus_tier) & coloc_locus_tier == "MID"]
snp_panel[, snp_smr         := !is.na(SMR_topSNP_pass) & SMR_topSNP_pass == TRUE]

snp_panel[, multi_support_count := rowSums(cbind(
  as.integer(snp_regstrong),
  as.integer(snp_coloc),
  as.integer(snp_smr)
), na.rm = TRUE)]

snp_panel[, strong_score_snp := multi_support_count]

snp_panel[, SNP_Evidence_Level := "Unclassified"]
snp_panel[multi_support_count >= 2, SNP_Evidence_Level := "MultiEvidence"]
snp_panel[SNP_Evidence_Level == "Unclassified" & snp_coloc,       SNP_Evidence_Level := "Colocalized"]
snp_panel[SNP_Evidence_Level == "Unclassified" & snp_smr,         SNP_Evidence_Level := "SMRSupported"]
snp_panel[SNP_Evidence_Level == "Unclassified" & snp_regstrong,   SNP_Evidence_Level := "RegulatoryStrong"]
snp_panel[SNP_Evidence_Level == "Unclassified" & snp_regmoderate, SNP_Evidence_Level := "RegulatoryModerate"]

snp_panel[, SNP_Evidence_Reason := {
  v <- character(0)
  if (snp_regstrong) v <- c(v, "GTEx_strong")
  if (snp_regmoderate && !snp_regstrong) v <- c(v, "GTEx_moderate")
  if (snp_coloc_STRONG) {
    v <- c(v, "coloc_STRONG")
  } else if (snp_coloc_MID) {
    v <- c(v, "coloc_MID")
  } else if (snp_coloc) {
    v <- c(v, "coloc")
  }
  if (snp_smr) v <- c(v, "SMR_support")
  if (length(v) == 0) NA_character_ else paste(v, collapse = "+")
}, by = seq_len(nrow(snp_panel))]

snp_panel[, GTEx_Evidence_Level := GTEx_evidence_class]
snp_panel[, Coloc_Evidence_Level := map_coloc_level(coloc_locus_tier)]
snp_panel[, Tissue_Abbrev := abbr_tissue(Tissue)]

# =========================================================
# 10) STRICT SNP FILTER
# =========================================================
snp_panel[, strict_snp_pass := is_strict_snp(A1, A2)]

strict_fail_dt <- snp_panel[strict_snp_pass == FALSE, .(
  SNP, SNP_type, CHR, BP, A1, A2, EA, NEA,
  Tissue, Tissue_Abbrev, Gene, Gene_GTEx,
  SNP_Evidence_Level, SNP_Evidence_Reason
)]

snp_panel_strict <- copy(snp_panel[strict_snp_pass == TRUE])
gene_panel_strict <- copy(gene_panel)

# =========================================================
# 11) QC SUMMARY
# =========================================================
diag <- data.table(
  n_gene_total_before_strict_filter = nrow(gene_panel),
  n_gene_total_after_strict_filter = nrow(gene_panel_strict),
  n_gene_HighConfidence = sum(gene_panel_strict$Gene_Evidence_Level == "HighConfidence", na.rm = TRUE),
  n_gene_ModerateConfidence = sum(gene_panel_strict$Gene_Evidence_Level == "ModerateConfidence", na.rm = TRUE),
  n_gene_Suggestive = sum(gene_panel_strict$Gene_Evidence_Level == "Suggestive", na.rm = TRUE),
  n_gene_Unclassified = sum(gene_panel_strict$Gene_Evidence_Level == "Unclassified", na.rm = TRUE),

  n_snp_total_before_strict_filter = nrow(snp_panel),
  n_snp_removed_by_strict_filter = sum(snp_panel$strict_snp_pass == FALSE, na.rm = TRUE),
  n_snp_total_after_strict_filter = nrow(snp_panel_strict),

  n_snp_MultiEvidence = sum(snp_panel_strict$SNP_Evidence_Level == "MultiEvidence", na.rm = TRUE),
  n_snp_SMRSupported = sum(snp_panel_strict$SNP_Evidence_Level == "SMRSupported", na.rm = TRUE),
  n_snp_Colocalized = sum(snp_panel_strict$SNP_Evidence_Level == "Colocalized", na.rm = TRUE),
  n_snp_RegulatoryStrong = sum(snp_panel_strict$SNP_Evidence_Level == "RegulatoryStrong", na.rm = TRUE),
  n_snp_RegulatoryModerate = sum(snp_panel_strict$SNP_Evidence_Level == "RegulatoryModerate", na.rm = TRUE),
  n_snp_Unclassified = sum(snp_panel_strict$SNP_Evidence_Level == "Unclassified", na.rm = TRUE)
)

gene_cnt <- gene_panel_strict[, .N, by = .(Gene_Evidence_Level)][order(Gene_Evidence_Level)]
snp_cnt  <- snp_panel_strict[, .N, by = .(SNP_Evidence_Level)][order(SNP_Evidence_Level)]
strict_fail_cnt <- snp_panel[strict_snp_pass == FALSE, .N, by = .(SNP_Evidence_Level)][order(-N)]

# =========================================================
# 12) AUXILIARY TABLES
# =========================================================
debug_overlap_genes <- gene_panel_strict[
  (coloc_tier %in% c("STRONG", "MID")) | (SMR_support == TRUE) | (GTEx_L1 == TRUE) | (GTEx_L2 == TRUE),
  .(
    Gene,
    Gene_Evidence_Level,
    gene_strong_score,
    Gene_Evidence_Reason,
    GTEx_Evidence_Level,
    GTEx_Best_Evidence_Level,
    bestTissue,
    GTEx_Best_Tissue_Abbrev,
    bestSNP,
    bestP_GTEx,
    Coloc_Evidence_Level,
    coloc_bestTissue,
    Coloc_Best_Tissue_Abbrev,
    coloc_leadSNP,
    PP_H4,
    p_SMR,
    p_HEIDI,
    SMR_support
  )
]

# =========================================================
# 13) FINAL OUTPUT TABLES
# =========================================================
final_gene_list <- gene_panel_strict[
  Gene_Evidence_Level %in% c("HighConfidence", "ModerateConfidence", "Suggestive"),
  .(
    Gene,
    Gene_Evidence_Level,
    Gene_Evidence_Reason,
    gene_strong_score,

    GTEx_Evidence_Level,
    GTEx_Best_Evidence_Level,
    GTEx_Best_P = bestP_GTEx,
    GTEx_Best_SNP = bestSNP,
    GTEx_Best_Tissue = bestTissue,
    GTEx_Best_Tissue_Abbrev,

    Coloc_Evidence_Level,
    PP_H4,
    Coloc_Best_Tissue = coloc_bestTissue,
    Coloc_Best_Tissue_Abbrev,
    Coloc_Lead_SNP = coloc_leadSNP,

    P_SMR = p_SMR,
    P_HEIDI = p_HEIDI,
    SMR_support,
    SMR_Top_SNP = SMR_topSNP,

    GTEx_Strong = GTEx_L1,
    GTEx_Moderate = GTEx_L2,
    Coloc_Strong = coloc_STRONG,
    Coloc_Moderate = coloc_MID
  )
]

final_snp_list <- snp_panel_strict[
  SNP_Evidence_Level %in% c("MultiEvidence", "SMRSupported", "Colocalized", "RegulatoryStrong", "RegulatoryModerate"),
  .(
    SNP,
    SNP_type,
    CHR,
    BP,
    A1,
    A2,
    EA,
    NEA,
    EAF,

    FRQ_A_gwas,
    FRQ_U_gwas,
    OR_gwas,
    SE_gwas,
    P_GWAS = P_gwas,
    INFO_gwas,

    Tissue,
    Tissue_Abbrev,
    Gene_GTEx,
    P_GTEx,
    GTEx_Evidence_Level,

    SMR_Top_SNP_Pass = SMR_topSNP_pass,
    SMR_support = snp_smr,
    SMR_Top_SNP_Any = SMR_topSNP_any,

    Colocalized = snp_coloc,
    Coloc_Evidence_Level,
    strong_score_snp,

    SNP_Evidence_Level,
    SNP_Evidence_Reason,

    Gene,
    Gene_Evidence_Level,
    Gene_Evidence_Reason
  )
]

top_multi_evidence_snps <- snp_panel_strict[
  SNP_Evidence_Level == "MultiEvidence",
  .(
    SNP, CHR, BP, A1, A2, EA, NEA, EAF,
    Tissue, Tissue_Abbrev,
    GTEx_Evidence_Level,
    Coloc_Evidence_Level,
    SMR_Top_SNP_Pass = SMR_topSNP_pass,
    SMR_support = snp_smr,
    strong_score_snp,
    SNP_Evidence_Level,
    SNP_Evidence_Reason,
    Gene,
    Gene_Evidence_Level,
    Gene_Evidence_Reason,
    FRQ_A_gwas,
    FRQ_U_gwas,
    OR_gwas,
    P_GWAS = P_gwas
  )
]

snp_panel_gwas <- snp_panel_strict[SNP %in% gwas_map$SNP]
snp_panel_gwas[, locus_bin := as.integer(BP / WINDOW_BIN)]
snp_panel_gwas[, locus_key := paste0(CHR, "_", locus_bin)]
snp_panel_gwas[, is_lead := as.integer(SNP_type == "LEAD")]
snp_panel_gwas[, P_gwas_num := suppressWarnings(as.numeric(P_gwas))]
snp_panel_gwas[, P_GTEx_num := suppressWarnings(as.numeric(P_GTEx))]
setorder(snp_panel_gwas, -is_lead, -strong_score_snp, P_gwas_num, P_GTEx_num)
snp_gwas_top <- snp_panel_gwas[, .SD[1], by = locus_key]

# =========================================================
# 14) SAVE INTEGRATION OUTPUTS
# =========================================================
fwrite(gene_panel,        file.path(dir_tables, "Step6_gene_panel_unified_raw.tsv"), sep = "\t")
fwrite(snp_panel,         file.path(dir_tables, "Step6_snp_panel_unified_raw.tsv"), sep = "\t")
fwrite(gene_panel_strict, file.path(dir_tables, "Step6_gene_panel_unified.tsv"), sep = "\t")
fwrite(snp_panel_strict,  file.path(dir_tables, "Step6_snp_panel_unified.tsv"), sep = "\t")

fwrite(final_gene_list,         file.path(dir_tables, "Step6_final_candidate_genes.tsv"), sep = "\t")
fwrite(final_snp_list,          file.path(dir_tables, "Step6_final_candidate_snps.tsv"), sep = "\t")
fwrite(top_multi_evidence_snps, file.path(dir_tables, "Step6_multi_evidence_snps.tsv"), sep = "\t")
fwrite(snp_gwas_top,            file.path(dir_tables, "Step6_gwas_matched_topSNP_per_500kb_bin.tsv"), sep = "\t")

fwrite(gene_cnt,        file.path(dir_summary, "Step6_counts_gene_panel.tsv"), sep = "\t")
fwrite(snp_cnt,         file.path(dir_summary, "Step6_counts_snp_panel.tsv"), sep = "\t")
fwrite(diag,            file.path(dir_summary, "Step6_integration_diagnostics.tsv"), sep = "\t")
fwrite(strict_fail_cnt, file.path(dir_summary, "Step6_strict_filter_removed_snp_counts.tsv"), sep = "\t")

fwrite(debug_overlap_genes, file.path(dir_debug, "Step6_debug_overlap_genes.tsv"), sep = "\t")
fwrite(strict_fail_dt,      file.path(dir_debug, "Step6_debug_removed_nonstandard_snps.tsv"), sep = "\t")

# =========================================================
# 15) STEP 6 MAIN TABLES
# =========================================================
step6_table_highconfidence_genes <- copy(gene_panel_strict[Gene_Evidence_Level == "HighConfidence"])

step6_table_highconfidence_genes_out <- step6_table_highconfidence_genes[, .(
  Gene,
  Gene_Evidence_Level = as.character(Gene_Evidence_Level),
  Gene_Evidence_Reason,

  GTEx_Evidence_Level,
  GTEx_Best_Tissue = GTEx_Best_Tissue_Abbrev,
  GTEx_Best_SNP = bestSNP,
  GTEx_Best_P = bestP_GTEx,

  Coloc_Evidence_Level,
  PP_H4 = ifelse(is.na(PP_H4), NA, signif(PP_H4, 3)),
  Coloc_Best_Tissue = Coloc_Best_Tissue_Abbrev,
  Coloc_Lead_SNP = coloc_leadSNP,

  P_SMR = fmt_p(p_SMR),
  P_HEIDI = fmt_p(p_HEIDI),
  SMR_support,
  SMR_Top_SNP = SMR_topSNP
)]
setorder(step6_table_highconfidence_genes_out, Gene)
fwrite(step6_table_highconfidence_genes_out, file.path(dir_tables, "Step6_highconfidence_genes_summary.tsv"), sep = "\t")

cand_snp <- copy(snp_panel_strict[
  SNP_Evidence_Level %in% c("MultiEvidence", "RegulatoryStrong", "RegulatoryModerate", "Colocalized", "SMRSupported")
])

cand_snp[, snp_class_rank := fifelse(SNP_Evidence_Level == "MultiEvidence", 1L,
                                     fifelse(SNP_Evidence_Level == "RegulatoryStrong", 2L,
                                             fifelse(SNP_Evidence_Level == "RegulatoryModerate", 3L,
                                                     fifelse(SNP_Evidence_Level == "Colocalized", 4L,
                                                             fifelse(SNP_Evidence_Level == "SMRSupported", 5L, 9L)))))]

if ("P_GTEx" %in% names(cand_snp)) cand_snp[, P_GTEx_num := suppressWarnings(as.numeric(P_GTEx))] else cand_snp[, P_GTEx_num := NA_real_]
if ("P_gwas" %in% names(cand_snp)) cand_snp[, P_GWAS_num := suppressWarnings(as.numeric(P_gwas))] else cand_snp[, P_GWAS_num := NA_real_]

setorderv(cand_snp, cols = c("snp_class_rank", "P_GWAS_num", "P_GTEx_num"), order = c(1, 1, 1))
step6_table_rep_snp <- cand_snp[, .SD[1], by = Gene]
step6_table_strong_snp <- cand_snp[SNP_Evidence_Level == "MultiEvidence"]
step6_table_representative_snp <- unique(rbindlist(list(step6_table_strong_snp, step6_table_rep_snp), fill = TRUE), by = "SNP")

step6_table_representative_snp_out <- step6_table_representative_snp[, .(
  SNP,
  SNP_Evidence_Level = as.character(SNP_Evidence_Level),
  SNP_Evidence_Reason,
  Gene,
  Gene_Evidence_Level = as.character(Gene_Evidence_Level),
  Tissue = Tissue_Abbrev,
  P_GTEx = fmt_p(P_GTEx_num),
  P_GWAS = fmt_p(P_GWAS_num),
  Coloc_Evidence_Level,
  SMR_support = snp_smr,
  A1,
  A2
)]
setorder(step6_table_representative_snp_out, Gene_Evidence_Level, SNP_Evidence_Level, Gene, SNP)
fwrite(step6_table_representative_snp_out, file.path(dir_tables, "Step6_representative_candidate_snps_summary.tsv"), sep = "\t")

# =========================================================
# 16) SUPPLEMENTARY TABLES
# =========================================================
supp_gene_cols <- intersect(c(
  "Gene", "Gene_Evidence_Level", "gene_strong_score", "Gene_Evidence_Reason",
  "GTEx_Evidence_Level", "GTEx_Best_Evidence_Level", "bestP_GTEx", "bestSNP",
  "GTEx_Best_Tissue_Abbrev", "PP_H4", "Coloc_Evidence_Level", "Coloc_Best_Tissue_Abbrev",
  "coloc_leadSNP", "p_SMR", "p_HEIDI", "SMR_support", "SMR_topSNP",
  "GTEx_L1", "GTEx_L2", "coloc_STRONG", "coloc_MID"
), names(gene_panel_strict))

supp_snp_cols <- intersect(c(
  "SNP", "SNP_type", "CHR", "BP", "A1", "A2", "EA", "NEA",
  "EAF", "FRQ_A_gwas", "FRQ_U_gwas", "OR_gwas", "SE_gwas", "P_gwas",
  "INFO_gwas", "Tissue_Abbrev", "Gene_GTEx", "P_GTEx", "GTEx_Evidence_Level",
  "SMR_topSNP_pass", "snp_smr", "SMR_topSNP_any", "snp_coloc",
  "Coloc_Evidence_Level", "strong_score_snp", "SNP_Evidence_Level", "SNP_Evidence_Reason",
  "Gene", "Gene_Evidence_Level", "Gene_Evidence_Reason"
), names(snp_panel_strict))

fwrite(
  gene_panel_strict[Gene_Evidence_Level %in% c("HighConfidence", "ModerateConfidence", "Suggestive"), ..supp_gene_cols],
  file.path(dir_tables, "Step6_supplementary_genes_prioritized.tsv"),
  sep = "\t"
)

fwrite(
  snp_panel_strict[SNP_Evidence_Level %in% c("MultiEvidence", "RegulatoryStrong", "RegulatoryModerate", "Colocalized", "SMRSupported"), ..supp_snp_cols],
  file.path(dir_tables, "Step6_supplementary_candidate_snps.tsv"),
  sep = "\t"
)

# =========================================================
# 17) STEP 6 FIGURE A
# =========================================================
count_gene <- gene_panel_strict[, .N, by = Gene_Evidence_Level]
count_gene[, Gene_Evidence_Level := factor(
  as.character(Gene_Evidence_Level),
  levels = c("HighConfidence", "ModerateConfidence", "Suggestive", "Unclassified")
)]
count_gene <- count_gene[order(Gene_Evidence_Level)]

pA <- ggplot(count_gene, aes(x = Gene_Evidence_Level, y = N, fill = Gene_Evidence_Level)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.3) +
  geom_text(aes(label = N), vjust = -0.35, size = 4.2) +
  scale_fill_manual(values = c(
    "HighConfidence" = "#D55E00",
    "ModerateConfidence" = "#E69F00",
    "Suggestive" = "#56B4E9",
    "Unclassified" = "grey75"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "A  Integrated gene evidence levels", x = NULL, y = "Number of genes") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 15),
    axis.text.x = element_text(size = 11, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.position = "none"
  )

# =========================================================
# 18) STEP 6 FIGURE B
# =========================================================
plot_gene <- copy(gene_panel_strict[Gene_Evidence_Level %in% c("HighConfidence", "ModerateConfidence", "Suggestive")])
setorderv(plot_gene, cols = c("Gene_Evidence_Level", "gene_strong_score", "Gene"), order = c(1, -1, 1))
plot_gene[, Gene_plot := factor(Gene, levels = rev(unique(Gene)))]

ev_long <- rbindlist(list(
  plot_gene[, .(Gene, Gene_plot, Gene_Evidence_Level, Evidence = "GTEx strong",    Present = as.integer(GTEx_L1))],
  plot_gene[, .(Gene, Gene_plot, Gene_Evidence_Level, Evidence = "Coloc strong",   Present = as.integer(coloc_STRONG))],
  plot_gene[, .(Gene, Gene_plot, Gene_Evidence_Level, Evidence = "SMR support",    Present = as.integer(SMR_support))],
  plot_gene[, .(Gene, Gene_plot, Gene_Evidence_Level, Evidence = "GTEx moderate",  Present = as.integer(GTEx_L2))],
  plot_gene[, .(Gene, Gene_plot, Gene_Evidence_Level, Evidence = "Coloc moderate", Present = as.integer(coloc_MID))]
), fill = TRUE)

ev_long[, Evidence := factor(
  Evidence,
  levels = c("GTEx strong", "Coloc strong", "SMR support", "GTEx moderate", "Coloc moderate")
)]
ev_long[, Present_lab := ifelse(Present == 1, "Present", "Absent")]

pB <- ggplot(ev_long, aes(x = Evidence, y = Gene_plot, fill = Present_lab)) +
  geom_tile(color = "white", linewidth = 0.35) +
  facet_grid(Gene_Evidence_Level ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Present" = "#009E73", "Absent" = "grey92")) +
  labs(title = "B  Evidence matrix of prioritized genes", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 15),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9),
    strip.text.y = element_text(face = "bold", size = 10),
    legend.title = element_blank(),
    legend.position = "top"
  )

# =========================================================
# 19) STEP 6 FIGURE C
# =========================================================
count_snp <- snp_panel_strict[, .N, by = SNP_Evidence_Level]
count_snp[, SNP_Evidence_Level := factor(
  as.character(SNP_Evidence_Level),
  levels = c("MultiEvidence", "RegulatoryStrong", "RegulatoryModerate", "Colocalized", "SMRSupported", "Unclassified")
)]
count_snp <- count_snp[order(SNP_Evidence_Level)]

pC <- ggplot(count_snp, aes(x = SNP_Evidence_Level, y = N, fill = SNP_Evidence_Level)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.3) +
  geom_text(aes(label = N), vjust = -0.35, size = 3.8) +
  scale_fill_manual(values = c(
    "MultiEvidence" = "#D55E00",
    "RegulatoryStrong" = "#E69F00",
    "RegulatoryModerate" = "#56B4E9",
    "Colocalized" = "#009E73",
    "SMRSupported" = "#CC79A7",
    "Unclassified" = "grey75"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "C  SNP evidence categories", x = NULL, y = "Number of SNPs") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 15),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )

# =========================================================
# 20) SAVE FIGURES
# =========================================================
top_row <- pA | pC
bottom_row <- pB
p_combined <- top_row / bottom_row + plot_layout(heights = c(1.1, 3.2))

ggsave(file.path(dir_figures, "Step6_gene_evidence_levels.png"), pA, width = 6.6, height = 4.8, dpi = 600, bg = "white")
ggsave(file.path(dir_figures, "Step6_gene_evidence_levels.pdf"), pA, width = 6.6, height = 4.8, bg = "white")
ggsave(file.path(dir_figures, "Step6_gene_evidence_matrix.png"), pB, width = 9.0, height = 9.4, dpi = 600, bg = "white")
ggsave(file.path(dir_figures, "Step6_gene_evidence_matrix.pdf"), pB, width = 9.0, height = 9.4, bg = "white")
ggsave(file.path(dir_figures, "Step6_snp_evidence_distribution.png"), pC, width = 8.0, height = 4.9, dpi = 600, bg = "white")
ggsave(file.path(dir_figures, "Step6_snp_evidence_distribution.pdf"), pC, width = 8.0, height = 4.9, bg = "white")
ggsave(file.path(dir_figures, "Step6_combined_integrated_evidence.png"), p_combined, width = 11.8, height = 12.8, dpi = 600, bg = "white")
ggsave(file.path(dir_figures, "Step6_combined_integrated_evidence.pdf"), p_combined, width = 11.8, height = 12.8, bg = "white")

# =========================================================
# 21) SUMMARY TABLES
# =========================================================
summary_gene <- data.table(
  Metric = c("HighConfidence genes", "ModerateConfidence genes", "Suggestive genes", "Unclassified genes"),
  N = c(
    sum(gene_panel_strict$Gene_Evidence_Level == "HighConfidence", na.rm = TRUE),
    sum(gene_panel_strict$Gene_Evidence_Level == "ModerateConfidence", na.rm = TRUE),
    sum(gene_panel_strict$Gene_Evidence_Level == "Suggestive", na.rm = TRUE),
    sum(gene_panel_strict$Gene_Evidence_Level == "Unclassified", na.rm = TRUE)
  )
)

summary_snp <- data.table(
  Metric = c("MultiEvidence SNPs", "RegulatoryStrong SNPs", "RegulatoryModerate SNPs", "Colocalized SNPs", "SMRSupported SNPs", "Unclassified SNPs"),
  N = c(
    sum(snp_panel_strict$SNP_Evidence_Level == "MultiEvidence", na.rm = TRUE),
    sum(snp_panel_strict$SNP_Evidence_Level == "RegulatoryStrong", na.rm = TRUE),
    sum(snp_panel_strict$SNP_Evidence_Level == "RegulatoryModerate", na.rm = TRUE),
    sum(snp_panel_strict$SNP_Evidence_Level == "Colocalized", na.rm = TRUE),
    sum(snp_panel_strict$SNP_Evidence_Level == "SMRSupported", na.rm = TRUE),
    sum(snp_panel_strict$SNP_Evidence_Level == "Unclassified", na.rm = TRUE)
  )
)

fwrite(summary_gene, file.path(dir_summary, "Step6_summary_gene_levels.tsv"), sep = "\t")
fwrite(summary_snp,  file.path(dir_summary, "Step6_summary_snp_levels.tsv"), sep = "\t")

# =========================================================
# 22) MANIFEST
# =========================================================
manifest <- data.table(
  category = c(rep("tables", 8), rep("figures", 8), rep("summary", 6), rep("debug", 2)),
  file = c(
    "Step6_gene_panel_unified_raw.tsv",
    "Step6_snp_panel_unified_raw.tsv",
    "Step6_gene_panel_unified.tsv",
    "Step6_snp_panel_unified.tsv",
    "Step6_final_candidate_genes.tsv",
    "Step6_final_candidate_snps.tsv",
    "Step6_highconfidence_genes_summary.tsv",
    "Step6_representative_candidate_snps_summary.tsv",
    "Step6_gene_evidence_levels.png",
    "Step6_gene_evidence_levels.pdf",
    "Step6_gene_evidence_matrix.png",
    "Step6_gene_evidence_matrix.pdf",
    "Step6_snp_evidence_distribution.png",
    "Step6_snp_evidence_distribution.pdf",
    "Step6_combined_integrated_evidence.png",
    "Step6_combined_integrated_evidence.pdf",
    "Step6_counts_gene_panel.tsv",
    "Step6_counts_snp_panel.tsv",
    "Step6_integration_diagnostics.tsv",
    "Step6_summary_gene_levels.tsv",
    "Step6_summary_snp_levels.tsv",
    "Step6_strict_filter_removed_snp_counts.tsv",
    "Step6_debug_overlap_genes.tsv",
    "Step6_debug_removed_nonstandard_snps.tsv"
  ),
  full_path = c(
    file.path(dir_tables, "Step6_gene_panel_unified_raw.tsv"),
    file.path(dir_tables, "Step6_snp_panel_unified_raw.tsv"),
    file.path(dir_tables, "Step6_gene_panel_unified.tsv"),
    file.path(dir_tables, "Step6_snp_panel_unified.tsv"),
    file.path(dir_tables, "Step6_final_candidate_genes.tsv"),
    file.path(dir_tables, "Step6_final_candidate_snps.tsv"),
    file.path(dir_tables, "Step6_highconfidence_genes_summary.tsv"),
    file.path(dir_tables, "Step6_representative_candidate_snps_summary.tsv"),
    file.path(dir_figures, "Step6_gene_evidence_levels.png"),
    file.path(dir_figures, "Step6_gene_evidence_levels.pdf"),
    file.path(dir_figures, "Step6_gene_evidence_matrix.png"),
    file.path(dir_figures, "Step6_gene_evidence_matrix.pdf"),
    file.path(dir_figures, "Step6_snp_evidence_distribution.png"),
    file.path(dir_figures, "Step6_snp_evidence_distribution.pdf"),
    file.path(dir_figures, "Step6_combined_integrated_evidence.png"),
    file.path(dir_figures, "Step6_combined_integrated_evidence.pdf"),
    file.path(dir_summary, "Step6_counts_gene_panel.tsv"),
    file.path(dir_summary, "Step6_counts_snp_panel.tsv"),
    file.path(dir_summary, "Step6_integration_diagnostics.tsv"),
    file.path(dir_summary, "Step6_summary_gene_levels.tsv"),
    file.path(dir_summary, "Step6_summary_snp_levels.tsv"),
    file.path(dir_summary, "Step6_strict_filter_removed_snp_counts.tsv"),
    file.path(dir_debug, "Step6_debug_overlap_genes.tsv"),
    file.path(dir_debug, "Step6_debug_removed_nonstandard_snps.tsv")
  )
)
fwrite(manifest, file.path(out_root, "Step6_integration_manifest.tsv"), sep = "\t")

message("06_integration_and_visualization.R finished successfully.")
