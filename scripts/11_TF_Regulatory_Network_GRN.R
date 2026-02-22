# ==============================================================================
# 11_TF_Regulatory_Network_GRN.R
# ------------------------------------------------------------------------------
# TF regulon GRN visualization, module-specific deep dive,
# and ASCL1 focal network (Astro-M7 focus)
# ==============================================================================
# Project: Heat_Call_Single_Nucleus_Multiome
# Author:  Prakrit Subba
# Date:    2026
#
# Paper: Subba et al. 2026 — Multiome Profiling Reveals Astrocyte and
#        Neuroendocrine Targets of Prenatal Acoustic Programming in
#        Zebra Finch Embryos
#
# Description:
#   Extends 09_hdWGCNA_Coexpression_Networks.R with per-module TF regulon
#   visualization, ASCL1 focal GRN, and GO enrichment gene-list export.
#   Loads pre-computed regulon scores and differential regulon tables from
#   results/09_hdWGCNA/ rather than recomputing.
#
#   Steps:
#     (1)  Load hdWGCNA object + pre-computed regulon scores (Script 09)
#     (2)  Transfer subcluster metadata + wnn_astro embedding
#     (3)  Set TFgroup metadata (Heat / Control)
#     (4)  Load JASPAR TRUE-TF filter
#     (5)  Load differential regulon results (Script 09 manual method)
#     (6)  Get module membership
#     (7)  Global quadrant plot + UMAP figures
#     (8)  Module-specific analysis loop (Astro-M7):
#            regulon bar plots, TF network plots,
#            UMAP visualizations, loose GO gene lists
#     (9)  ASCL1 focal network (strategy C, depth 2, Cor edges)  → Fig 3D
#     (10) TOM-scaffold GRN (ggraph)
#
# Inputs:  (all from results/09_hdWGCNA/)
#   AstroTFNetworkCheckpoint.rds
#   subclusters.rds
#   AstroPositiveRegulonScoresAstroNovTruth.rds
#   AstroNegativeRegulonScoresAstroNovTruth.rds
#   AstroDiffRegulonsManualHeatVsControlAllTFsNovTruth.csv
#   data/jaspar_mapping.rds                (optional)
#
# Outputs:  (written to results/11_TF_GRN/)
#   0GlobalAnalysis/DiffRegulonsQuadrantGlobal[Labeled].pdf
#   0GlobalAnalysis/UMAPSubclusters.pdf
#   0GlobalAnalysis/UMAPTFgroup.pdf
#   Astro-M7Analysis/{1–6}/ per-TF figures and GO lists
#   ASCL1_TFNetwork_strategyCdepth2Cor_labeled.pdf   ← Figure 3D
#   ASCL1_TOMScaffoldNetwork.pdf                     ← Figure 3E
#
# Key naming conventions (EXACT — match Script 09):
#   wgcna_name        = "AstroAnnotation"
#   Module name       = "Astro-M7"     (hyphen)
#   TFgroup values    = "Heat" / "Control"   (NOT "Heat Call")
#   Conditionnum      :  0 = Heat Call  →  TFgroup = "Heat"
#                       >0 = Control    →  TFgroup = "Control"
#   Diff regulon cols : avglog2FCpositive / pvaladjpositive   (NO underscores)
#                       avglog2FCnegative / pvaladjnegative
#   Subcluster col    = "astrosubclusterwnn"
#   Reduction name    = "wnn_astro"
#
# Previous Script: 10_Cicero_Chromatin_Rewiring.R
# Next Script:     12_Pseudotime_Trajectory_Analysis.R
# ==============================================================================


# ── 0. LIBRARIES ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(hdWGCNA)
  library(WGCNA)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(dplyr)
  library(tibble)
  library(Matrix)
  library(ggrepel)
  library(igraph)
  library(ggraph)
})
set.seed(12345)
theme_set(theme_cowplot())


# ── 1. CONFIGURATION ──────────────────────────────────────────────────────────
# Input paths — all from Script 09
SCRIPT09_DIR      <- "results/09_hdWGCNA"
TF_CHECKPOINT_RDS <- "AstroTFNetworkCheckpoint.rds"
SUBCLUSTERS_RDS   <- "subclusters.rds"
POS_SCORES_RDS    <- file.path(SCRIPT09_DIR,
                                "AstroPositiveRegulonScoresAstroNovTruth.rds")
NEG_SCORES_RDS    <- file.path(SCRIPT09_DIR,
                                "AstroNegativeRegulonScoresAstroNovTruth.rds")
DREGS_CSV         <- file.path(SCRIPT09_DIR,
                                "AstroDiffRegulonsManualHeatVsControlAllTFsNovTruth.csv")
JASPAR_MAPPING_RDS <- "data/jaspar_mapping.rds"

# Output
OUT_DIR <- "results/11_TF_GRN"
DIR_FIG <- file.path(OUT_DIR, "FIGURES")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, showWarnings = FALSE)

# EXACT WGCNA slot name — must match Script 09
wgcna_name <- "AstroAnnotation"

# Modules to analyze in detail
MODULES_OF_INTEREST <- c("Astro-M7")

# Minimum targets for GO list export
GO_MIN_TARGETS <- 10

cat("Script: 11_TF_Regulatory_Network_GRN.R\n")
cat("OUT_DIR:", OUT_DIR, "\n")
cat("wgcna_name:", wgcna_name, "\n\n")


# ── 2. LOAD HDWGCNA OBJECT AND PRE-COMPUTED REGULON SCORES ───────────────────
cat("=== STEP 1: LOADING HDWGCNA OBJECT ===\n")

# Load the TF network checkpoint (combined_wnn with hdWGCNA slot populated)
combined_wnn <- readRDS(TF_CHECKPOINT_RDS)
combined_wnn@misc$active_wgcna <- wgcna_name
DefaultAssay(combined_wnn)     <- "RNA"

stopifnot(
  "cluster_annotation" %in% colnames(combined_wnn@meta.data),
  "Condition"          %in% colnames(combined_wnn@meta.data)
)

cat("Cells in combined_wnn:", ncol(combined_wnn), "\n")

# Load pre-computed regulon scores from Script 09
# (avoids re-running RegulonScores — computationally expensive)
cat("Loading pre-computed regulon scores from Script 09...\n")
stopifnot(file.exists(POS_SCORES_RDS),
          file.exists(NEG_SCORES_RDS))

pos_scores <- readRDS(POS_SCORES_RDS)
neg_scores <- readRDS(NEG_SCORES_RDS)

cat("Positive regulon scores dims:",
    paste(dim(pos_scores), collapse = " x "), "\n")
cat("Negative regulon scores dims:",
    paste(dim(neg_scores), collapse = " x "), "\n\n")

# Load differential regulon results from Script 09
cat("Loading differential regulon results from Script 09...\n")
stopifnot(file.exists(DREGS_CSV))
dregs_manual <- read.csv(DREGS_CSV, stringsAsFactors = FALSE)
cat("TFs in dregs_manual:", nrow(dregs_manual), "\n")

# Verify exact column names match Script 09 (NO underscores)
expected_cols <- c("tf", "avglog2FCpositive", "pvaladjpositive",
                   "avglog2FCnegative", "pvaladjnegative", "ntargets")
missing_cols  <- setdiff(expected_cols, colnames(dregs_manual))
if (length(missing_cols) > 0)
  warning("Missing expected columns in dregs_manual: ",
          paste(missing_cols, collapse = ", "))


# ── 3. LOAD ASTROCYTE SUBCLUSTER OBJECT ──────────────────────────────────────
cat("\n=== STEP 2: LOADING ASTROCYTE SUBCLUSTERS ===\n")

stopifnot(file.exists(SUBCLUSTERS_RDS))
astrosub <- readRDS(SUBCLUSTERS_RDS)

stopifnot(
  "wsnn_res.0.06"    %in% colnames(astrosub@meta.data),
  all(colnames(astrosub) %in% colnames(combined_wnn))
)

# Copy subcluster labels to combined_wnn
# Column name: "astrosubclusterwnn" (exact from Script 09 violin plot)
combined_wnn$astrosubclusterwnn <- NA_character_
combined_wnn$astrosubclusterwnn[colnames(astrosub)] <-
  as.character(astrosub$wsnn_res.0.06)

cat("Subclusters:", length(unique(astrosub$wsnn_res.0.06)), "\n")
print(table(combined_wnn$astrosubclusterwnn, useNA = "ifany"))

# Extract wnn_astro embedding (reduction name = "wnn_astro", exact from Script 09)
if ("wnn_astro" %in% Reductions(astrosub)) {
  emb <- Embeddings(astrosub, "wnn_astro")
  cat("wnn_astro embedding found:", nrow(emb), "cells\n\n")
} else {
  emb <- NULL
  cat("NOTE: wnn_astro reduction not found in astrosub.\n\n")
}


# ── 4. SET TFGROUP METADATA ───────────────────────────────────────────────────
cat("=== STEP 3: SETTING TFgroup METADATA ===\n")
# TFgroup values: "Heat" / "Control" (NOT "Heat Call")
# Conditionnum encoding (from Script 09):
#   0  = Heat Call  →  TFgroup = "Heat"
#  >0  = Control    →  TFgroup = "Control"

astro_cells <- colnames(combined_wnn)[combined_wnn$cluster_annotation == "Astrocyte"]
stopifnot(length(astro_cells) > 0)

combined_wnn$TFgroup <- NA_character_

if ("Conditionnum" %in% colnames(combined_wnn@meta.data)) {
  combined_wnn$TFgroup[astro_cells] <- ifelse(
    combined_wnn$Conditionnum[astro_cells] > 0,
    "Control",
    "Heat"
  )
} else if ("Condition" %in% colnames(combined_wnn@meta.data)) {
  combined_wnn$TFgroup[astro_cells] <-
    as.character(combined_wnn$Condition[astro_cells])
} else {
  stop("No Conditionnum or Condition column found in metadata.")
}

cat("TFgroup counts (full object):\n")
print(table(combined_wnn$TFgroup, useNA = "ifany"))
cat("\n")

# Build astrocyte-only object
astroobj <- subset(combined_wnn, cells = astro_cells)
DefaultAssay(astroobj)         <- "RNA"
astroobj@misc$active_wgcna    <- wgcna_name

# Copy WGCNA slot (using exact misc slot structure from Script 09)
if (!is.null(combined_wnn@misc[[wgcna_name]])) {
  astroobj@misc[[wgcna_name]] <- combined_wnn@misc[[wgcna_name]]
}

# Attach wnn_astro reduction to astroobj
if (!is.null(emb)) {
  emb2 <- emb[colnames(astroobj), , drop = FALSE]
  astroobj[["wnn_astro"]] <- CreateDimReducObject(
    embeddings = emb2,
    key        = "wnnastro_",
    assay      = DefaultAssay(astroobj)
  )
}

cat("Astrocyte cells:", ncol(astroobj), "\n")
cat("TFgroup in astroobj:\n")
print(table(astroobj$TFgroup, useNA = "ifany"))
cat("\n")


# ── 5. ENSURE TF REGULONS PRESENT ────────────────────────────────────────────
cat("=== STEP 4: CHECKING TF REGULONS ===\n")

tf_regulons_full <- tryCatch(
  GetTFRegulons(combined_wnn),
  error = function(e) NULL
)

if (is.null(tf_regulons_full) || nrow(tf_regulons_full) == 0) {
  cat("TF regulons missing — running AssignTFRegulons (strategy A)...\n")
  combined_wnn <- AssignTFRegulons(
    combined_wnn,
    strategy  = "A",
    regthresh = 0.01,
    ntfs      = 10
  )
  tf_regulons_full <- GetTFRegulons(combined_wnn)
  if (is.null(tf_regulons_full) || nrow(tf_regulons_full) == 0)
    stop("AssignTFRegulons failed — cannot retrieve TF regulons.")
  # Propagate to astroobj
  astroobj@misc[[wgcna_name]]$TFregulons <- tf_regulons_full
}

cat("TF regulon pairs:", nrow(tf_regulons_full), "\n\n")


# ── 6. LOAD JASPAR TRUE-TF FILTER ─────────────────────────────────────────────
cat("=== STEP 5: JASPAR TRUE-TF FILTER ===\n")

if (file.exists(JASPAR_MAPPING_RDS)) {
  tf_mapping       <- readRDS(JASPAR_MAPPING_RDS)
  true_tfs         <- unique(tf_mapping$tf_gene)   # column = tf_gene (Script 09)
  cat("TRUE TFs from JASPAR:", length(true_tfs), "\n\n")
} else {
  cat("NOTE: jaspar_mapping.rds not found — proceeding without TRUE TF filter.\n\n")
  true_tfs <- NULL
}


# ── 7. GET MODULE MEMBERSHIP ──────────────────────────────────────────────────
cat("=== STEP 6: MODULE MEMBERSHIP ===\n")

modules        <- GetModules(astroobj, wgcnaname = wgcna_name)
modules_simple <- unique(modules[, c("gene_name", "module")])
modules_simple <- as.data.frame(modules_simple)
cat("Module summary:\n")
print(table(modules$module))
cat("\n")


# ── 8. GLOBAL DIFFERENTIAL REGULON PLOTS ─────────────────────────────────────
cat("=== STEP 7: GLOBAL DIFFERENTIAL REGULON PLOTS ===\n")

global_dir <- file.path(OUT_DIR, "0GlobalAnalysis")
dir.create(global_dir, showWarnings = FALSE)

dregs_plot <- merge(dregs_manual, modules_simple,
                    by.x = "tf", by.y = "gene_name", all.x = TRUE)
dregs_plot <- as.data.frame(dregs_plot)

# Quadrant plot: positive vs negative regulon log2FC
p_dregs <- ggplot(dregs_plot,
                  aes(x = avglog2FCpositive,
                      y = avglog2FCnegative,
                      color = module)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_discrete(na.value = "grey70") +
  theme_cowplot() +
  labs(x     = "log2FC positive regulon (Heat vs Control)",
       y     = "log2FC negative regulon (Heat vs Control)",
       title = "Differential TF regulon activity in astrocytes",
       color = "Module")

ggsave(file.path(global_dir, "DiffRegulonsQuadrantGlobal.pdf"),
       p_dregs, width = 8, height = 7)

# Labeled version: top 25 TFs by adjusted p-value
top_idx <- order(dregs_plot$pvaladjpositive,
                 na.last = TRUE)[seq_len(min(25, nrow(dregs_plot)))]
top_sig <- dregs_plot[top_idx, ]

p_dregs_labeled <- p_dregs +
  geom_text_repel(data = top_sig, aes(label = tf),
                  size = 3, max.overlaps = Inf)
ggsave(file.path(global_dir, "DiffRegulonsQuadrantGlobalLabeled.pdf"),
       p_dregs_labeled, width = 8, height = 7)

# Global UMAPs
if ("wnn_astro" %in% Reductions(astroobj)) {
  p_sc <- DimPlot(astroobj, reduction = "wnn_astro",
                  group.by = "astrosubclusterwnn", label = TRUE,
                  pt.size  = 0.5) +
    ggtitle("Astro subclusters")
  p_grp <- DimPlot(astroobj, reduction = "wnn_astro",
                   group.by = "TFgroup", pt.size = 0.5) +
    ggtitle("TFgroup: Heat vs Control")
  ggsave(file.path(global_dir, "UMAPSubclusters.pdf"), p_sc,  width = 8, height = 7)
  ggsave(file.path(global_dir, "UMAPTFgroup.pdf"),     p_grp, width = 8, height = 7)
}
cat("Global plots saved.\n")


# ── 9. MODULE-SPECIFIC ANALYSES ───────────────────────────────────────────────
cat("\n=== STEP 8: MODULE-SPECIFIC ANALYSES ===\n")

for (current_module in MODULES_OF_INTEREST) {

  cat(strrep("=", 70), "\n")
  cat("ANALYZING MODULE:", current_module, "\n")
  cat(strrep("=", 70), "\n")

  moddir <- file.path(OUT_DIR, paste0(current_module, "Analysis"))
  for (d in c("1TrueTFsFiltered", "2AllGenesUnfiltered",
              "3RegulonBarPlots",  "4TFNetworkPlots",
              "5UMAPVisualizations", "6GORegulonListsLoose")) {
    dir.create(file.path(moddir, d), recursive = TRUE, showWarnings = FALSE)
  }

  # ── 9a. Module-specific TF filtering ───────────────────────────────────────
  dregs_with_module <- merge(dregs_manual, modules_simple,
                             by.x = "tf", by.y = "gene_name", all.x = TRUE)
  dregs_with_module <- as.data.frame(dregs_with_module)

  mod_tfs_full <- dregs_with_module[
    !is.na(dregs_with_module$module) &
    dregs_with_module$module == current_module, ]
  mod_tfs_full <- mod_tfs_full[order(mod_tfs_full$pvaladjpositive,
                                     -abs(mod_tfs_full$avglog2FCpositive)), ]
  cat(" -", current_module, "TFs with DREG evidence:", nrow(mod_tfs_full), "\n")

  # Heat-like: positive regulon up AND negative regulon down
  mod_heat_like <- mod_tfs_full[
    mod_tfs_full$avglog2FCpositive > 0 &
    !is.na(mod_tfs_full$avglog2FCnegative) &
    mod_tfs_full$avglog2FCnegative < 0, ]
  mod_heat_like <- mod_heat_like[order(mod_heat_like$pvaladjpositive,
                                       -abs(mod_heat_like$avglog2FCpositive)), ]
  cat(" -", current_module, "Heat-like TFs (pos↑ neg↓):", nrow(mod_heat_like), "\n")

  # TRUE TF filter (JASPAR-validated)
  if (!is.null(true_tfs)) {
    mod_true_tfs <- mod_heat_like[mod_heat_like$tf %in% true_tfs, ]
    mod_true_tfs <- mod_true_tfs[order(mod_true_tfs$pvaladjpositive,
                                       -abs(mod_true_tfs$avglog2FCpositive)), ]
    cat(" - Filtered TRUE TFs:", nrow(mod_true_tfs), "\n")
  } else {
    mod_true_tfs <- data.frame()
    cat(" - Skipping TRUE TF filter (not available)\n")
  }

  # Save CSVs
  write.csv(mod_tfs_full,
            file.path(moddir, "2AllGenesUnfiltered",
                      paste0(current_module, "_DREGTFsAll.csv")),
            row.names = FALSE)
  write.csv(mod_heat_like,
            file.path(moddir, "2AllGenesUnfiltered",
                      paste0(current_module, "_HeatLikeUnfiltered.csv")),
            row.names = FALSE)
  if (nrow(mod_true_tfs) > 0)
    write.csv(mod_true_tfs,
              file.path(moddir, "1TrueTFsFiltered",
                        paste0(current_module, "_HeatLikeTrueTFs.csv")),
              row.names = FALSE)

  # Select key TFs for visualization (TRUE > unfiltered fallback)
  top_true_tfs   <- if (nrow(mod_true_tfs)  > 0)
                      head(mod_true_tfs$tf,  5) else character(0)
  top_unfiltered <- if (nrow(mod_heat_like) > 0)
                      head(mod_heat_like$tf, 5) else character(0)
  key_tfs        <- if (length(top_true_tfs) > 0) top_true_tfs else top_unfiltered

  if (length(key_tfs) == 0) {
    cat("No TFs to plot for", current_module, "— skipping visualization.\n")
    next
  }
  cat("Key TFs for", current_module, ":", paste(key_tfs, collapse = ", "), "\n")

  # ── 9b. Regulon Bar Plots ─────────────────────────────────────────────────
  cat("STEP: Regulon bar plots...\n")
  for (tf in key_tfs) {
    tryCatch({
      p <- RegulonBarPlot(astroobj, selectedtf = tf,
                          cutoff    = 0,
                          wgcnaname = wgcna_name) +
        ggtitle(paste(current_module, "TF regulon:", tf))
      ggsave(file.path(moddir, "3RegulonBarPlots",
                       paste0("RegulonBar_", tf, ".pdf")),
             p, width = 6, height = 5)
    }, error = function(e) cat("  - Failed", tf, ":", e$message, "\n"))
  }

  # ── 9c. TF Network Plots (Gain-weighted) ──────────────────────────────────
  cat("STEP: TF network plots...\n")
  for (tf in key_tfs) {
    tryCatch({
      p <- TFNetworkPlot(astroobj,
                         selectedtfs = tf,
                         depth       = 2,
                         targettype  = "positive",
                         edgeweight  = "Gain",
                         cutoff      = 0.05,
                         wgcnaname   = wgcna_name) +
        ggtitle(paste(current_module, "TF network:", tf))
      ggsave(file.path(moddir, "4TFNetworkPlots",
                       paste0("TFNetwork_", tf, "_Depth2.pdf")),
             p, width = 9, height = 7)
    }, error = function(e) {
      cat("  - Depth 2 failed for", tf, "— trying depth 1...\n")
      tryCatch({
        p <- TFNetworkPlot(astroobj,
                           selectedtfs = tf,
                           depth       = 1,
                           targettype  = "positive",
                           edgeweight  = "Gain",
                           cutoff      = 0.05,
                           wgcnaname   = wgcna_name) +
          ggtitle(paste(current_module, "TF network:", tf, "(depth1)"))
        ggsave(file.path(moddir, "4TFNetworkPlots",
                         paste0("TFNetwork_", tf, "_Depth1.pdf")),
               p, width = 9, height = 7)
      }, error = function(e2)
        cat("  - Also failed depth 1:", e2$message, "\n"))
    })
  }

  # ── 9d. UMAP Visualizations ───────────────────────────────────────────────
  cat("STEP: UMAP visualizations...\n")
  if ("wnn_astro" %in% Reductions(astroobj)) {
    for (tf in key_tfs) {
      tryCatch({
        astroobj$pos_regulon_score <- pos_scores[, tf]

        p_expr <- FeaturePlot(astroobj, features  = tf,
                              reduction = "wnn_astro",
                              cols      = c("lightgrey", "red"),
                              pt.size   = 0.5) +
          ggtitle(paste(tf, "expression"))
        p_pos  <- FeaturePlot(astroobj, features  = "pos_regulon_score",
                              reduction = "wnn_astro",
                              cols      = c("lightgrey", "red"),
                              pt.size   = 0.5) +
          ggtitle(paste(tf, "positive regulon score"))

        if (tf %in% colnames(neg_scores)) {
          astroobj$neg_regulon_score <- neg_scores[, tf]
          p_neg <- FeaturePlot(astroobj, features  = "neg_regulon_score",
                               reduction = "wnn_astro",
                               cols      = c("lightgrey", "seagreen"),
                               pt.size   = 0.5) +
            ggtitle(paste(tf, "negative regulon score"))
          p_combined <- p_expr | p_pos | p_neg
          ggsave(file.path(moddir, "5UMAPVisualizations",
                           paste0("UMAP_", tf, "_ExprPosNegRegulons.pdf")),
                 p_combined, width = 15, height = 5)
        } else {
          p_combined <- p_expr | p_pos
          ggsave(file.path(moddir, "5UMAPVisualizations",
                           paste0("UMAP_", tf, "_ExprPosRegulon.pdf")),
                 p_combined, width = 10, height = 5)
        }
      }, error = function(e) cat("  - Failed", tf, ":", e$message, "\n"))
    }
  }

  # ── 9e. Loose Regulons for GO Enrichment (strategy C, regthresh 0.005) ───
  cat("STEP: Loose regulons for GO enrichment...\n")
  astroobj_loose    <- AssignTFRegulons(astroobj,
                                        strategy  = "C",
                                        regthresh = 0.005,
                                        wgcnaname = wgcna_name)
  tf_regulons_loose <- GetTFRegulons(astroobj_loose, wgcnaname = wgcna_name)

  for (tf in key_tfs) {
    cur <- tf_regulons_loose[tf_regulons_loose$tf == tf, ]
    if (nrow(cur) == 0) { cat("  - No loose regulon for", tf, "\n"); next }
    pos_targets <- unique(cur$gene[cur$Cor > 0])
    neg_targets <- unique(cur$gene[cur$Cor < 0])
    cat(sprintf("  - %s: %d pos, %d neg targets\n",
                tf, length(pos_targets), length(neg_targets)))
    if (length(pos_targets) >= GO_MIN_TARGETS)
      write.table(data.frame(gene = pos_targets),
                  file.path(moddir, "6GORegulonListsLoose",
                            paste0(tf, "_positive_targets.txt")),
                  quote = FALSE, row.names = FALSE,
                  col.names = TRUE, sep = "\t")
    if (length(neg_targets) >= GO_MIN_TARGETS)
      write.table(data.frame(gene = neg_targets),
                  file.path(moddir, "6GORegulonListsLoose",
                            paste0(tf, "_negative_targets.txt")),
                  quote = FALSE, row.names = FALSE,
                  col.names = TRUE, sep = "\t")
  }

  # ── 9f. Analysis summary text file ───────────────────────────────────────
  summary_file <- file.path(moddir, paste0(current_module, "_AnalysisSummary.txt"))
  sink(summary_file)
  cat("ASTROCYTE TF ANALYSIS SUMMARY —", current_module, "\n\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("STATISTICS\n")
  cat(" - Total module TFs with DREG evidence:", nrow(mod_tfs_full), "\n")
  cat(" - Heat-like TFs (unfiltered):          ", nrow(mod_heat_like), "\n")
  cat(" - TRUE TFs (filtered):                 ", nrow(mod_true_tfs), "\n")
  cat("KEY TFs FOR PLOTTING\n")
  for (i in seq_along(key_tfs)) cat(i, ".", key_tfs[i], "\n")
  cat("\nOUTPUT FILES\n")
  cat(" - CSVs in 1TrueTFsFiltered and 2AllGenesUnfiltered\n")
  cat(" - Regulon bar plots in 3RegulonBarPlots\n")
  cat(" - TF network plots in 4TFNetworkPlots\n")
  cat(" - UMAP visualizations in 5UMAPVisualizations\n")
  cat(" - GO target lists in 6GORegulonListsLoose\n")
  sink()
  cat("Summary saved:", summary_file, "\n")
}


# ── 10. ASCL1 FOCAL NETWORK — FIGURE 3D ──────────────────────────────────────
cat("\n=== STEP 9: ASCL1 FOCAL NETWORK (Figure 3D) ===\n")

# Use strategy C loose regulons for ASCL1
# tf_regulons_loose from the last module loop (Astro-M7); reassign if needed
astrotmp <- astroobj
if (is.null(astrotmp@misc[[wgcna_name]])) {
  astrotmp@misc[[wgcna_name]] <- astroobj@misc[[wgcna_name]]
}
astrotmp@misc[[wgcna_name]]$TFregulons <- tf_regulons_loose

cur_tf      <- "ASCL1"
ascl1_loose <- tf_regulons_loose[tf_regulons_loose$tf == cur_tf, ]
ascl1_top   <- ascl1_loose[order(-abs(ascl1_loose$Cor)), ]

# Label top 5 targets by |Cor| + always include NOTCH1
top_targets_to_label <- head(ascl1_top$gene, 5)
label_genes          <- unique(c("NOTCH1", top_targets_to_label))

p_ascl1 <- tryCatch(
  TFNetworkPlot(
    astrotmp,
    selectedtfs  = cur_tf,
    depth        = 2,
    keep         = "NOTCH1",       # preserve NOTCH1 at depth 2
    targettype   = "both",
    edgeweight   = "Cor",
    cutoff       = 0.05,
    colorcutoff  = 0.75,
    highcolor    = "firebrick3",
    midcolor     = "grey95",
    lowcolor     = "steelblue3",
    wgcnaname    = wgcna_name,
    labelTFs     = 0,
    labelgenes   = NULL
  ) + ggtitle("ASCL1 TF network (strategy C, depth 2, Cor edges)"),
  error = function(e) {
    cat("ASCL1 TFNetworkPlot failed:", e$message, "\n")
    NULL
  }
)

if (!is.null(p_ascl1)) {
  node_data   <- p_ascl1$data
  label_nodes <- node_data %>% dplyr::filter(name %in% label_genes)
  p_ascl1_final <- p_ascl1 +
    ggrepel::geom_label_repel(
      data         = label_nodes,
      aes(x = x, y = y, label = name),
      size         = 3,
      max.overlaps = Inf,
      fill         = "white",
      label.size   = 0.1
    )
  ggsave(file.path(OUT_DIR, "ASCL1_TFNetwork_strategyCdepth2Cor_labeled.pdf"),
         p_ascl1_final, width = 8, height = 5, units = "in", dpi = 300)
  cat("Saved: ASCL1_TFNetwork_strategyCdepth2Cor_labeled.pdf\n")
}


# ── 11. TOM-SCAFFOLD GRN (ggraph) — FIGURE 3E ────────────────────────────────
cat("\n=== STEP 10: TOM-SCAFFOLD GRN PLOT (Figure 3E) ===\n")

# Build TF-target regulon edges and TOM scaffold edges for ggraph
for (tf in c("ASCL1")) {

  tf_edges <- tf_regulons_loose %>%
    filter(.data$tf == tf) %>%
    transmute(from = tf, to = gene, Cor = Cor, Gain = Gain,
              edgeclass = "TFreg")

  # TOM scaffold: gene–gene edges from WGCNA TOM for the same gene set
  tom_edges <- tryCatch({
    all_nodes <- unique(c(tf, tf_edges$to))
    tom_mat   <- GetTOM(astroobj, wgcnaname = wgcna_name)
    if (!is.null(tom_mat)) {
      shared <- intersect(all_nodes, rownames(tom_mat))
      tm     <- tom_mat[shared, shared]
      tom_long <- as.data.frame(as.table(as.matrix(tm))) %>%
        filter(as.character(Var1) < as.character(Var2), Freq > 0.05) %>%
        transmute(from = as.character(Var1),
                  to   = as.character(Var2),
                  Cor  = NA_real_, Gain = NA_real_,
                  edgeclass = "TOMscaffold",
                  tom = Freq)
      tom_long
    } else NULL
  }, error = function(e) NULL)

  edges_all <- if (!is.null(tom_edges)) {
    bind_rows(
      tf_edges   %>% mutate(tom = NA_real_),
      tom_edges
    )
  } else {
    tf_edges %>% mutate(tom = NA_real_)
  }

  g <- igraph::graph_from_data_frame(edges_all, directed = TRUE)

  set.seed(1)
  p_tom <- ggraph(g, layout = "stress") +
    geom_edge_link(
      data      = function(x) subset(x, edgeclass == "TOMscaffold"),
      aes(alpha = tom),
      color = "grey80", width = 0.3
    ) +
    geom_edge_link(
      data      = function(x) subset(x, edgeclass == "TFreg"),
      aes(edge_colour = Cor, edge_width = Gain, edge_alpha = abs(Cor)),
      arrow = arrow(length = unit(1, "mm"), type = "closed")
    ) +
    scale_edge_colour_gradient2(low = "steelblue3", mid = "grey95",
                                high = "firebrick3") +
    scale_edge_width_continuous(range = c(0.2, 1.4)) +
    scale_edge_alpha_continuous(range = c(0.1, 0.9), guide = "none") +
    geom_node_point(size = 2.2, color = "black", fill = "white") +
    theme_void() +
    labs(edge_colour = "TF→target Cor", edge_width = "Gain") +
    ggtitle(paste("TF network with TOM scaffold —", tf))

  ggsave(file.path(OUT_DIR,
                   paste0(tf, "_TOMScaffoldNetwork.pdf")),
         p_tom, width = 9, height = 7)
  cat("Saved:", tf, "TOM scaffold network\n")
}


# ── 12. VERIFICATION CHECKLIST ───────────────────────────────────────────────
cat("\n=== FILE VERIFICATION ===\n")

files_expected <- c(
  file.path(OUT_DIR,     "0GlobalAnalysis", "DiffRegulonsQuadrantGlobal.pdf"),
  file.path(OUT_DIR,     "0GlobalAnalysis", "DiffRegulonsQuadrantGlobalLabeled.pdf"),
  file.path(OUT_DIR,     "0GlobalAnalysis", "UMAPSubclusters.pdf"),
  file.path(OUT_DIR,     "0GlobalAnalysis", "UMAPTFgroup.pdf"),
  file.path(OUT_DIR,     "ASCL1_TFNetwork_strategyCdepth2Cor_labeled.pdf"),
  file.path(OUT_DIR,     "ASCL1_TOMScaffoldNetwork.pdf"),
  file.path(OUT_DIR,     "Astro-M7Analysis",
            "1TrueTFsFiltered", "Astro-M7_HeatLikeTrueTFs.csv"),
  file.path(OUT_DIR,     "Astro-M7Analysis",
            "2AllGenesUnfiltered", "Astro-M7_HeatLikeUnfiltered.csv")
)

for (f in files_expected)
  cat(sprintf("  %-65s %s\n", basename(f),
              ifelse(file.exists(f), "EXISTS", "MISSING")))


# ── 13. FINAL SUMMARY ────────────────────────────────────────────────────────
cat("\n========== 11_TF_Regulatory_Network_GRN.R COMPLETE ==========\n")
cat("wgcna_name:                    ", wgcna_name, "\n")
cat("TF regulon pairs:              ", nrow(tf_regulons_full), "\n")
cat("TFs in dregs_manual:           ", nrow(dregs_manual), "\n")
cat("Modules analyzed:              ",
    paste(MODULES_OF_INTEREST, collapse = ", "), "\n")
cat("\nKey outputs:\n")
cat("  - 0GlobalAnalysis/DiffRegulonsQuadrantGlobal[Labeled].pdf\n")
cat("  - 0GlobalAnalysis/UMAP*.pdf\n")
for (m in MODULES_OF_INTEREST)
  cat(sprintf("  - %sAnalysis/ (bar plots, network plots, GO lists)\n", m))
cat("  - ASCL1_TFNetwork_strategyCdepth2Cor_labeled.pdf   (Figure 3D)\n")
cat("  - ASCL1_TOMScaffoldNetwork.pdf                     (Figure 3E)\n")
cat("\nNext step: run 12_Pseudotime_Trajectory_Analysis.R\n")
