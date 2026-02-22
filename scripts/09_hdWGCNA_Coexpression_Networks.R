################################################################################
# Script:  09_hdWGCNA_Coexpression_Networks.R
# Project: Heat_Call_Single_Nucleus_Multiome
# Author:  Prakrit Subba
# Date:    2026
#
# Paper:   Subba et al. 2026 — Multiome Profiling Reveals Astrocyte and
#          Neuroendocrine Targets of Prenatal Acoustic Programming in
#          Zebra Finch Embryos
#
# Description:
#   Weighted gene co-expression network analysis (hdWGCNA) on the astrocyte
#   subset, followed by TF regulon scoring and differential regulon analysis
#   between Heat Call and Control conditions.
#
#   Steps:
#     (1) hdWGCNA network construction on astrocytes
#         SetupForWGCNA → MetacellsByGroups → NormalizeMetacells →
#         SetDatExpr → TestSoftPowers → ConstructNetwork →
#         ModuleEigengenes → ModuleConnectivity
#     (2) Module trait correlation → AstroModuleTraitCorrelationsRNA.csv
#     (3) TF regulon assignment (AssignTFRegulons, strategy = "A")
#         and regulon scoring (positive + negative)
#     (4) Differential regulon analysis: Heat vs Control
#         Method 1 (manual): FindMarkers on regulon score matrix
#         Method 2 (built-in): FindDifferentialRegulons
#     (5) Figure 4 hdWGCNA TF-Module networks:
#         Panel A — ModuleRegulatoryHeatmap (delta)
#         Panel B — ModuleRegulatoryNetworkPlot (delta)
#     (6) Astrocyte subcluster Astro-M7 score violin plot
#
# Inputs:
#   - results/08_ChromVAR/08_astroobj_chromvar.rds
#   - AstrohdWGCNARNAWorkspaceNov2025.RData  (pre-computed hdWGCNA workspace)
#   - AstroTFNetworkCheckpoint.rds           (TF network checkpoint for Fig 4)
#   - subclusters.rds                        (astrocyte sub-clustering object)
#   - data/jaspar_mapping.rds                (optional JASPAR TF-gene mapping)
#
# Outputs:
#   Module trait correlation:
#     AstroModuleTraitCorrelationsRNA.csv
#     AstroModuleTraitHeatmapRNA.pdf
#     AstroModuleTraitHeatmapRNA.svg
#
#   Regulon scores:
#     AstroPositiveRegulonScoresAstroNovTruth.rds
#     AstroNegativeRegulonScoresAstroNovTruth.rds
#
#   Differential regulons:
#     AstroDiffRegulonsManualHeatVsControlAllTFsNovTruth.csv
#     AstroDiffRegulonshdWGCNAHeatVsControlNovTruth.csv
#
#   Figure 4:
#     FIGURES/Fig4ATFModuleDeltaHeatmap.pdf
#     FIGURES/Fig4BTFModuleNetwork.pdf
#
#   Astro-M7 module score:
#     AstroM7ScoreViolinBySubcluster.pdf
#     AstroM7ScoreViolinBySubcluster.png
#
# NOTE on naming:
#   wgcna_name              = "AstroAnnotation"       (no spaces)
#   Module name             = "Astro-M7"              (hyphen, not underscore)
#   Column TFgroup values   = "Heat" / "Control"      (NOT "Heat Call")
#   Conditionnum encoding   : 0 = Heat Call, 1 = Control
#   Diff regulon columns    : avglog2FCpositive / avglog2FCnegative
#                             pvaladjpositive  / pvaladjnegative
#                             (NO underscores between term and direction)
#   Violin PDF device       : cairo_pdf                (exact)
#
# Previous Script: 08_ChromVAR_Motif_Analysis.R
# Next Script:     10_Cicero_Chromatin_Rewiring.R
################################################################################


################################################################################
# 0. CONFIGURATION
################################################################################

ASTRO_RDS      <- "results/08_ChromVAR/08_astroobj_chromvar.rds"
WNN_RDS        <- "results/06_DA/06_combined_wnn_DA.rds"

# Pre-computed hdWGCNA workspace (constructed Nov 2025)
WGCNA_WORKSPACE   <- "AstrohdWGCNARNAWorkspaceNov2025.RData"
TF_CHECKPOINT_RDS <- "AstroTFNetworkCheckpoint.rds"
SUBCLUSTERS_RDS   <- "subclusters.rds"

# Optional JASPAR TF-gene mapping
JASPAR_MAPPING_RDS <- "data/jaspar_mapping.rds"

OUT_DIR  <- "results/09_hdWGCNA"
DIR_FIG  <- file.path(OUT_DIR, "FIGURES")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG, showWarnings = FALSE)

# hdWGCNA parameters
wgcna_name   <- "AstroAnnotation"    # Exact WGCNA slot name throughout
SOFT_POWER   <- 12                   # Soft-thresholding power from TestSoftPowers
FDR_THRESHOLD <- 0.05

set.seed(12345)


################################################################################
# 1. LOAD LIBRARIES
################################################################################

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
})

theme_set(theme_cowplot())


################################################################################
# 2. LOAD OBJECTS
################################################################################

cat("=== Loading objects ===\n")

combined_wnn <- readRDS(WNN_RDS)
astroobj     <- readRDS(ASTRO_RDS)

cat("combined_wnn — cells:", ncol(combined_wnn), "\n")
cat("astroobj — cells    :", ncol(astroobj),     "\n\n")

stopifnot(
  "cluster_annotation" %in% colnames(combined_wnn@meta.data),
  "Condition"          %in% colnames(combined_wnn@meta.data)
)

# Set active WGCNA slot name on both objects
combined_wnn@misc$active_wgcna  <- wgcna_name
DefaultAssay(combined_wnn)      <- "RNA"
DefaultAssay(astroobj)          <- "RNA"


################################################################################
# 3. LOAD PRE-COMPUTED hdWGCNA WORKSPACE
#
# The hdWGCNA network (metacells, TOM, modules) was constructed in a prior
# session and saved as AstrohdWGCNARNAWorkspaceNov2025.RData.
# Loading this workspace defines `combined_wnn` with the hdWGCNA slot already
# populated (modules, eigengenes, connectivity).
# If the workspace does not exist, run Section 4 (network construction) first.
################################################################################

if (file.exists(WGCNA_WORKSPACE)) {

  cat("=== Loading pre-computed hdWGCNA workspace ===\n")
  cat("  File:", WGCNA_WORKSPACE, "\n")

  load(WGCNA_WORKSPACE)    # defines combined_wnn with hdWGCNA slot populated

  stopifnot(exists("combined_wnn"))
  combined_wnn@misc$active_wgcna <- wgcna_name
  DefaultAssay(combined_wnn)     <- "RNA"

  cat("  Workspace loaded — cells:", ncol(combined_wnn), "\n\n")

} else {
  cat("NOTE: Workspace not found. Running network construction (Section 4).\n\n")
}


################################################################################
# 4. HDWGCNA NETWORK CONSTRUCTION (run if workspace absent)
#
# Standard hdWGCNA pipeline on astrocyte subset only.
# Metacells are constructed grouped by cluster_annotation × Condition to
# preserve condition-level variation. Network is built from SCT-normalized
# RNA assay. SOFT_POWER selected from TestSoftPowers elbow plot.
#
# Skip this section if AstrohdWGCNARNAWorkspaceNov2025.RData already exists.
################################################################################

if (!file.exists(WGCNA_WORKSPACE)) {

  cat("=== Running hdWGCNA network construction ===\n")
  cat("  (Computationally intensive — run on HPC)\n\n")

  # Subset to astrocytes for construction
  astro_seurat <- subset(combined_wnn, subset = cluster_annotation == "Astrocyte")
  DefaultAssay(astro_seurat) <- "RNA"

  # 4a. Setup
  astro_seurat <- SetupForWGCNA(
    astro_seurat,
    gene_select = "fraction",
    fraction    = 0.05,
    wgcna_name  = wgcna_name
  )

  # 4b. Metacells
  astro_seurat <- MetacellsByGroups(
    seurat_obj      = astro_seurat,
    group.by        = c("cluster_annotation", "Condition"),
    reduction       = "pca",
    k               = 25,
    target_metacells = 250,
    assay           = "RNA",
    wgcna_name      = wgcna_name
  )

  astro_seurat <- NormalizeMetacells(astro_seurat)

  # 4c. Expression matrix
  astro_seurat <- SetDatExpr(
    astro_seurat,
    group_name  = "Astrocyte",
    group.by    = "cluster_annotation",
    assay       = "RNA",
    slot        = "data",
    wgcna_name  = wgcna_name
  )

  # 4d. Soft power selection
  astro_seurat <- TestSoftPowers(
    astro_seurat,
    networkType = "signed",
    wgcna_name  = wgcna_name
  )

  sp_table   <- GetPowerTable(astro_seurat)
  sp_plot    <- PlotSoftPowers(astro_seurat)

  pdf(file.path(OUT_DIR, "AstroSoftPowerPlot.pdf"), width = 10, height = 6)
  print(sp_plot)
  dev.off()

  cat("  Soft power selected:", SOFT_POWER,
      "(verify from AstroSoftPowerPlot.pdf)\n\n")

  # 4e. Construct network
  astro_seurat <- ConstructNetwork(
    astro_seurat,
    soft_power   = SOFT_POWER,
    networkType  = "signed",
    tom_name     = "Astro",
    wgcna_name   = wgcna_name
  )

  # 4f. Module eigengenes
  astro_seurat <- ModuleEigengenes(
    astro_seurat,
    group.by.vars = "Condition",
    wgcna_name    = wgcna_name
  )

  # 4g. Module connectivity
  astro_seurat <- ModuleConnectivity(
    astro_seurat,
    group.by    = "cluster_annotation",
    group_name  = "Astrocyte",
    wgcna_name  = wgcna_name
  )

  # Copy hdWGCNA slot back to combined_wnn
  combined_wnn@misc[[wgcna_name]] <- astro_seurat@misc[[wgcna_name]]
  combined_wnn@misc$active_wgcna  <- wgcna_name

  # Save workspace
  save(combined_wnn,
       file = WGCNA_WORKSPACE)

  cat("Workspace saved:", WGCNA_WORKSPACE, "\n\n")
}


################################################################################
# 5. MODULE EIGENGENES AND TRAIT CORRELATION
#
# ModuleTraitCorrelation() correlates harmonized module eigengenes (hMEs)
# with sample-level traits (Condition, Sex).
# PlotModuleTraitCorrelation() generates the heatmap figure.
#
# Outputs:
#   AstroModuleTraitCorrelationsRNA.csv   (correlation matrix + p-values)
#   AstroModuleTraitHeatmapRNA.pdf
#   AstroModuleTraitHeatmapRNA.svg
################################################################################

cat("=== Module eigengenes and trait correlation ===\n")

hMEs <- GetMEs(combined_wnn)

cat("hMEs dimensions:", paste(dim(hMEs), collapse = " x "), "\n\n")

# Traits: Condition and Sex (numeric encoding required by ModuleTraitCorrelation)
combined_wnn$ConditionNum <- ifelse(combined_wnn$Condition == "Heat Call", 0, 1)

combined_wnn <- ModuleTraitCorrelation(
  combined_wnn,
  traits     = c("ConditionNum"),
  group.by   = "cluster_annotation",
  wgcna_name = wgcna_name
)

# Retrieve and save correlation table
mt_results <- GetModuleTraitCorrelation(combined_wnn)

write.csv(
  mt_results$cor,
  file.path(OUT_DIR, "AstroModuleTraitCorrelationsRNA.csv")
)

cat("Saved: AstroModuleTraitCorrelationsRNA.csv\n\n")

# Heatmap figure
p_mt <- PlotModuleTraitCorrelation(
  combined_wnn,
  label     = "fdr",
  label_symbol = "stars",
  trait_name = "ConditionNum",
  combine    = TRUE,
  wgcna_name = wgcna_name
)

pdf(file.path(OUT_DIR, "AstroModuleTraitHeatmapRNA.pdf"), width = 10, height = 8)
print(p_mt)
dev.off()

svg(file.path(OUT_DIR, "AstroModuleTraitHeatmapRNA.svg"), width = 10, height = 8)
print(p_mt)
dev.off()

cat("Saved: AstroModuleTraitHeatmapRNA.pdf / .svg\n\n")


################################################################################
# 6. SET UP TFGROUP METADATA COLUMN
#
# TFgroup encodes the condition for use in FindDifferentialRegulons and
# regulon score comparisons. Values are "Heat" and "Control" (NOT "Heat Call").
#
# Encoding (exact from original code):
#   Conditionnum = 0 → Heat Call → TFgroup = "Heat"
#   Conditionnum = 1 → Control   → TFgroup = "Control"
#
# NOTE: ifelse(Conditionnum > 0, "Control", "Heat") — Control when > 0.
################################################################################

cat("=== Setting up TFgroup metadata column ===\n")

astrocells <- colnames(combined_wnn)[combined_wnn$cluster_annotation == "Astrocyte"]
stopifnot(length(astrocells) > 0)

combined_wnn$TFgroup <- NA_character_

if ("Conditionnum" %in% colnames(combined_wnn@meta.data)) {
  combined_wnn$TFgroup[astrocells] <- ifelse(
    combined_wnn$Conditionnum[astrocells] > 0,
    "Control",
    "Heat"
  )
} else if ("Condition" %in% colnames(combined_wnn@meta.data)) {
  combined_wnn$TFgroup[astrocells] <- as.character(
    combined_wnn$Condition[astrocells]
  )
} else {
  stop("No Conditionnum or Condition column found in metadata.")
}

cat("TFgroup counts (full object):\n")
print(table(combined_wnn$TFgroup, useNA = "ifany"))
cat("\n")

# Propagate TFgroup to astroobj subset
astroobj <- subset(combined_wnn, cells = astrocells)
DefaultAssay(astroobj)          <- "RNA"
astroobj@misc$WGCNA$name        <- wgcna_name
astroobj@misc$active_wgcna      <- wgcna_name

cat("Astrocyte subset — cells:", ncol(astroobj), "\n")
cat("TFgroup in astroobj:\n")
print(table(astroobj$TFgroup, useNA = "ifany"))
cat("\n")


################################################################################
# 7. LOAD ASTROCYTE SUBCLUSTER OBJECT
#
# astrosub contains astrocyte sub-clustering results (wsnn_res.0.06 resolution)
# and a wnn_astro reduction for UMAP visualization. Subcluster labels are
# copied to combined_wnn as the $astrosubcluster metadata column.
################################################################################

cat("=== Loading astrocyte subcluster object ===\n")

if (file.exists(SUBCLUSTERS_RDS)) {

  astrosub <- readRDS(SUBCLUSTERS_RDS)

  stopifnot(
    "wsnn_res.0.06" %in% colnames(astrosub@meta.data),
    all(colnames(astrosub) %in% colnames(combined_wnn))
  )

  # Copy subcluster labels to combined_wnn
  combined_wnn$astrosubcluster <- NA_character_
  combined_wnn$astrosubcluster[colnames(astrosub)] <-
    as.character(astrosub$wsnn_res.0.06)

  cat("Subclusters:", length(unique(astrosub$wsnn_res.0.06)), "\n")
  print(table(combined_wnn$astrosubcluster, useNA = "ifany"))
  cat("\n")

  # Copy wnn_astro UMAP embedding to astrosub if present
  if ("wnn_astro" %in% Reductions(astrosub)) {
    emb <- Embeddings(astrosub, "wnn_astro")
    cat("wnn_astro embedding found:", nrow(emb), "cells\n\n")
  }

} else {
  cat("NOTE: subclusters.rds not found — skipping subcluster steps.\n")
  astrosub <- NULL
}


################################################################################
# 8. TF REGULON ASSIGNMENT AND SCORING
#
# hdWGCNA TF regulon analysis integrates WGCNA co-expression modules with
# TF-target gene relationships from motif data.
#
# AssignTFRegulons parameters (exact from original code):
#   strategy  = "A"      (TF-module assignment by correlation)
#   regthresh = 0.01     (minimum regulon correlation threshold)
#   ntfs      = 10       (maximum TFs per module)
#
# RegulonScores computed for both positive (activation) and negative
# (repression) targets. ncores = 1 for reproducibility.
#
# Outputs:
#   AstroPositiveRegulonScoresAstroNovTruth.rds
#   AstroNegativeRegulonScoresAstroNovTruth.rds
################################################################################

cat("=== TF regulon assignment and scoring ===\n")

# Check for existing TF regulons
tf_regulons_full <- tryCatch(
  GetTFRegulons(combined_wnn),
  error = function(e) NULL
)

if (is.null(tf_regulons_full) || nrow(tf_regulons_full) == 0) {
  cat("  TF regulons absent — running AssignTFRegulons...\n")

  combined_wnn <- AssignTFRegulons(
    combined_wnn,
    strategy  = "A",
    regthresh = 0.01,
    ntfs      = 10
  )

  tf_regulons_full <- GetTFRegulons(combined_wnn)

  if (is.null(tf_regulons_full) || nrow(tf_regulons_full) == 0) {
    stop("AssignTFRegulons failed — cannot retrieve TF regulons.")
  }
}

cat("TF regulon pairs:", nrow(tf_regulons_full), "\n\n")

# Propagate WGCNA slot and regulons to astroobj
if (is.null(astroobj@misc$WGCNA)) {
  astroobj@misc$WGCNA <- list()
}

astroobj@misc$WGCNA[[wgcna_name]]               <- combined_wnn@misc[[wgcna_name]]
astroobj@misc$WGCNA[[wgcna_name]]$TFregulons    <- tf_regulons_full
astroobj@misc$active_wgcna                       <- wgcna_name

# Compute regulon scores: positive and negative targets
cat("  Computing positive regulon scores (ncores = 1)...\n")
astroobj <- RegulonScores(astroobj, targettype = "positive", ncores = 1)

cat("  Computing negative regulon scores (ncores = 1, corthresh = -0.05)...\n")
astroobj <- RegulonScores(astroobj, targettype = "negative",
                          corthresh = -0.05, ncores = 1)

pos_scores <- GetRegulonScores(astroobj, targettype = "positive")
neg_scores <- GetRegulonScores(astroobj, targettype = "negative")

cat("Positive regulon scores dimensions:",
    paste(dim(pos_scores), collapse = " x "), "\n")
cat("Negative regulon scores dimensions:",
    paste(dim(neg_scores), collapse = " x "), "\n\n")

saveRDS(pos_scores,
        file.path(OUT_DIR, "AstroPositiveRegulonScoresAstroNovTruth.rds"))
saveRDS(neg_scores,
        file.path(OUT_DIR, "AstroNegativeRegulonScoresAstroNovTruth.rds"))

cat("Saved: AstroPositiveRegulonScoresAstroNovTruth.rds\n")
cat("Saved: AstroNegativeRegulonScoresAstroNovTruth.rds\n\n")


################################################################################
# 9. DIFFERENTIAL REGULON ANALYSIS — HEAT VS CONTROL
#
# Two complementary methods:
#
# METHOD 1 (manual — fast):
#   Creates temporary Seurat assays from regulon score matrices (t(pos_scores),
#   t(neg_scores)), then runs FindMarkers (Wilcoxon) comparing Heat vs Control
#   astrocytes. cells.1 = heat_astro, cells.2 = control_astro.
#
#   Parameters: slot = "counts", test.use = "wilcox",
#               logfc.threshold = 0, min.pct = 0
#
#   Output column names (NO underscores between term and direction):
#     dregs_posmanual : tf, avglog2FCpositive, pvalpositive, pvaladjpositive
#     dregs_negmanual : tf, avglog2FCnegative, pvalnegative, pvaladjnegative
#     dregs_manual    : merged; includes ntargets from regulon summary
#
# METHOD 2 (built-in hdWGCNA):
#   FindDifferentialRegulons(group.by = "TFgroup", group1 = "Heat",
#                             group2 = "Control", test.use = "wilcox")
#
# Outputs:
#   AstroDiffRegulonsManualHeatVsControlAllTFsNovTruth.csv
#   AstroDiffRegulonshdWGCNAHeatVsControlNovTruth.csv
################################################################################

cat("=== Differential regulon analysis — Heat vs Control ===\n")

astro_cells_obj <- colnames(astroobj)
control_astro   <- astro_cells_obj[astroobj$TFgroup == "Control"]
heat_astro      <- astro_cells_obj[astroobj$TFgroup == "Heat"]

cat("Controls:", length(control_astro), "| Heat:", length(heat_astro), "\n\n")

# --- METHOD 1: Manual FindMarkers on regulon score matrices ---

cat("  Method 1: Manual FindMarkers on regulon score matrix...\n")

# Positive regulon differential activity
pos_assay <- CreateAssayObject(counts = t(pos_scores))

dregs_posmanual <- FindMarkers(
  object          = pos_assay,
  cells.1         = heat_astro,
  cells.2         = control_astro,
  slot            = "counts",
  test.use        = "wilcox",
  logfc.threshold = 0,
  min.pct         = 0
)

# Add tf column and rename (exact column names — NO underscores between term+direction)
dregs_posmanual$tf  <- rownames(dregs_posmanual)
dregs_posmanual     <- dregs_posmanual[, c("tf", "avg_log2FC", "p_val", "p_val_adj")]
colnames(dregs_posmanual) <- c("tf", "avglog2FCpositive", "pvalpositive", "pvaladjpositive")

cat("  Positive regulons done. Running negative regulons...\n")

# Negative regulon differential activity
neg_assay <- CreateAssayObject(counts = t(neg_scores))

dregs_negmanual <- FindMarkers(
  object          = neg_assay,
  cells.1         = heat_astro,
  cells.2         = control_astro,
  slot            = "counts",
  test.use        = "wilcox",
  logfc.threshold = 0,
  min.pct         = 0
)

dregs_negmanual$tf  <- rownames(dregs_negmanual)
dregs_negmanual     <- dregs_negmanual[, c("tf", "avg_log2FC", "p_val", "p_val_adj")]
colnames(dregs_negmanual) <- c("tf", "avglog2FCnegative", "pvalnegative", "pvaladjnegative")

cat("  Negative regulons done. Combining results...\n")

# Merge positive + negative
dregs_manual <- merge(dregs_posmanual, dregs_negmanual, by = "tf", all = TRUE)

# Add regulon summary (ntargets, mean_gain, mean_cor)
reg_summary <- tf_regulons_full %>%
  group_by(tf) %>%
  summarise(
    ntargets  = n(),
    mean_gain = mean(Gain,   na.rm = TRUE),
    mean_cor  = mean(Cor,    na.rm = TRUE),
    .groups   = "drop"
  )

dregs_manual <- as.data.frame(dregs_manual)
reg_summary  <- as.data.frame(reg_summary)
dregs_manual <- merge(
  dregs_manual,
  reg_summary[, c("tf", "ntargets")],
  by     = "tf",
  all.x  = TRUE
)

write.csv(
  dregs_manual,
  file.path(OUT_DIR, "AstroDiffRegulonsManualHeatVsControlAllTFsNovTruth.csv"),
  row.names = FALSE
)

cat("  TFs analyzed (manual method):", nrow(dregs_manual), "\n")
cat("  Saved: AstroDiffRegulonsManualHeatVsControlAllTFsNovTruth.csv\n\n")

# Preview top 20 hits
dregs_manual_sorted <- dregs_manual[order(dregs_manual$pvaladjpositive,
                                          na.last = TRUE), ]
top_hits <- head(dregs_manual_sorted, 20)
cat("Top 20 differential TFs (manual method):\n")
print(top_hits[, c("tf", "avglog2FCpositive", "pvaladjpositive", "ntargets")])
cat("\n")


# --- METHOD 2: Built-in hdWGCNA FindDifferentialRegulons ---

cat("  Method 2: Built-in hdWGCNA FindDifferentialRegulons...\n")

dregs_hdwgcna <- tryCatch(
  FindDifferentialRegulons(
    astroobj,
    group.by   = "TFgroup",
    group1     = "Heat",
    group2     = "Control",
    test.use   = "wilcox",
    wgcna_name = wgcna_name
  ),
  error = function(e) {
    cat("  Built-in method failed:", conditionMessage(e), "\n")
    cat("  Using manual method results only.\n\n")
    NULL
  }
)

if (!is.null(dregs_hdwgcna)) {
  cat("  Built-in method TFs analyzed:", nrow(dregs_hdwgcna), "\n")

  write.csv(
    dregs_hdwgcna,
    file.path(OUT_DIR, "AstroDiffRegulonshdWGCNAHeatVsControlNovTruth.csv"),
    row.names = FALSE
  )

  cat("  Saved: AstroDiffRegulonshdWGCNAHeatVsControlNovTruth.csv\n\n")
}


################################################################################
# 10. OPTIONAL: JASPAR TRUE-TF FILTER
#
# If jaspar_mapping.rds exists, loads a TF-gene-motif mapping table to
# identify which genes in the dataset are bona fide TFs with JASPAR motifs.
# Used to filter regulon results to known TFs only.
################################################################################

if (file.exists(JASPAR_MAPPING_RDS)) {

  cat("=== Loading JASPAR TF filter ===\n")

  tf_mapping       <- readRDS(JASPAR_MAPPING_RDS)
  true_tfs         <- unique(tf_mapping$tf_gene)

  cat("TRUE TFs from JASPAR:", length(true_tfs), "\n")

} else {
  cat("NOTE: jaspar_mapping.rds not found — proceeding without TF filter.\n")
  true_tfs <- NULL
}

cat("\n")


################################################################################
# 11. FIGURE 4 — HDWGCNA TF-MODULE NETWORKS
#
# Loads AstroTFNetworkCheckpoint.rds (pre-computed with RegulonScores for
# all cells) and generates two panels:
#
# Panel A — ModuleRegulatoryHeatmap (delta feature):
#   Heatmap of TF-module correlation differences (HC − Control).
#   scale_fill_gradient2: blue (negative) → white → red (positive)
#   midpoint = 0, limits = c(-0.5, 0.5)
#   Output: FIGURES/Fig4ATFModuleDeltaHeatmap.pdf (10 × 8 inches)
#
# Panel B — ModuleRegulatoryNetworkPlot (delta feature):
#   Network graph of TF-module connections, cutoff = 0.3, maxval = 1.0
#   Output: FIGURES/Fig4BTFModuleNetwork.pdf (10 × 8 inches)
################################################################################

cat("=== Generating Figure 4 — hdWGCNA TF-Module Networks ===\n")

if (file.exists(TF_CHECKPOINT_RDS)) {

  cat("  Loading:", TF_CHECKPOINT_RDS, "\n")

  astroobj_fig4 <- readRDS(TF_CHECKPOINT_RDS)

  wgcna_name                          <- "AstroAnnotation"
  DefaultAssay(astroobj_fig4)         <- "RNA"
  astroobj_fig4@misc$WGCNA$name      <- wgcna_name

  # Subset to Astrocyte cells only for TF network plots
  astrocells_fig4   <- colnames(astroobj_fig4)[
    astroobj_fig4$cluster_annotation == "Astrocyte"
  ]
  astro_only_wgcna  <- subset(astroobj_fig4, cells = astrocells_fig4)
  astro_only_wgcna@misc$WGCNA$name <- wgcna_name

  # TFgroup for astro_only_wgcna
  if (!"TFgroup" %in% colnames(astro_only_wgcna@meta.data)) {
    astro_only_wgcna$TFgroup <- ifelse(
      astro_only_wgcna$Conditionnum > 0,
      "Control",
      "Heat"
    )
  }

  # Ensure TF regulons are present
  tf_regs <- GetTFRegulons(astro_only_wgcna)
  if (is.null(tf_regs)) {
    astro_only_wgcna <- AssignTFRegulons(
      astro_only_wgcna,
      strategy  = "A",
      regthresh = 0.01,
      ntfs      = 10
    )
  }

  # Ensure positive regulon scores are present
  pos_check <- tryCatch(
    GetRegulonScores(astro_only_wgcna, targettype = "positive"),
    error = function(e) NULL
  )
  if (is.null(pos_check)) {
    astro_only_wgcna <- RegulonScores(
      astro_only_wgcna,
      targettype = "positive",
      ncores     = 4
    )
  }

  # Panel A: TF-module regulatory heatmap (delta: HC − Control)
  p4a <- ModuleRegulatoryHeatmap(
    astro_only_wgcna,
    feature       = "delta",
    TFsonly       = TRUE,
    dendrogram    = FALSE,
    wgcna_name    = wgcna_name,
    show_rownames = TRUE,
    show_colnames = TRUE
  ) +
    scale_fill_gradient2(
      low      = "blue",
      mid      = "white",
      high     = "red",
      midpoint = 0,
      name     = "Correlation\nHC-Control",
      limits   = c(-0.5, 0.5)
    ) +
    ggtitle("Panel A: TF-Module Regulation Shifts") +
    theme(
      plot.title  = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1)
    )

  ggsave(
    file.path(DIR_FIG, "Fig4ATFModuleDeltaHeatmap.pdf"),
    p4a, width = 10, height = 8
  )

  cat("  Saved: FIGURES/Fig4ATFModuleDeltaHeatmap.pdf\n")

  # Panel B: TF-module network graph (delta feature, cutoff = 0.3)
  p4b <- tryCatch(
    ModuleRegulatoryNetworkPlot(
      astro_only_wgcna,
      feature    = "delta",
      cutoff     = 0.3,
      maxval     = 1.0,
      TFsonly    = TRUE,
      wgcna_name = wgcna_name
    ) +
      ggtitle("Panel B: TF-Module Network"),
    error = function(e) {
      cat("  Panel B failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(p4b)) {
    ggsave(
      file.path(DIR_FIG, "Fig4BTFModuleNetwork.pdf"),
      p4b, width = 10, height = 8
    )
    cat("  Saved: FIGURES/Fig4BTFModuleNetwork.pdf\n\n")
  }

} else {
  cat("NOTE:", TF_CHECKPOINT_RDS, "not found — skipping Figure 4.\n\n")
}


################################################################################
# 12. ASTRO-M7 MODULE SCORE VIOLIN PLOT
#
# Adds the Astro-M7 module eigengene score to astrosub metadata using the
# hMEs matrix (cells × modules). common_cells ensures only cells present
# in both objects are used.
#
# Module name: "Astro-M7" (hyphen — exact column name in hMEs matrix)
# Metadata column added: $AstroM7score
#
# Violin plot saved as:
#   AstroM7ScoreViolinBySubcluster.pdf  (device = cairo_pdf — exact)
#   AstroM7ScoreViolinBySubcluster.png  (dpi = 300)
################################################################################

if (!is.null(astrosub)) {

  cat("=== Adding Astro-M7 module score to subcluster object ===\n")

  common_cells <- intersect(colnames(astrosub), rownames(hMEs))

  cat("Common cells between astrosub and hMEs:", length(common_cells), "\n")

  # Add Astro-M7 score (module name uses hyphen: "Astro-M7")
  astrosub$AstroM7score                   <- NA_real_
  astrosub$AstroM7score[common_cells]     <- hMEs[common_cells, "Astro-M7"]

  cat("Astro-M7 score range:",
      paste(round(range(astrosub$AstroM7score, na.rm = TRUE), 3),
            collapse = " – "), "\n\n")

  # Subcluster UMAP
  if ("wnn_astro" %in% Reductions(astrosub)) {

    n_clusters    <- length(unique(astrosub$astrosubclusterwnn))
    cluster_colors <- scales::hue_pal()(n_clusters)

    p_subclusters <- DimPlot(
      astrosub,
      reduction  = "wnn_astro",
      group.by   = "astrosubclusterwnn",
      label      = TRUE,
      label.size = 6,
      pt.size    = 0.5,
      repel      = TRUE
    ) +
      ggtitle("Astrocyte Subclusters") +
      theme_cowplot(font_size = 14) +
      theme(
        plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "right",
        axis.line       = element_line(linewidth = 0.5),
        axis.text       = element_text(size = 12),
        axis.title      = element_text(size = 14)
      ) +
      labs(color = "Subcluster")

    ggsave(
      file.path(DIR_FIG, "AstroSubclusterUMAP.pdf"),
      p_subclusters, width = 8, height = 7
    )

    cat("  Saved: FIGURES/AstroSubclusterUMAP.pdf\n")
  }

  # Violin plot: Astro-M7 score by subcluster
  p_violin <- ggplot(
    astrosub@meta.data,
    aes(x = astrosubclusterwnn, y = AstroM7score, fill = astrosubclusterwnn)
  ) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_cowplot(font_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 0, hjust = 0.5),
      plot.title      = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      x     = "Astrocyte Subcluster",
      y     = "Astro-M7 Module Score",
      title = "Astro-M7 Module Score Distribution by Subcluster"
    ) +
    stat_summary(
      fun   = mean,
      geom  = "point",
      shape = 23,
      size  = 3,
      fill  = "black"
    )

  # PDF: device = cairo_pdf (exact from original code)
  ggsave(
    file.path(OUT_DIR, "AstroM7ScoreViolinBySubcluster.pdf"),
    p_violin,
    width  = 10,
    height = 6,
    device = cairo_pdf
  )

  # PNG: dpi = 300
  ggsave(
    file.path(OUT_DIR, "AstroM7ScoreViolinBySubcluster.png"),
    p_violin,
    width  = 10,
    height = 6,
    dpi    = 300
  )

  cat("  Saved: AstroM7ScoreViolinBySubcluster.pdf (cairo_pdf)\n")
  cat("  Saved: AstroM7ScoreViolinBySubcluster.png\n\n")
}


################################################################################
# 13. VERIFICATION CHECKLIST
################################################################################

cat("=== File verification ===\n")

files_expected <- c(
  file.path(OUT_DIR, "AstroModuleTraitCorrelationsRNA.csv"),
  file.path(OUT_DIR, "AstroModuleTraitHeatmapRNA.pdf"),
  file.path(OUT_DIR, "AstroModuleTraitHeatmapRNA.svg"),
  file.path(OUT_DIR, "AstroPositiveRegulonScoresAstroNovTruth.rds"),
  file.path(OUT_DIR, "AstroNegativeRegulonScoresAstroNovTruth.rds"),
  file.path(OUT_DIR, "AstroDiffRegulonsManualHeatVsControlAllTFsNovTruth.csv"),
  file.path(DIR_FIG, "Fig4ATFModuleDeltaHeatmap.pdf"),
  file.path(DIR_FIG, "Fig4BTFModuleNetwork.pdf"),
  file.path(OUT_DIR, "AstroM7ScoreViolinBySubcluster.pdf"),
  file.path(OUT_DIR, "AstroM7ScoreViolinBySubcluster.png")
)

for (f in files_expected) {
  cat(sprintf("  %-60s %s\n",
              basename(f),
              ifelse(file.exists(f), "EXISTS", "MISSING")))
}


################################################################################
# 14. SUMMARY
################################################################################

cat("\n========================================\n")
cat("SCRIPT 09 COMPLETE — hdWGCNA SUMMARY\n")
cat("========================================\n")
cat("WGCNA slot name          :", wgcna_name, "\n")
cat("Modules (Astro-M7 etc.)  : see AstroModuleTraitCorrelationsRNA.csv\n")
cat("TF regulon pairs         :", nrow(tf_regulons_full), "\n")
cat("Differential TFs (manual):", nrow(dregs_manual), "\n")
cat("\nOutput files:\n")
cat("  AstrohdWGCNARNAWorkspaceNov2025.RData  (pre-computed workspace)\n")
cat("  AstroModuleTraitCorrelationsRNA.csv\n")
cat("  AstroModuleTraitHeatmapRNA.pdf / .svg\n")
cat("  AstroPositiveRegulonScoresAstroNovTruth.rds\n")
cat("  AstroNegativeRegulonScoresAstroNovTruth.rds\n")
cat("  AstroDiffRegulonsManualHeatVsControlAllTFsNovTruth.csv\n")
cat("  FIGURES/Fig4ATFModuleDeltaHeatmap.pdf\n")
cat("  FIGURES/Fig4BTFModuleNetwork.pdf\n")
cat("  AstroM7ScoreViolinBySubcluster.pdf / .png\n\n")
cat("Next: Run 10_Cicero_Chromatin_Rewiring.R\n\n")

