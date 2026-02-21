################################################################################
# Script:  03_WNN_Integration_and_Clustering.R
# Project: Heat_Call_Single_Nucleus_Multiome
# Author:  Prakrit Subba
# Date:    2026
#
# Paper:   Subba et al. 2026 — Multiome Profiling Reveals Astrocyte and
#          Neuroendocrine Targets of Prenatal Acoustic Programming in
#          Zebra Finch Embryos
#
# Description:
#   Integrates the QC-filtered, normalized snRNA-seq and snATAC-seq objects
#   from Scripts 01 and 02 into a joint multiome object. Performs:
#     (1) RNA-only clustering: PCA → UMAP → Leiden (47 clusters, res = 1.1)
#     (2) ATAC-only integration via reciprocal LSI (RLSI) and clustering
#         (57 clusters, res = 1.3)
#     (3) Weighted Nearest Neighbor (WNN) joint clustering across both
#         modalities (49 clusters, res = 0.6)
#
# Inputs:
#   - results/01_RNA_QC/01_merged_seurat_filt.rds
#       SCTransform-normalized RNA Seurat object (from Script 01)
#   - results/02_ATAC_QC/02_combined_atac_filt.rds
#       TF-IDF normalized, LSI-reduced ATAC Signac object (from Script 02)
#
# Outputs:
#   - combined_wnn : Joint RNA+ATAC Seurat object with WNN clustering
#   - 03_combined_wnn.rds : Saved for all downstream scripts (04–13)
#   - UMAP plots: RNA-only, ATAC-only (integrated), WNN
#
# Manuscript Reference:
#   Methods — WNN Integration and Clustering
#   Figure 1B : RNA UMAP (47 clusters), ATAC UMAP (57 clusters)
#   Figure 1C : WNN UMAP (49 clusters)
#   Results   — 55,950 nuclei analyzed across both modalities
#
# Previous Scripts: 01_RNA_QC_and_Processing.R, 02_ATAC_QC_and_Processing.R
# Next Script:      04_Cell_Type_Annotation.R
################################################################################


################################################################################
# 0. CONFIGURATION
################################################################################

RNA_RDS  <- "results/01_RNA_QC/01_merged_seurat_filt.rds"
ATAC_RDS <- "results/02_ATAC_QC/02_combined_atac_filt.rds"

OUT_DIR  <- "results/03_WNN_Clustering"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# RNA clustering parameters (manuscript: 47 clusters)
RNA_PCA_DIMS       <- 1:35
RNA_CLUSTER_RES    <- 1.1

# ATAC integration and clustering parameters (manuscript: 57 clusters)
ATAC_LSI_DIMS      <- 2:50          # LSI dim 1 excluded (depth-correlated)
ATAC_CLUSTER_RES   <- 1.3
ATAC_CLUSTER_ALG   <- 3             # Leiden algorithm

# WNN clustering parameters (manuscript: 49 clusters)
WNN_CLUSTER_RES    <- 0.6
WNN_CLUSTER_ALG    <- 3             # Leiden algorithm

set.seed(1234)


################################################################################
# 1. LOAD LIBRARIES
################################################################################

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(patchwork)


################################################################################
# 2. LOAD PROCESSED OBJECTS FROM SCRIPTS 01 AND 02
################################################################################

cat("=== Loading processed objects ===\n")

merged_seurat_filt <- readRDS(RNA_RDS)
combined_filt      <- readRDS(ATAC_RDS)

cat("RNA object  — cells:", ncol(merged_seurat_filt),
    "| assays:", paste(names(merged_seurat_filt@assays), collapse = ", "), "\n")
cat("ATAC object — cells:", ncol(combined_filt),
    "| reductions:", paste(names(combined_filt@reductions), collapse = ", "), "\n\n")


################################################################################
# 3. RNA-ONLY CLUSTERING
#
# PCA is run on the SCTransform-normalized assay. FindNeighbors and
# FindClusters identify 47 clusters at resolution = 1.1 (as reported in
# manuscript Figure 1B). RunUMAP projects cells into 2D for visualization.
################################################################################

cat("=== RNA-only: PCA, clustering, UMAP ===\n")

# PCA on SCT assay
merged_seurat_filt <- RunPCA(merged_seurat_filt, assay = "SCT")

# Graph-based clustering
merged_seurat_filt <- FindNeighbors(merged_seurat_filt, dims = RNA_PCA_DIMS)
merged_seurat_filt <- FindClusters(merged_seurat_filt, resolution = RNA_CLUSTER_RES)

# UMAP for visualization
merged_seurat_filt <- RunUMAP(
  merged_seurat_filt,
  dims           = RNA_PCA_DIMS,
  reduction.name = "rna_umap",
  reduction.key  = "RNAUMAP_"
)

cat("RNA clusters (res =", RNA_CLUSTER_RES, "):",
    length(unique(Idents(merged_seurat_filt))), "\n\n")


################################################################################
# 4. ATAC-ONLY UMAP (PRE-INTEGRATION)
#
# An initial UMAP is computed directly on the LSI reduction for visualization
# purposes before cross-sample integration. Component 1 is excluded because
# it correlates with sequencing depth rather than biological variation.
################################################################################

cat("=== ATAC-only: initial LSI UMAP ===\n")

combined_filt <- RunUMAP(
  combined_filt,
  dims           = ATAC_LSI_DIMS,
  reduction      = "lsi",
  reduction.name = "atac_umap_preintegration",
  reduction.key  = "ATACUMAP_"
)

cat("Initial ATAC UMAP computed (pre-integration).\n\n")


################################################################################
# 5. ATAC CROSS-SAMPLE INTEGRATION (RECIPROCAL LSI)
#
# Because ATAC data can have stronger batch effects than RNA, samples are
# integrated using Reciprocal LSI (RLSI) — the ATAC equivalent of Reciprocal
# PCA integration. Steps:
#   (1) SplitObject() by sample (dataset)
#   (2) FindIntegrationAnchors() with reduction = "rlsi", dims = 2:50
#   (3) IntegrateEmbeddings() produces a harmonized "integrated_lsi" reduction
#   (4) RunUMAP / FindNeighbors / FindClusters on integrated_lsi
#
# Note: anchor.features = rownames(combined_filt) uses all consensus peaks
# as anchoring features.
################################################################################

cat("=== ATAC cross-sample integration (RLSI) ===\n")
cat("  Splitting object by sample...\n")

object.list <- SplitObject(combined_filt, split.by = "dataset")

cat("  Finding integration anchors (reduction = 'rlsi', dims = 2:50)...\n")
cat("  (This may take several minutes)\n")

integration.anchors <- FindIntegrationAnchors(
  object.list     = object.list,
  anchor.features = rownames(combined_filt),
  reduction       = "rlsi",
  dims            = ATAC_LSI_DIMS
)

cat("  Integrating embeddings...\n")

integrated <- IntegrateEmbeddings(
  anchorset        = integration.anchors,
  reductions       = combined_filt[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate  = 1:50
)

# UMAP on integrated LSI
integrated <- RunUMAP(
  integrated,
  reduction      = "integrated_lsi",
  dims           = ATAC_LSI_DIMS,
  reduction.name = "atac_umap",
  reduction.key  = "ATACUMAP_"
)

# Graph construction and clustering on integrated LSI
integrated <- FindNeighbors(
  integrated,
  reduction = "integrated_lsi",
  dims      = ATAC_LSI_DIMS
)

integrated <- FindClusters(
  integrated,
  verbose     = FALSE,
  algorithm   = ATAC_CLUSTER_ALG,
  resolution  = ATAC_CLUSTER_RES
)

cat("ATAC clusters (res =", ATAC_CLUSTER_RES, "):",
    length(unique(Idents(integrated))), "\n\n")


################################################################################
# 6. ASSEMBLE JOINT RNA + ATAC MULTIOME OBJECT
#
# Cells present in both the RNA-filtered and ATAC-filtered/integrated objects
# are retained. The ATAC ChromatinAssay and the integrated_lsi reduction are
# transferred into the RNA Seurat object to create a single joint object.
# This is the standard Signac 10X Multiome workflow for objects processed
# independently per modality.
#
# Manuscript: 55,950 nuclei retained across both modalities
################################################################################

cat("=== Assembling joint RNA + ATAC object ===\n")

# Identify barcodes present in both filtered modalities
common_cells <- intersect(
  colnames(merged_seurat_filt),
  colnames(integrated)
)

cat("  RNA cells (post-QC) :", ncol(merged_seurat_filt), "\n")
cat("  ATAC cells (post-QC):", ncol(integrated), "\n")
cat("  Common cells        :", length(common_cells), "\n\n")

# Subset both objects to shared barcodes
rna_common  <- merged_seurat_filt[, common_cells]
atac_common <- integrated[, common_cells]

# Transfer ATAC ChromatinAssay into RNA object
rna_common[["ATAC"]] <- atac_common[["ATAC"]]

# Transfer integrated_lsi reduction (required for WNN)
rna_common[["lsi"]] <- atac_common[["integrated_lsi"]]

# Transfer ATAC UMAP for visualization
rna_common[["atac_umap"]] <- atac_common[["atac_umap"]]

# Propagate ATAC metadata columns
rna_common$dataset           <- atac_common$dataset[common_cells]
rna_common$nCount_ATAC       <- atac_common$nCount_ATAC[common_cells]
rna_common$pct_reads_in_peaks <- atac_common$pct_reads_in_peaks[common_cells]
rna_common$nucleosome_signal  <- atac_common$nucleosome_signal[common_cells]
rna_common$TSS.enrichment     <- atac_common$TSS.enrichment[common_cells]

# Rename to combined_wnn (working name used throughout downstream scripts)
combined_wnn <- rna_common

cat("Joint object assembled.\n")
cat("Assays  :", paste(names(combined_wnn@assays), collapse = ", "), "\n")
cat("Reductions:", paste(names(combined_wnn@reductions), collapse = ", "), "\n\n")


################################################################################
# 7. WNN JOINT CLUSTERING
#
# Weighted Nearest Neighbor (WNN) analysis (Hao et al. 2021, Cell) integrates
# RNA and ATAC modalities by learning per-cell modality weights. Cells where
# ATAC provides more information receive a higher ATAC weight, and vice versa.
#
# Parameters:
#   reduction.list        : list("pca", "lsi")
#   dims.list             : list(1:35, 2:50)   — RNA PCA dims; ATAC LSI dims
#   modality.weight.name  : "RNA.weight"
#   resolution            : 0.6
#   algorithm             : 3 (Leiden)
#
# Result: 49 WNN clusters (wsnn_res.0.6), as reported in manuscript Figure 1C
################################################################################

cat("=== Running Weighted Nearest Neighbor (WNN) integration ===\n")

combined_wnn <- FindMultiModalNeighbors(
  object             = combined_wnn,
  reduction.list     = list("pca", "lsi"),
  dims.list          = list(RNA_PCA_DIMS, ATAC_LSI_DIMS),
  modality.weight.name = "RNA.weight",
  verbose            = TRUE
)

cat("\n=== WNN UMAP ===\n")

combined_wnn <- RunUMAP(
  object         = combined_wnn,
  nn.name        = "weighted.nn",
  reduction.name = "wnn_umap",
  reduction.key  = "WNNUMAP_"
)

cat("\n=== WNN clustering (res =", WNN_CLUSTER_RES, ", algorithm =",
    WNN_CLUSTER_ALG, ") ===\n")

combined_wnn <- FindClusters(
  object     = combined_wnn,
  resolution = WNN_CLUSTER_RES,
  graph.name = "wsnn",
  algorithm  = WNN_CLUSTER_ALG
)

cat("WNN clusters (res =", WNN_CLUSTER_RES, "):",
    length(unique(combined_wnn$wsnn_res.0.6)), "\n\n")


################################################################################
# 8. CLUSTERING SUMMARY VISUALIZATION
################################################################################

cat("=== Generating clustering UMAP plots ===\n")

# RNA-only UMAP
p_rna <- DimPlot(
  combined_wnn,
  reduction = "rna_umap",
  group.by  = "seurat_clusters",
  label     = TRUE,
  label.size = 3,
  pt.size   = 0.1,
  raster    = FALSE
) +
  ggtitle("RNA-only Clustering (47 clusters, res = 1.1)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

# ATAC-only (integrated LSI) UMAP
p_atac <- DimPlot(
  combined_wnn,
  reduction = "atac_umap",
  group.by  = "dataset",
  pt.size   = 0.1,
  raster    = FALSE
) +
  ggtitle("ATAC Integration by Sample") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

# WNN UMAP by cluster
p_wnn_cluster <- DimPlot(
  combined_wnn,
  reduction = "wnn_umap",
  group.by  = "wsnn_res.0.6",
  label     = TRUE,
  label.size = 3,
  pt.size   = 0.1,
  raster    = FALSE
) +
  ggtitle("WNN Clustering (49 clusters, res = 0.6)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

# WNN UMAP split by condition
p_wnn_condition <- DimPlot(
  combined_wnn,
  reduction  = "wnn_umap",
  group.by   = "wsnn_res.0.6",
  split.by   = "Condition",
  pt.size    = 0.1,
  label      = FALSE,
  raster     = FALSE
) +
  ggtitle("WNN UMAP — Heat Call vs Control") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

# WNN UMAP colored by modality weights
p_rna_weight <- FeaturePlot(
  combined_wnn,
  features  = "RNA.weight",
  reduction = "wnn_umap",
  pt.size   = 0.1,
  raster    = FALSE
) +
  ggtitle("RNA Modality Weight per Cell") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

# Save composite panel
p_composite <- (p_rna | p_atac) / (p_wnn_cluster | p_wnn_condition | p_rna_weight)

ggsave(
  filename = file.path(OUT_DIR, "03_Clustering_UMAP_Panel.pdf"),
  plot     = p_composite,
  width    = 22, height   = 14
)

ggsave(
  filename = file.path(OUT_DIR, "03_WNN_UMAP_bycluster.pdf"),
  plot     = p_wnn_cluster,
  width    = 10, height   = 8
)

ggsave(
  filename = file.path(OUT_DIR, "03_WNN_UMAP_byCondition.pdf"),
  plot     = p_wnn_condition,
  width    = 14, height   = 7
)

cat("UMAP plots saved to:", OUT_DIR, "\n\n")


################################################################################
# 9. CLUSTER CELL COUNT SUMMARY
################################################################################

cat("=== Cluster cell count summary ===\n")

cluster_summary <- combined_wnn@meta.data %>%
  group_by(wsnn_res.0.6, Condition) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(as.numeric(as.character(wsnn_res.0.6)))

write.csv(
  cluster_summary,
  file.path(OUT_DIR, "03_WNN_cluster_cell_counts.csv"),
  row.names = FALSE
)

cat("WNN clusters × Condition cell counts:\n")
print(
  tidyr::pivot_wider(cluster_summary,
    names_from  = Condition,
    values_from = n_cells,
    values_fill = 0
  )
)
cat("\n")


################################################################################
# 10. SUMMARY STATISTICS
################################################################################

cat("========================================\n")
cat("SCRIPT 03 COMPLETE — CLUSTERING SUMMARY\n")
cat("========================================\n")
cat("Total nuclei in joint object :", ncol(combined_wnn), "\n")
cat("RNA-only clusters  (res=1.1) :", length(unique(combined_wnn$RNA_snn_res.1.1)), "\n")
cat("WNN clusters       (res=0.6) :", length(unique(combined_wnn$wsnn_res.0.6)), "\n")
cat("\nAssays in combined_wnn :", paste(names(combined_wnn@assays), collapse = ", "), "\n")
cat("Reductions          :", paste(names(combined_wnn@reductions), collapse = ", "), "\n")
cat("\nCondition breakdown:\n")
print(table(combined_wnn$Condition, combined_wnn$Sex))
cat("\n")


################################################################################
# 11. SAVE OUTPUT
#
# combined_wnn is the central object used by ALL downstream scripts (04–13).
# It contains:
#   Assays      : RNA (counts/data), SCT (Pearson residuals), ATAC (counts)
#   Reductions  : pca, rna_umap, lsi (integrated), atac_umap, wnn_umap
#   Graphs      : RNA_snn, ATAC_nn, wsnn (WNN graph)
#   Metadata    : Condition, Sex, dataset, QC metrics, wsnn_res.0.6
################################################################################

cat("=== Saving combined_wnn ===\n")

saveRDS(
  combined_wnn,
  file = file.path(OUT_DIR, "03_combined_wnn.rds")
)

cat("Saved: 03_combined_wnn.rds\n")
cat("Next:  Run 04_Cell_Type_Annotation.R\n\n")


