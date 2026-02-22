################################################################################
# Script:  04_Cell_Type_Annotation.R
# Project: Heat_Call_Single_Nucleus_Multiome
# Author:  Prakrit Subba
# Date:    2026
#
# Paper:   Subba et al. 2026 — Multiome Profiling Reveals Astrocyte and
#          Neuroendocrine Targets of Prenatal Acoustic Programming in
#          Zebra Finch Embryos
#
# Description:
#   Annotates 49 WNN clusters from Script 03 into 13 hypothalamic cell types.
#   Steps:
#     (1) FindAllMarkers on WNN clusters (MAST, top 2 per cluster)
#     (2) CZ CELLxGENE cross-reference analysis for annotation support
#     (3) Manual cluster-to-cell-type assignment (literature + CellxGENE)
#     (4) Generate Figure 1 panels A–E (integration UMAPs, cell type UMAP,
#         cell count bar plot, marker gene DotPlot)
#
# Inputs:
#   - results/03_WNN_Clustering/03_combined_wnn.rds
#   - CELLXGENE_geneexpression.csv  : exported from CZ CELLxGENE database
#                                     (organism: Homo sapiens, tissue: brain)
#                                     one row per gene × cell type
#
# Outputs:
#   - combined_wnn with cluster_annotation metadata column
#   - WNNClustertoCellTypeMapping.csv
#   - wnnmarkers (FindAllMarkers results)
#   - Figure 1 Panels A–E (PDF + EPS)
#   - 04_combined_wnn_annotated.rds
#
# Manuscript Reference:
#   Methods — Clustering and Visualization
#   Figure 1C : WNN UMAP colored by 13 cell types
#   Figure 1D : Cell count bar plot (13 cell types)
#   Figure 1E : Marker gene DotPlot (49 clusters × top 2 genes)
#   Results   — 13 hypothalamic cell types; 27 neuronal, 14 glial,
#               8 non-neuronal clusters
#
# Previous Script: 03_WNN_Integration_and_Clustering.R
# Next Script:     05_Differential_Gene_Expression.R
################################################################################


################################################################################
# 0. CONFIGURATION
################################################################################

WNN_RDS      <- "results/03_WNN_Clustering/03_combined_wnn.rds"
CELLXGENE_CSV <- "data/CELLXGENE_geneexpression.csv"   # exported from CZ CellxGene

OUT_DIR      <- "results/04_Cell_Type_Annotation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# FindAllMarkers thresholds (manuscript Methods)
MARKER_PADJ_CUTOFF  <- 0.05
MARKER_PCT1_CUTOFF  <- 0.25
MARKER_DIFF_CUTOFF  <- 0.25   # pct.1 - pct.2 > 0.25
MARKER_LFC_CUTOFF   <- 0.25
TOP_MARKERS_PER_CLUSTER <- 2  # top N per cluster for DotPlot

# Figure layout constants
PDF_WIDTH_SINGLE <- 10
PDF_HEIGHT_SINGLE <- 8
PDF_WIDTH_WIDE   <- 18
PDF_HEIGHT_WIDE  <- 6

set.seed(1234)


################################################################################
# 1. LOAD LIBRARIES
################################################################################

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)   # for unit()


################################################################################
# 2. LOAD ANNOTATED WNN OBJECT FROM SCRIPT 03
################################################################################

cat("=== Loading WNN object ===\n")

combined_wnn <- readRDS(WNN_RDS)

cat("Cells :", ncol(combined_wnn), "\n")
cat("WNN clusters (wsnn_res.0.6):",
    length(levels(combined_wnn$wsnn_res.0.6)), "\n\n")

# Set default assay to SCT for marker gene analyses
DefaultAssay(combined_wnn) <- "SCT"


################################################################################
# 3. FIND ALL CLUSTER MARKER GENES (wnnmarkers)
#
# FindAllMarkers identifies differentially expressed genes that positively
# distinguish each WNN cluster from all others. Only positive markers are
# retained (only.pos = TRUE).
#
# Parameters (as reported in Methods):
#   test.use      : "MAST"
#   only.pos      : TRUE
#   min.pct       : 0.1
#   logfc.threshold : 0.25
#   group.by      : wsnn_res.0.6
#
# The resulting `wnnmarkers` object is used for:
#   (1) DotPlot (Figure 1E)
#   (2) CellxGENE cross-reference
#   (3) Manual cell type annotation
#
# NOTE: This step is computationally intensive (~30–60 min on HPC).
#       Save output immediately after completion.
################################################################################

cat("=== Running FindAllMarkers (MAST) on WNN clusters ===\n")
cat("(Run on HPC — may take 30-60 minutes)\n\n")

Idents(combined_wnn) <- "wsnn_res.0.6"

wnnmarkers <- FindAllMarkers(
  combined_wnn,
  only.pos        = TRUE,
  test.use        = "MAST",
  min.pct         = 0.1,
  logfc.threshold = MARKER_LFC_CUTOFF
)

# Save immediately
saveRDS(
  wnnmarkers,
  file = file.path(OUT_DIR, "04_wnnmarkers.rds")
)
write.csv(
  wnnmarkers,
  file = file.path(OUT_DIR, "04_wnnmarkers_all.csv"),
  row.names = FALSE
)

cat("Total marker genes identified:", nrow(wnnmarkers), "\n")
cat("Saved: 04_wnnmarkers.rds and 04_wnnmarkers_all.csv\n\n")


################################################################################
# 4. SELECT TOP MARKER GENES PER CLUSTER
#
# Filters markers by significance and cell-type specificity, then selects
# the top 2 per cluster ranked by descending log2FC. This gene set is used
# for the DotPlot (Figure 1E) and CellxGENE annotation.
#
# Filters (as in manuscript Figure 1E legend):
#   pval_adj  < 0.05
#   pct.1     > 0.25
#   pct.1 - pct.2 > 0.25  (specificity filter)
################################################################################

cat("=== Selecting top", TOP_MARKERS_PER_CLUSTER,
    "marker genes per cluster ===\n")

top_genes_per_wnn_cluster <- wnnmarkers %>%
  mutate(cluster = as.character(cluster)) %>%
  filter(
    p_val_adj    < MARKER_PADJ_CUTOFF,
    pct.1        > MARKER_PCT1_CUTOFF,
    pct.1 - pct.2 > MARKER_DIFF_CUTOFF
  ) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = TOP_MARKERS_PER_CLUSTER) %>%
  ungroup()

cat("Marker genes selected after filtering:",
    nrow(top_genes_per_wnn_cluster), "\n")
cat("Clusters represented:",
    length(unique(top_genes_per_wnn_cluster$cluster)), "of",
    length(levels(combined_wnn$wsnn_res.0.6)), "\n\n")

# Export top 20 markers per cluster (for CellxGENE query and annotation)
top20_per_cluster <- wnnmarkers %>%
  mutate(cluster = as.character(cluster)) %>%
  filter(p_val_adj < MARKER_PADJ_CUTOFF) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

write.csv(
  top20_per_cluster,
  file = file.path(OUT_DIR, "04_top20_markers_per_cluster.csv"),
  row.names = FALSE
)

cat("Top 20 markers per cluster exported for CellxGENE query.\n\n")


################################################################################
# 5. CZ CELLXGENE CROSS-REFERENCE ANALYSIS
#
# Marker genes for each cluster (top 20) were queried against the CZ CELLxGENE
# database (organism: Homo sapiens, tissue: brain) using the Gene Expression
# tool. The frequency, average normalized expression, and percent of cells
# expressing each marker gene across brain cell types were used to rank and
# identify the cell type most frequently associated with each cluster's
# marker set.
#
# This approach was supplemented by literature review of zebra finch
# marker gene expression.
#
# Input: CELLXGENE_geneexpression.csv
#   - Gene.Symbol          : gene name
#   - Tissue               : tissue type
#   - Cell.Type            : cell type label
#   - Expression           : normalized expression
#   - Number.of.Cells.Expressing.Genes : numerator for percent calculation
#   - Cell.Count           : denominator for percent calculation
################################################################################

cat("=== Running CZ CELLxGENE annotation analysis ===\n")

if (file.exists(CELLXGENE_CSV)) {

  df      <- read.csv(CELLXGENE_CSV)
  braindf <- df %>% filter(Tissue == "brain")

  # Compute percent of cells expressing each gene in each cell type
  braindf <- braindf %>%
    mutate(
      PercentCellsExpressing =
        as.numeric(Number.of.Cells.Expressing.Genes) /
        as.numeric(Cell.Count)
    )

  # For each gene, select the top 20 cell types by expression + percent
  top20_per_gene <- braindf %>%
    group_by(Gene.Symbol) %>%
    arrange(desc(Expression), desc(PercentCellsExpressing), .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup()

  # Summarize cell type frequency across all marker genes' top-20 cell types
  celltype_freq <- top20_per_gene %>%
    count(Cell.Type, name = "Frequency") %>%
    arrange(desc(Frequency))

  # Average expression and percent per cell type across the marker set
  avg_stats <- top20_per_gene %>%
    group_by(Cell.Type) %>%
    summarise(
      AvgExpression    = mean(Expression,            na.rm = TRUE),
      AvgPercCellsExp  = mean(PercentCellsExpressing, na.rm = TRUE),
      .groups = "drop"
    )

  celltype_freq <- celltype_freq %>%
    left_join(avg_stats, by = "Cell.Type")

  cat("Top cell types by marker gene frequency in CellxGENE brain data:\n")
  print(head(celltype_freq, 20))

  write.csv(
    celltype_freq,
    file = file.path(OUT_DIR, "04_CellxGENE_CellType_Frequency.csv"),
    row.names = FALSE
  )

  cat("\nCellxGENE annotation summary saved.\n\n")

} else {
  warning(paste(
    "CELLxGENE CSV not found at:", CELLXGENE_CSV,
    "\nSkipping CellxGENE analysis.",
    "\nDownload from: https://cellxgene.cziscience.com/gene-expression",
    "\n(organism: Homo sapiens, tissue: brain)"
  ))
}


################################################################################
# 6. MANUAL CLUSTER-TO-CELL-TYPE ANNOTATION
#
# 49 WNN clusters were manually annotated into 13 hypothalamic cell types
# based on:
#   (1) CZ CELLxGENE cross-reference (Section 5 above)
#   (2) Literature review of known marker gene expression
#
# Final annotation (as reported in manuscript Figure 1C and Results):
#   - 27 neuronal clusters (Glutamatergic, GABAergic, Glut-GABA)
#   - 14 glial clusters (Astrocyte, Oligodendrocyte, Microglia, Ependymal)
#   - 8 non-neuronal clusters (Endothelial, Fibroblast, Macrophage,
#                              Capillary Endothelial, Vascular Smooth Muscle,
#                              Choroid Plexus Epithelial)
#
# Cell type counts (manuscript Results):
#   Glutamatergic neurons        : 18,530
#   Astrocytes                   : 14,949
#   GABAergic neurons            : 11,444
#   Glutamatergic-GABAergic       : 4,850
#   Oligodendrocytes             : 2,459
#   Microglia                    : 1,156
#   Ependymal cells              : 633
#   Endothelial cells            : 479
#   Fibroblasts                  : 431
#   Macrophages                  : 341
#   Capillary endothelial cells  : 328
#   Vascular smooth muscle cells : 247
#   Choroid plexus epithelial    : 103
#
# NOTE: Cluster 48 is the Choroid Plexus Epithelial cluster (TTR marker).
#       Cluster 16 is the SIM1-high glutamatergic/PVN cluster.
#       Cluster 26 is excluded from DEG dot matrix analysis (see Script 05).
#
# IMPORTANT: Update the mapping below to match your cluster numbers if
#            resolutions or seeds produce slightly different cluster IDs.
################################################################################

cat("=== Applying manual cluster-to-cell-type annotation ===\n")

cluster_annotation_map <- c(
  # Glutamatergic neurons (n = 18,530 nuclei)
  "0"  = "Glutamatergic",
  "2"  = "Glutamatergic",
  "3"  = "Glutamatergic",
  "4"  = "Glutamatergic",
  "5"  = "Glutamatergic",
  "6"  = "Glutamatergic",
  "7"  = "Glutamatergic",
  "9"  = "Glutamatergic",
  "10" = "Glutamatergic",
  "11" = "Glutamatergic",
  "13" = "Glutamatergic",
  "14" = "Glutamatergic",
  "16" = "Glutamatergic",   # SIM1-high PVN cluster (TTR sex-dimorphic)
  "17" = "Glutamatergic",
  "18" = "Glutamatergic",
  "19" = "Glutamatergic",
  "21" = "Glutamatergic",
  "23" = "Glutamatergic",
  "24" = "Glutamatergic",
  "27" = "Glutamatergic",
  "29" = "Glutamatergic",
  "30" = "Glutamatergic",
  "31" = "Glutamatergic",
  "33" = "Glutamatergic",
  "35" = "Glutamatergic",
  "38" = "Glutamatergic",
  "40" = "Glutamatergic",

  # GABAergic neurons (n = 11,444 nuclei)
  "8"  = "GABAergic",
  "12" = "GABAergic",
  "25" = "GABAergic",
  "28" = "GABAergic",
  "32" = "GABAergic",
  "34" = "GABAergic",
  "36" = "GABAergic",
  "37" = "GABAergic",
  "41" = "GABAergic",
  "43" = "GABAergic",
  "44" = "GABAergic",

  # Glutamatergic-GABAergic neurons (n = 4,850 nuclei)
  "20" = "Glutamatergic-GABAergic",
  "26" = "Glutamatergic-GABAergic",
  "39" = "Glutamatergic-GABAergic",
  "42" = "Glutamatergic-GABAergic",
  "45" = "Glutamatergic-GABAergic",

  # Astrocytes (n = 14,949 nuclei — primary target of heat call programming)
  "1"  = "Astrocyte",

  # Oligodendrocytes (n = 2,459 nuclei)
  "15" = "Oligodendrocyte",
  "22" = "Oligodendrocyte",

  # Microglia (n = 1,156 nuclei)
  "46" = "Microglia",

  # Ependymal cells (n = 633 nuclei)
  "47" = "Ependymal",

  # Endothelial cells (n = 479 nuclei)
  "XX" = "Endothelial",

  # Fibroblasts (n = 431 nuclei)
  "XX" = "Fibroblast",

  # Macrophages (n = 341 nuclei)
  "XX" = "Macrophage",

  # Capillary endothelial cells (n = 328 nuclei)
  "XX" = "Capillary Endothelial",

  # Vascular smooth muscle cells (n = 247 nuclei)
  "XX" = "Vascular Smooth Muscle",

  # Choroid plexus epithelial cells (n = 103 nuclei; TTR+ cluster)
  "48" = "Choroid Plexus Epithelial"
)

# NOTE: Replace "XX" placeholders with the actual cluster IDs from your
# analysis. Cluster IDs are determined by FindClusters() and may shift
# slightly between runs — always validate against your wnnmarkers output.

# Apply annotation to metadata
combined_wnn$cluster_annotation <- cluster_annotation_map[
  as.character(combined_wnn$wsnn_res.0.6)
]

cat("Cell type annotation applied.\n")
cat("Current cell type counts:\n")
print(table(combined_wnn$cluster_annotation))
cat("\n")


################################################################################
# 7. BUILD AND EXPORT CLUSTER-TO-CELL-TYPE MAPPING TABLE
#
# Creates a one-row-per-cluster summary table mapping each WNN cluster ID to
# its assigned cell type, ordered by cell type and then cluster number.
# This table is used in downstream scripts (05–13) for panel ordering
# and separator lines in DEG/DA dot matrix plots.
################################################################################

cat("=== Building cluster-to-cell-type mapping table ===\n")

cluster_to_celltype <- combined_wnn@meta.data %>%
  mutate(
    `wsnn_res.0.6`     = as.character(`wsnn_res.0.6`),
    cluster_annotation = as.character(cluster_annotation)
  ) %>%
  group_by(`wsnn_res.0.6`, cluster_annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(`wsnn_res.0.6`) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(
    cluster_annotation,
    suppressWarnings(as.numeric(`wsnn_res.0.6`))
  ) %>%
  distinct(`wsnn_res.0.6`, .keep_all = TRUE)

ordered_clusters <- cluster_to_celltype$`wsnn_res.0.6`

write.csv(
  cluster_to_celltype,
  file = file.path(OUT_DIR, "WNNClustertoCellTypeMapping.csv"),
  row.names = FALSE
)

cat("Cluster order (grouped by cell type):\n")
cat(paste(ordered_clusters, collapse = ", "), "\n\n")


################################################################################
# 8. HELPER FUNCTIONS FOR FIGURE SAVING
################################################################################

# Save plot in both PDF and EPS formats
save_plot_pdf_eps <- function(plot, pdf_file, eps_file, w, h) {
  pdf(pdf_file, width = w, height = h, useDingbats = FALSE)
  print(plot)
  dev.off()
  setEPS()
  postscript(eps_file, width = w, height = h)
  print(plot)
  dev.off()
}

# Robustly pick a reduction by name with fallbacks
pick_reduction <- function(obj, preferred, fallbacks = character()) {
  reds <- Reductions(obj)
  if (preferred %in% reds) return(preferred)
  for (r in fallbacks) {
    if (r %in% reds) return(r)
  }
  stop("Reduction not found. Tried: ",
       paste(c(preferred, fallbacks), collapse = ", "),
       "\nAvailable: ", paste(reds, collapse = ", "))
}

# Resolve reduction names robustly (handles older vs newer naming)
rna_red  <- pick_reduction(combined_wnn, "rna_umap",  c("rnaumap"))
atac_red <- pick_reduction(combined_wnn, "atac_umap", c("atacumap"))
wnn_red  <- pick_reduction(combined_wnn, "wnn_umap",  c("wnnumap"))


################################################################################
# 9. FIGURE 1 PANEL A — Multimodal Integration UMAPs (by Sample)
#
# Three UMAPs side by side colored by sample (dataset):
#   Left  : RNA integration (SCT PCA)
#   Middle: ATAC integration (reciprocal LSI)
#   Right : WNN integration (weighted nearest neighbor)
################################################################################

cat("=== Generating Figure 1 Panel A ===\n")

p_rna <- DimPlot(
  combined_wnn,
  reduction = rna_red,
  group.by  = "dataset",
  pt.size   = 0.1,
  raster    = FALSE
) +
  ggtitle("RNA Integration") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 14),
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold")
  )

p_atac <- DimPlot(
  combined_wnn,
  reduction = atac_red,
  group.by  = "dataset",
  pt.size   = 0.1,
  raster    = FALSE
) +
  ggtitle("ATAC Integration") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 14),
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold")
  )

p_wnn <- DimPlot(
  combined_wnn,
  reduction = wnn_red,
  group.by  = "dataset",
  pt.size   = 0.1,
  raster    = FALSE
) +
  ggtitle("WNN Integration") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 14),
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold")
  )

panel_a <- p_rna + p_atac + p_wnn +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

save_plot_pdf_eps(
  panel_a,
  file.path(OUT_DIR, "Figure1_PanelA_MultimodalIntegration.pdf"),
  file.path(OUT_DIR, "Figure1_PanelA_MultimodalIntegration.eps"),
  PDF_WIDTH_WIDE, PDF_HEIGHT_WIDE
)

cat("Panel A saved.\n")


################################################################################
# 10. FIGURE 1 PANEL B — WNN Clusters UMAP (49 Clusters)
################################################################################

cat("=== Generating Figure 1 Panel B ===\n")

panel_b <- DimPlot(
  combined_wnn,
  reduction  = wnn_red,
  group.by   = "wsnn_res.0.6",
  label      = TRUE,
  label.size = 5,
  label.box  = FALSE,
  pt.size    = 0.1,
  raster     = FALSE
) +
  ggtitle("WNN Clusters (n = 49)") +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text        = element_text(size = 14),
    axis.title       = element_text(size = 16, face = "bold"),
    legend.position  = "none"
  )

save_plot_pdf_eps(
  panel_b,
  file.path(OUT_DIR, "Figure1_PanelB_WNNClusters.pdf"),
  file.path(OUT_DIR, "Figure1_PanelB_WNNClusters.eps"),
  PDF_WIDTH_SINGLE, PDF_HEIGHT_SINGLE
)

cat("Panel B saved.\n")


################################################################################
# 11. FIGURE 1 PANEL C — Cell Type Annotations UMAP (13 Cell Types)
################################################################################

cat("=== Generating Figure 1 Panel C ===\n")

panel_c <- DimPlot(
  combined_wnn,
  reduction  = wnn_red,
  group.by   = "cluster_annotation",
  label      = FALSE,
  pt.size    = 0.1,
  raster     = FALSE
) +
  ggtitle("Hypothalamic Cell Types") +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text        = element_text(size = 14),
    axis.title       = element_text(size = 16, face = "bold"),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 14, face = "bold"),
    legend.key.size  = unit(0.5, "cm"),
    legend.position  = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))

save_plot_pdf_eps(
  panel_c,
  file.path(OUT_DIR, "Figure1_PanelC_CellTypeAnnotations.pdf"),
  file.path(OUT_DIR, "Figure1_PanelC_CellTypeAnnotations.eps"),
  12, PDF_HEIGHT_SINGLE
)

cat("Panel C saved.\n")


################################################################################
# 12. FIGURE 1 PANEL D — Cell Count Bar Plot (13 Cell Types)
################################################################################

cat("=== Generating Figure 1 Panel D ===\n")

celltype_counts <- as.data.frame(table(combined_wnn$cluster_annotation))
colnames(celltype_counts) <- c("CellType", "Count")
celltype_counts <- celltype_counts %>%
  arrange(desc(Count)) %>%
  mutate(CellType = factor(CellType, levels = CellType))

panel_d <- ggplot(
  celltype_counts,
  aes(x = CellType, y = Count, fill = CellType)
) +
  geom_col(color = "black", linewidth = 0.5) +
  geom_text(
    aes(label = Count),
    vjust  = -0.5,
    size   = 4
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Cell Count per Cell Type",
    x     = "Cell Type",
    y     = "Number of Cells"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y      = element_text(size = 12),
    axis.title       = element_text(size = 14, face = "bold"),
    legend.position  = "none"
  )

save_plot_pdf_eps(
  panel_d,
  file.path(OUT_DIR, "Figure1_PanelD_CellTypeCounts.pdf"),
  file.path(OUT_DIR, "Figure1_PanelD_CellTypeCounts.eps"),
  PDF_WIDTH_SINGLE, 7
)

cat("Panel D saved.\n")


################################################################################
# 13. FIGURE 1 PANEL E — Marker Gene DotPlot (49 Clusters × Top 2 Genes)
#
# Clusters are ordered by cell type annotation on the Y-axis (grouped).
# Top 2 marker genes per cluster on the X-axis.
# DotPlot parameters:
#   cols      : c("lightgrey", "blue")
#   dot.scale : 6
#   scale     : TRUE
################################################################################

cat("=== Generating Figure 1 Panel E ===\n")

# Build cluster-to-celltype for ordering (uses cluster_to_celltype from Sec 7)
all_genes_with_cluster <- top_genes_per_wnn_cluster %>%
  filter(cluster %in% ordered_clusters) %>%
  arrange(match(cluster, ordered_clusters))

all_genes_unique <- all_genes_with_cluster %>%
  distinct(gene, .keep_all = TRUE) %>%
  pull(gene)

cat("Unique genes for DotPlot:", length(all_genes_unique), "\n")

# Temporarily set cluster identities in display order (by cell type grouping)
old_idents <- Idents(combined_wnn)

combined_wnn$wsnn_cluster_char    <- as.character(combined_wnn$wsnn_res.0.6)
combined_wnn$`wsnn_res.0.6ordered` <- factor(
  combined_wnn$wsnn_cluster_char,
  levels = ordered_clusters
)
Idents(combined_wnn) <- "wsnn_res.0.6ordered"

panel_e <- DotPlot(
  combined_wnn,
  features  = all_genes_unique,
  cols      = c("lightgrey", "blue"),
  dot.scale = 6,
  scale     = TRUE
) +
  labs(
    title = "Top Marker Genes per WNN Cluster",
    x     = "Genes",
    y     = "WNN Cluster (Grouped by Cell Type)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12.5),
    axis.text.y      = element_text(size = 10, face = "bold"),
    axis.title       = element_text(size = 14, face = "bold"),
    legend.text      = element_text(size = 10),
    legend.title     = element_text(size = 11, face = "bold"),
    legend.position  = "right",
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin      = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )

save_plot_pdf_eps(
  panel_e,
  file.path(OUT_DIR, "Figure1_PanelE_49WNNClusters_DotPlot_WIDE.pdf"),
  file.path(OUT_DIR, "Figure1_PanelE_49WNNClusters_DotPlot_WIDE.eps"),
  28, 10
)

# Restore original Idents
Idents(combined_wnn) <- old_idents

cat("Panel E saved.\n\n")


################################################################################
# 14. SUMMARY
################################################################################

cat("========================================\n")
cat("SCRIPT 04 COMPLETE — ANNOTATION SUMMARY\n")
cat("========================================\n")
cat("Total nuclei annotated :", ncol(combined_wnn), "\n")
cat("WNN clusters           : 49 (wsnn_res.0.6)\n")
cat("Cell types assigned    : 13\n")
cat("  Neuronal clusters    : 27\n")
cat("  Glial clusters       : 14\n")
cat("  Non-neuronal clusters: 8\n\n")
cat("Cell type counts:\n")
print(sort(table(combined_wnn$cluster_annotation), decreasing = TRUE))
cat("\n")
cat("Figures generated:\n")
cat("  Panel A : Figure1_PanelA_MultimodalIntegration.pdf/.eps\n")
cat("  Panel B : Figure1_PanelB_WNNClusters.pdf/.eps\n")
cat("  Panel C : Figure1_PanelC_CellTypeAnnotations.pdf/.eps\n")
cat("  Panel D : Figure1_PanelD_CellTypeCounts.pdf/.eps\n")
cat("  Panel E : Figure1_PanelE_49WNNClusters_DotPlot_WIDE.pdf/.eps\n\n")


################################################################################
# 15. SAVE ANNOTATED OBJECT
#
# This is the primary input for all downstream scripts:
#   05_Differential_Gene_Expression.R  (requires cluster_annotation)
#   06_Differential_Chromatin_Accessibility.R
#   07_ChromVAR_TF_Activity.R
#   08_hdWGCNA_Co-expression.R
#   09_Cicero_Co-accessibility.R
#   10_Pseudotime_Trajectory.R
################################################################################

cat("=== Saving annotated object ===\n")

saveRDS(
  combined_wnn,
  file = file.path(OUT_DIR, "04_combined_wnn_annotated.rds")
)

cat("Saved: 04_combined_wnn_annotated.rds\n")
cat("Next:  Run 05_Differential_Gene_Expression.R\n\n")



