################################################################################
# Script:  04_Cell_Type_Annotation.R
# Project: Heat_Call_Single_Nucleus_Multiome
# Author:  Prakrit Subba
# Date:    2026
#
# Description:
#   Annotates 49 WNN clusters into 13 hypothalamic cell types using a
#   pre-curated cluster annotation table (WNN_Cell_Annotation_Post_Processing
#   spreadsheet, included in repository). Runs FindAllMarkers for marker gene
#   identification and generates Figure 1 panels A–E.
#
# Inputs:
#   - results/03_WNN_Clustering/03_combined_wnn.rds
#   - data/WNN_Cell_Annotation_Post_Processing.csv
#       Columns: Cell.Cluster, Cluster (broad cell type), Cell.Type (CellxGene),
#                Frequency, Avg_Expression, Avg_Perc_Cells_Exp
#
# Outputs:
#   - combined_wnn with cluster_annotation metadata column
#   - 04_wnnmarkers.rds / 04_wnnmarkers_all.csv
#   - Figure 1 Panels A–E (PDF + EPS)
#   - 04_combined_wnn_annotated.rds
#
# Previous Script: 03_WNN_Integration_and_Clustering.R
# Next Script:     05_Differential_Gene_Expression.R
################################################################################


################################################################################
# 0. CONFIGURATION
################################################################################

WNN_RDS        <- "results/03_WNN_Clustering/03_combined_wnn.rds"
ANNOTATION_CSV <- "data/WNN_Cell_Annotation_Post_Processing.csv"

OUT_DIR <- "results/04_Cell_Type_Annotation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# FindAllMarkers thresholds (manuscript Methods)
MARKER_PADJ_CUTOFF       <- 0.05
MARKER_PCT1_CUTOFF       <- 0.25
MARKER_DIFF_CUTOFF       <- 0.25    # pct.1 - pct.2
MARKER_LFC_CUTOFF        <- 0.25
TOP_MARKERS_PER_CLUSTER  <- 2

# Figure layout
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
library(grid)


################################################################################
# 2. LOAD OBJECTS
################################################################################

cat("=== Loading objects ===\n")

combined_wnn    <- readRDS(WNN_RDS)
annotation_df   <- read.csv(ANNOTATION_CSV)

cat("WNN object — cells:", ncol(combined_wnn), "\n")
cat("Annotation table — clusters:", nrow(annotation_df), "\n\n")

# Expected columns in annotation CSV
stopifnot(all(c("Cell.Cluster", "Cluster") %in% colnames(annotation_df)))

DefaultAssay(combined_wnn) <- "SCT"
Idents(combined_wnn)       <- "wsnn_res.0.6"


################################################################################
# 3. APPLY CLUSTER ANNOTATIONS
#
# Maps each WNN cluster (0–48) to a broad cell type using the pre-curated
# annotation spreadsheet (WNN_Cell_Annotation_Post_Processing.csv).
# See spreadsheet for CellxGENE cross-reference details and marker gene
# evidence used to assign each cluster.
################################################################################

cat("=== Applying cluster annotations from spreadsheet ===\n")

annotation_map <- setNames(
  as.character(annotation_df$Cluster),
  as.character(annotation_df$Cell.Cluster)
)

combined_wnn$cluster_annotation <- annotation_map[
  as.character(combined_wnn$wsnn_res.0.6)
]

cat("Cell type counts:\n")
print(sort(table(combined_wnn$cluster_annotation), decreasing = TRUE))
cat("\n")


################################################################################
# 4. FIND ALL CLUSTER MARKER GENES
#
# Parameters (manuscript Methods):
#   test.use        : "MAST"
#   only.pos        : TRUE
#   min.pct         : 0.1
#   logfc.threshold : 0.25
#
# NOTE: Computationally intensive — run on HPC (~30–60 min).
#       Output saved immediately to avoid re-running.
################################################################################

cat("=== Running FindAllMarkers (MAST) ===\n")
cat("(Run on HPC — may take 30–60 minutes)\n\n")

wnnmarkers <- FindAllMarkers(
  combined_wnn,
  only.pos        = TRUE,
  test.use        = "MAST",
  min.pct         = 0.1,
  logfc.threshold = MARKER_LFC_CUTOFF
)

saveRDS(wnnmarkers, file.path(OUT_DIR, "04_wnnmarkers.rds"))
write.csv(wnnmarkers, file.path(OUT_DIR, "04_wnnmarkers_all.csv"),
          row.names = FALSE)

cat("Total markers found:", nrow(wnnmarkers), "\n\n")


################################################################################
# 5. SELECT TOP MARKER GENES PER CLUSTER (FOR DOTPLOT)
#
# Filters: pval_adj < 0.05, pct.1 > 0.25, pct.1 - pct.2 > 0.25
# Ranking: top 2 per cluster by descending avg_log2FC
################################################################################

cat("=== Selecting top", TOP_MARKERS_PER_CLUSTER, "markers per cluster ===\n")

top_genes_per_wnn_cluster <- wnnmarkers %>%
  mutate(cluster = as.character(cluster)) %>%
  filter(
    p_val_adj     < MARKER_PADJ_CUTOFF,
    pct.1         > MARKER_PCT1_CUTOFF,
    pct.1 - pct.2 > MARKER_DIFF_CUTOFF
  ) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = TOP_MARKERS_PER_CLUSTER) %>%
  ungroup()

cat("Marker genes for DotPlot:", nrow(top_genes_per_wnn_cluster), "\n\n")


################################################################################
# 6. BUILD CLUSTER ORDERING TABLE
#
# Orders clusters by cell type (alphabetical) then by cluster number.
# This ordering is used for the DotPlot Y-axis and DEG/DA dot matrix
# separators in downstream scripts (05, 06).
################################################################################

cat("=== Building cluster ordering ===\n")

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
  file.path(OUT_DIR, "WNNClustertoCellTypeMapping.csv"),
  row.names = FALSE
)

cat("Cluster order (by cell type):\n")
cat(paste(ordered_clusters, collapse = ", "), "\n\n")


################################################################################
# 7. HELPER: SAVE PDF + EPS
################################################################################

save_plot_pdf_eps <- function(plot, pdf_file, eps_file, w, h) {
  pdf(pdf_file, width = w, height = h, useDingbats = FALSE)
  print(plot); dev.off()
  setEPS()
  postscript(eps_file, width = w, height = h)
  print(plot); dev.off()
}

pick_reduction <- function(obj, preferred, fallbacks = character()) {
  reds <- Reductions(obj)
  if (preferred %in% reds) return(preferred)
  for (r in fallbacks) if (r %in% reds) return(r)
  stop("Reduction not found. Tried: ",
       paste(c(preferred, fallbacks), collapse = ", "))
}

rna_red  <- pick_reduction(combined_wnn, "rna_umap",  c("rnaumap"))
atac_red <- pick_reduction(combined_wnn, "atac_umap", c("atacumap"))
wnn_red  <- pick_reduction(combined_wnn, "wnn_umap",  c("wnnumap"))


################################################################################
# 8. FIGURE 1 PANEL A — Multimodal Integration UMAPs (RNA | ATAC | WNN)
################################################################################

cat("=== Panel A ===\n")

make_integration_umap <- function(obj, red, title) {
  DimPlot(obj, reduction = red, group.by = "dataset",
          pt.size = 0.1, raster = FALSE) +
    ggtitle(title) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text    = element_text(size = 12),
      axis.title   = element_text(size = 14),
      legend.text  = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold")
    )
}

panel_a <- make_integration_umap(combined_wnn, rna_red,  "RNA Integration") |
           make_integration_umap(combined_wnn, atac_red, "ATAC Integration") |
           make_integration_umap(combined_wnn, wnn_red,  "WNN Integration") +
  plot_layout(guides = "collect") & theme(legend.position = "right")

save_plot_pdf_eps(panel_a,
  file.path(OUT_DIR, "Figure1_PanelA_MultimodalIntegration.pdf"),
  file.path(OUT_DIR, "Figure1_PanelA_MultimodalIntegration.eps"),
  PDF_WIDTH_WIDE, PDF_HEIGHT_WIDE)
cat("Panel A saved.\n")


################################################################################
# 9. FIGURE 1 PANEL B — WNN Clusters UMAP (49 clusters)
################################################################################

cat("=== Panel B ===\n")

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
    plot.title      = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text       = element_text(size = 14),
    axis.title      = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )

save_plot_pdf_eps(panel_b,
  file.path(OUT_DIR, "Figure1_PanelB_WNNClusters.pdf"),
  file.path(OUT_DIR, "Figure1_PanelB_WNNClusters.eps"),
  PDF_WIDTH_SINGLE, PDF_HEIGHT_SINGLE)
cat("Panel B saved.\n")


################################################################################
# 10. FIGURE 1 PANEL C — Cell Type Annotations UMAP (13 cell types)
################################################################################

cat("=== Panel C ===\n")

panel_c <- DimPlot(
  combined_wnn,
  reduction = wnn_red,
  group.by  = "cluster_annotation",
  label     = FALSE,
  pt.size   = 0.1,
  raster    = FALSE
) +
  ggtitle("Hypothalamic Cell Types") +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text       = element_text(size = 14),
    axis.title      = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 12),
    legend.title    = element_text(size = 14, face = "bold"),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))

save_plot_pdf_eps(panel_c,
  file.path(OUT_DIR, "Figure1_PanelC_CellTypeAnnotations.pdf"),
  file.path(OUT_DIR, "Figure1_PanelC_CellTypeAnnotations.eps"),
  12, PDF_HEIGHT_SINGLE)
cat("Panel C saved.\n")


################################################################################
# 11. FIGURE 1 PANEL D — Cell Count Bar Plot (13 cell types)
################################################################################

cat("=== Panel D ===\n")

celltype_counts <- as.data.frame(table(combined_wnn$cluster_annotation)) %>%
  setNames(c("CellType", "Count")) %>%
  arrange(desc(Count)) %>%
  mutate(CellType = factor(CellType, levels = CellType))

panel_d <- ggplot(celltype_counts, aes(x = CellType, y = Count, fill = CellType)) +
  geom_col(color = "black", linewidth = 0.5) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Cell Count per Cell Type", x = "Cell Type", y = "Number of Cells") +
  theme_classic(base_size = 16) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y     = element_text(size = 12),
    axis.title      = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

save_plot_pdf_eps(panel_d,
  file.path(OUT_DIR, "Figure1_PanelD_CellTypeCounts.pdf"),
  file.path(OUT_DIR, "Figure1_PanelD_CellTypeCounts.eps"),
  PDF_WIDTH_SINGLE, 7)
cat("Panel D saved.\n")


################################################################################
# 12. FIGURE 1 PANEL E — Marker Gene DotPlot
#
# Clusters ordered by cell type (Y-axis), top 2 marker genes per cluster.
# cols = c("lightgrey", "blue"), dot.scale = 6, scale = TRUE
################################################################################

cat("=== Panel E ===\n")

all_genes_unique <- top_genes_per_wnn_cluster %>%
  filter(cluster %in% ordered_clusters) %>%
  arrange(match(cluster, ordered_clusters)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  pull(gene)

cat("Genes in DotPlot:", length(all_genes_unique), "\n")

# Temporarily set ordered cluster identities (grouped by cell type)
old_idents <- Idents(combined_wnn)

combined_wnn$wsnn_cluster_char <- as.character(combined_wnn$wsnn_res.0.6)
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
    plot.title      = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12.5),
    axis.text.y     = element_text(size = 10, face = "bold"),
    axis.title      = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 10),
    legend.title    = element_text(size = 11, face = "bold"),
    legend.position = "right",
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin     = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )

save_plot_pdf_eps(panel_e,
  file.path(OUT_DIR, "Figure1_PanelE_49WNNClusters_DotPlot_WIDE.pdf"),
  file.path(OUT_DIR, "Figure1_PanelE_49WNNClusters_DotPlot_WIDE.eps"),
  28, 10)

Idents(combined_wnn) <- old_idents    # restore
cat("Panel E saved.\n\n")


################################################################################
# 13. SAVE OUTPUT
################################################################################

saveRDS(combined_wnn, file.path(OUT_DIR, "04_combined_wnn_annotated.rds"))

cat("========================================\n")
cat("SCRIPT 04 COMPLETE\n")
cat("========================================\n")
cat("Nuclei annotated :", ncol(combined_wnn), "\n")
cat("Cell types        : 13 (from annotation spreadsheet)\n")
cat("WNN clusters      : 49\n\n")
cat("Cell type counts:\n")
print(sort(table(combined_wnn$cluster_annotation), decreasing = TRUE))

writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "session_info_04.txt"))
cat("\nNext: Run 05_Differential_Gene_Expression.R\n")
