################################################################################
# Script:  05_Differential_Gene_Expression.R
# Project: Heat_Call_Single_Nucleus_Multiome
# Author:  Prakrit Subba
# Date:    2026
#
# Paper:   Subba et al. 2026 — Multiome Profiling Reveals Astrocyte and
#          Neuroendocrine Targets of Prenatal Acoustic Programming in
#          Zebra Finch Embryos
#
# Description:
#   Runs cluster-by-cluster differential gene expression (DEG) analysis
#   between Heat Call and Control conditions across all 49 WNN clusters,
#   plus sex-stratified DEG for Cluster 16 (SIM1+ PVN glutamatergic).
#   Generates the DEG dot matrix figure and TTR volcano plots.
#
#   Analyses:
#     (1) All-cluster DEG: Heat Call vs Control (MAST, per WNN cluster)
#     (2) Filtering to padj < 0.05 → FilteredAllClustersHeatCallvsControlDEG.csv
#     (3) Sex-stratified DEG for Cluster 16 (Male and Female separately)
#     (4) DEG dot matrix (Figure 2A): Up/Down tiles per cluster, grouped
#         by cell type, excluding Cluster 26
#     (5) TTR volcano plots (Figure 2B): Male vs Female Cluster 16
#
# Inputs:
#   - results/04_Cell_Type_Annotation/04_combined_wnn_annotated.rds
#
# Outputs (exact file names used by downstream scripts):
#   - Per-cluster:   Cluster_<N>_DEG.csv              (one per cluster)
#   - Combined raw:  AllClustersHeatCallvsControlDEG.csv
#   - Filtered:      FilteredAllClustersHeatCallvsControlDEG.csv
#   - Sex-stratified Cluster 16:
#       Cluster16MaleDEG.csv
#       Cluster16FemaleDEG.csv
#   - Figures:
#       FigureDEGDotMatrixHeatCallvsControl.pdf/.eps
#       FigureDEGDotMatrixWithSeparators.pdf/.eps
#       FigureTTRVolcanoMaleFemaleHQ_sideBySide_sameY.pdf
#
# Column names in output CSVs (from FindMarkers + added columns):
#   p_val, avg_log2FC, pct.1, pct.2, p_val_adj, Cluster, Gene
#
# NOTE on naming:
#   - cluster_condition   : metadata column (underscore between words)
#   - AllClusters...DEG   : combined output files (NO underscores between words)
#   - Cluster_<N>_DEG.csv : per-cluster files (underscores around cluster number)
#   - Cluster16MaleDEG    : sex-stratified (NO underscores — exact name required
#                           by figure scripts)
#
# Previous Script: 04_Cell_Type_Annotation.R
# Next Script:     06_Differential_Chromatin_Accessibility.R
################################################################################


################################################################################
# 0. CONFIGURATION
################################################################################

WNN_RDS <- "results/04_Cell_Type_Annotation/04_combined_wnn_annotated.rds"

OUT_DIR <- "results/05_DEG"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# FindMarkers parameters (manuscript Methods)
DEG_TEST          <- "MAST"
DEG_MIN_PCT       <- 0.25
DEG_LFC_THRESHOLD <- 0.25
DEG_RECORRECT_UMI <- FALSE

# Significance threshold for filtered output
PADJ_CUTOFF <- 0.05

# Cluster excluded from dot matrix (Glutamatergic-GABAergic; ambiguous)
EXCLUDE_CLUSTER <- "26"

set.seed(1234)


################################################################################
# 1. LOAD LIBRARIES
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)


################################################################################
# 2. LOAD ANNOTATED WNN OBJECT
################################################################################

cat("=== Loading annotated WNN object ===\n")

combined_wnn <- readRDS(WNN_RDS)

cat("Cells           :", ncol(combined_wnn), "\n")
cat("WNN clusters    :", length(levels(combined_wnn$wsnn_res.0.6)), "\n")
cat("Conditions      :", paste(unique(combined_wnn$Condition), collapse = ", "), "\n")
cat("Sex categories  :", paste(unique(combined_wnn$Sex), collapse = ", "), "\n\n")

# Validate required metadata columns
stopifnot(
  "wsnn_res.0.6"      %in% colnames(combined_wnn@meta.data),
  "Condition"         %in% colnames(combined_wnn@meta.data),
  "Sex"               %in% colnames(combined_wnn@meta.data),
  "cluster_annotation" %in% colnames(combined_wnn@meta.data)
)


################################################################################
# 3. SET UP IDENTITIES FOR DEG
#
# Create cluster_condition metadata column by pasting WNN cluster and
# Condition with underscore separator. This produces identities like:
#   "0_Heat Call", "0_Control", "1_Heat Call", "1_Control", ...
#
# These exact strings are used as ident.1 and ident.2 in FindMarkers.
################################################################################

cat("=== Setting up cluster_condition identities ===\n")

DefaultAssay(combined_wnn) <- "SCT"

combined_wnn$cluster_condition <- paste(
  combined_wnn$wsnn_res.0.6,
  combined_wnn$Condition,
  sep = "_"
)

Idents(combined_wnn) <- "cluster_condition"

cat("Example identities (first 6):\n")
print(head(sort(unique(combined_wnn$cluster_condition))))
cat("\n")


################################################################################
# 4. ALL-CLUSTER DEG: HEAT CALL VS CONTROL
#
# Loops over all 49 WNN clusters. For each cluster, runs FindMarkers
# comparing Heat Call vs Control cells within that cluster only.
#
# Parameters:
#   ident.1         = paste0(cluster, "_Heat Call")
#   ident.2         = paste0(cluster, "_Control")
#   test.use        = "MAST"
#   recorrect_umi   = FALSE
#   min.pct         = 0.25
#   logfc.threshold = 0.25
#
# Adds $Cluster and $Gene columns before saving each per-cluster CSV.
# Per-cluster output: Cluster_<N>_DEG.csv
################################################################################

cat("=== Running all-cluster DEG (Heat Call vs Control) ===\n")
cat("(MAST test — computationally intensive, run on HPC)\n\n")

clusters  <- levels(combined_wnn$wsnn_res.0.6)
DEresults <- list()

for (cluster in clusters) {

  ident_hc   <- paste0(cluster, "_Heat Call")
  ident_ctrl <- paste0(cluster, "_Control")

  # Skip if either identity is absent (too few cells)
  if (!ident_hc   %in% Idents(combined_wnn) ||
      !ident_ctrl %in% Idents(combined_wnn)) {
    cat("  Cluster", cluster, ": skipped (identity not found)\n")
    next
  }

  cat("  Cluster", cluster, "...")

  DEresults[[cluster]] <- FindMarkers(
    object          = combined_wnn,
    ident.1         = ident_hc,
    ident.2         = ident_ctrl,
    test.use        = DEG_TEST,
    recorrect_umi   = DEG_RECORRECT_UMI,
    min.pct         = DEG_MIN_PCT,
    logfc.threshold = DEG_LFC_THRESHOLD
  )

  if (nrow(DEresults[[cluster]]) > 0) {
    DEresults[[cluster]]$Cluster <- cluster
    DEresults[[cluster]]$Gene    <- rownames(DEresults[[cluster]])

    write.csv(
      DEresults[[cluster]],
      file = file.path(OUT_DIR, paste0("Cluster_", cluster, "_DEG.csv"))
    )
    cat(" done —", nrow(DEresults[[cluster]]), "genes tested\n")
  } else {
    cat(" no results\n")
  }
}

cat("\nAll-cluster DEG complete.\n\n")


################################################################################
# 5. COMBINE AND SAVE ALL-CLUSTER DEG RESULTS
#
# Combines all per-cluster results into a single data frame.
# Output file: AllClustersHeatCallvsControlDEG.csv
# (No underscores between words — exact name used by downstream scripts)
################################################################################

cat("=== Combining all DEG results ===\n")

combined_de <- do.call(rbind, DEresults)

write.csv(
  combined_de,
  file      = file.path(OUT_DIR, "AllClustersHeatCallvsControlDEG.csv"),
  row.names = FALSE
)

cat("Total genes tested across all clusters:", nrow(combined_de), "\n")
cat("Saved: AllClustersHeatCallvsControlDEG.csv\n\n")


################################################################################
# 6. FILTER TO SIGNIFICANT DEGs
#
# Filters combined_de to p_val_adj < 0.05 to produce the filtered table
# used by the DEG dot matrix figure panel.
# Output: FilteredAllClustersHeatCallvsControlDEG.csv
################################################################################

cat("=== Filtering to significant DEGs (p_val_adj <", PADJ_CUTOFF, ") ===\n")

combined_de_filtered <- combined_de %>%
  filter(p_val_adj < PADJ_CUTOFF)

write.csv(
  combined_de_filtered,
  file      = file.path(OUT_DIR, "FilteredAllClustersHeatCallvsControlDEG.csv"),
  row.names = FALSE
)

cat("Significant DEGs:", nrow(combined_de_filtered), "\n")
cat("Clusters with DEGs:", length(unique(combined_de_filtered$Cluster)), "\n")
cat("Saved: FilteredAllClustersHeatCallvsControlDEG.csv\n\n")

# Summary by cluster
cat("Top clusters by DEG count:\n")
print(
  combined_de_filtered %>%
    count(Cluster, name = "n_DEGs") %>%
    arrange(desc(n_DEGs)) %>%
    head(10)
)
cat("\n")


################################################################################
# 7. SEX-STRATIFIED DEG — CLUSTER 16 (SIM1+ PVN GLUTAMATERGIC)
#
# Cluster 16 shows sex-dimorphic TTR regulation (manuscript Figure 2B):
#   - Males:   TTR upregulated in Heat Call
#   - Females: TTR downregulated in Heat Call
#
# Objects are subset by Sex, then FindMarkers is run for Cluster 16 only.
# Parameters are identical to the all-cluster analysis.
#
# Output files (exact names required by TTR volcano figure code):
#   Cluster16MaleDEG.csv      (no underscores — required by Figure 2B code)
#   Cluster16FemaleDEG.csv    (no underscores — required by Figure 2B code)
################################################################################

cat("=== Sex-stratified DEG — Cluster 16 (PVN) ===\n")

for (sex in c("Male", "Female")) {

  cat("  Sex:", sex, "...\n")

  # Subset to this sex
  sex_obj <- subset(combined_wnn, subset = Sex == sex)

  # Re-create cluster_condition on the subset
  sex_obj$cluster_condition <- paste(
    sex_obj$wsnn_res.0.6,
    sex_obj$Condition,
    sep = "_"
  )
  Idents(sex_obj) <- "cluster_condition"

  ident_hc   <- "16_Heat Call"
  ident_ctrl <- "16_Control"

  if (!ident_hc %in% Idents(sex_obj) || !ident_ctrl %in% Idents(sex_obj)) {
    cat("    Skipped — identity not found for", sex, "in Cluster 16\n")
    next
  }

  de_sex_cl16 <- FindMarkers(
    object          = sex_obj,
    ident.1         = ident_hc,
    ident.2         = ident_ctrl,
    test.use        = DEG_TEST,
    recorrect_umi   = DEG_RECORRECT_UMI,
    min.pct         = DEG_MIN_PCT,
    logfc.threshold = DEG_LFC_THRESHOLD
  )

  de_sex_cl16$Gene    <- rownames(de_sex_cl16)
  de_sex_cl16$Cluster <- "16"
  de_sex_cl16$Sex     <- sex

  # File names: Cluster16MaleDEG.csv / Cluster16FemaleDEG.csv
  # NOTE: No underscores — these exact names are required by the
  # TTR volcano figure code (Figure 2B)
  outfile <- file.path(OUT_DIR, paste0("Cluster16", sex, "DEG.csv"))

  write.csv(de_sex_cl16, file = outfile)

  cat("    Genes tested:", nrow(de_sex_cl16), "\n")

  # Check TTR
  if ("TTR" %in% de_sex_cl16$Gene) {
    ttr_row <- de_sex_cl16[de_sex_cl16$Gene == "TTR", ]
    cat("    TTR — avg_log2FC:", round(ttr_row$avg_log2FC, 3),
        "| p_val_adj:", signif(ttr_row$p_val_adj, 3), "\n")
  } else {
    cat("    TTR not in results for", sex, "\n")
  }
  cat("    Saved:", basename(outfile), "\n")
}

cat("\n")


################################################################################
# 8. DEG DOT MATRIX — FIGURE 2A
#
# Tile matrix showing Up (red) / Down (blue) DEGs per cluster × gene.
# Clusters ordered by cell type annotation on Y-axis.
# Cluster 26 excluded.
# Cell-type separator lines added in green.
#
# Input:  FilteredAllClustersHeatCallvsControlDEG.csv
# Output: FigureDEGDotMatrixHeatCallvsControl.pdf/.eps
#         FigureDEGDotMatrixWithSeparators.pdf/.eps
################################################################################

cat("=== Building DEG dot matrix (Figure 2A) ===\n")

# Validate filtered DEG file exists
deg_file <- file.path(OUT_DIR, "FilteredAllClustersHeatCallvsControlDEG.csv")
if (!file.exists(deg_file)) {
  stop("Cannot find: ", deg_file)
}

dedotmatrixplot <- read.csv(deg_file)

# Build cluster ordering (grouped by cell type, cluster 26 excluded)
cluster_to_celltype2 <- combined_wnn@meta.data %>%
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

ordered_clusters <- cluster_to_celltype2$`wsnn_res.0.6`
ordered_clusters <- ordered_clusters[ordered_clusters != EXCLUDE_CLUSTER]

# Prepare plot data
dedotmatrixplot <- dedotmatrixplot %>%
  mutate(
    Cluster = as.character(Cluster),
    Gene    = as.character(Gene)
  )

filtereddata <- dedotmatrixplot %>%
  filter(p_val_adj < PADJ_CUTOFF, Cluster != EXCLUDE_CLUSTER)

allgenes <- sort(unique(filtereddata$Gene))

cat("Significant DEGs (excl. cluster 26):", length(allgenes), "unique genes\n")
cat("Clusters with DEGs:", length(unique(filtereddata$Cluster)), "\n\n")

# Build complete grid (all cluster × gene combinations)
completegrid <- expand.grid(
  Cluster = ordered_clusters,
  Gene    = allgenes,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Cluster = as.character(Cluster),
    Gene    = as.character(Gene)
  )

plotdata <- left_join(
  completegrid,
  filtereddata %>% select(Cluster, Gene, avg_log2FC),
  by = c("Cluster", "Gene")
)

plotdata$Cluster <- factor(plotdata$Cluster, levels = ordered_clusters)

# Build dot matrix plot
dotmatrixdeg <- ggplot(
  plotdata,
  aes(x = factor(Gene, levels = allgenes), y = Cluster)
) +
  geom_tile(fill = "white", color = "black", linewidth = 0.3) +
  geom_tile(
    data = plotdata %>% filter(!is.na(avg_log2FC)),
    aes(fill = ifelse(avg_log2FC > 0, "Up", "Down")),
    color = "black", linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c("Up" = "red", "Down" = "blue"),
    name   = "Regulation"
  ) +
  labs(
    title = "Heat Call vs Control — Clusters Grouped by Cell Type",
    x     = "Gene",
    y     = "WNN Cluster (Grouped by Cell Type)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y     = element_text(size = 11, face = "bold"),
    axis.title      = element_text(size = 14, face = "bold"),
    panel.grid      = element_blank(),
    legend.text     = element_text(size = 12),
    legend.title    = element_text(size = 13, face = "bold"),
    legend.position = "right"
  )

# Save without separators
pdf(
  file.path(OUT_DIR, "FigureDEGDotMatrixHeatCallvsControl.pdf"),
  width = 16, height = 10, useDingbats = FALSE
)
print(dotmatrixdeg)
dev.off()
setEPS()
postscript(
  file.path(OUT_DIR, "FigureDEGDotMatrixHeatCallvsControl.eps"),
  width = 16, height = 10
)
print(dotmatrixdeg)
dev.off()

# Add cell-type separator lines (green, between cell type groups)
celltypeinfo <- cluster_to_celltype2 %>%
  filter(`wsnn_res.0.6` != EXCLUDE_CLUSTER) %>%
  mutate(clusterpos = match(`wsnn_res.0.6`, ordered_clusters)) %>%
  group_by(cluster_annotation) %>%
  summarise(lastpos = max(clusterpos), .groups = "drop")

dotmatrixdeg_separated <- dotmatrixdeg +
  geom_hline(
    yintercept = celltypeinfo$lastpos + 0.5,
    color      = "darkgreen",
    linewidth  = 1
  )

pdf(
  file.path(OUT_DIR, "FigureDEGDotMatrixWithSeparators.pdf"),
  width = 16, height = 10, useDingbats = FALSE
)
print(dotmatrixdeg_separated)
dev.off()
setEPS()
postscript(
  file.path(OUT_DIR, "FigureDEGDotMatrixWithSeparators.eps"),
  width = 16, height = 10
)
print(dotmatrixdeg_separated)
dev.off()

cat("DEG dot matrix figures saved.\n\n")


################################################################################
# 9. TTR VOLCANO PLOTS — FIGURE 2B
#
# Side-by-side volcano plots for Cluster 16 (SIM1+ PVN glutamatergic):
#   Left  : Male   Heat Call vs Control  (TTR colored red3)
#   Right : Female Heat Call vs Control  (TTR colored dodgerblue3)
#
# Shared Y-axis across both panels (same y-limits).
# LFC cutoff lines at ±1; significance line at padj = 0.05.
#
# Input files (must exist):
#   Cluster16MaleDEG.csv
#   Cluster16FemaleDEG.csv
################################################################################

cat("=== Building TTR volcano plots (Figure 2B) ===\n")

male_file   <- file.path(OUT_DIR, "Cluster16MaleDEG.csv")
female_file <- file.path(OUT_DIR, "Cluster16FemaleDEG.csv")

if (!file.exists(male_file))   stop("Missing: ", male_file)
if (!file.exists(female_file)) stop("Missing: ", female_file)

male_df   <- read.csv(male_file)
female_df <- read.csv(female_file)

# Helper: build TTR-highlighted volcano plot
make_ttr_volcano <- function(df, title,
                              ttr_color = "red3",
                              lfc_cut   = 1,
                              p_cut     = 0.05,
                              xlim      = c(-2, 1.2),
                              ylim      = NULL) {

  if (!all(c("Gene", "avg_log2FC") %in% colnames(df))) {
    stop("DEG table must contain 'Gene' and 'avg_log2FC'.")
  }

  # Use p_val_adj if available, else p_val
  pcol <- if ("p_val_adj" %in% colnames(df)) "p_val_adj" else
          if ("p_val"     %in% colnames(df)) "p_val"     else
          stop("Need p_val or p_val_adj column.")

  df          <- df[!is.na(df[[pcol]]) & !is.na(df$avg_log2FC), ]
  df$neglog10p <- -log10(df[[pcol]])

  ttr_row <- df[df$Gene == "TTR", , drop = FALSE]

  if (is.null(ylim)) {
    ymax <- max(df$neglog10p, na.rm = TRUE)
    ylim <- c(0, ymax * 1.05)
  }

  ggplot(df, aes(x = avg_log2FC, y = neglog10p)) +
    geom_point(color = "grey65", size = 1.8, alpha = 0.75) +
    geom_hline(
      yintercept = -log10(p_cut),
      linetype = "dotted", linewidth = 0.8, color = "grey20"
    ) +
    geom_vline(
      xintercept = c(-lfc_cut, lfc_cut),
      linetype = "dotted", linewidth = 0.8, color = "grey20"
    ) +
    geom_point(data = ttr_row, color = ttr_color, size = 3.2) +
    geom_text(
      data     = ttr_row,
      aes(label = Gene),
      color    = ttr_color,
      fontface = "bold",
      size     = 5,
      nudge_x  = ifelse(nrow(ttr_row) == 0, 0,
                        ifelse(ttr_row$avg_log2FC > 0, 0.18, -0.18)),
      nudge_y  = 0.6
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim, clip = "off") +
    labs(
      title = title,
      x     = "avg_log2FC",
      y     = "-log10(p-value)"
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title  = element_text(face = "bold", size = 18, hjust = 0.5),
      axis.title  = element_text(face = "bold", size = 16),
      axis.text   = element_text(size = 14, color = "black"),
      axis.line   = element_line(linewidth = 1.1, color = "black"),
      axis.ticks  = element_line(linewidth = 1,   color = "black"),
      plot.margin = margin(8, 14, 8, 8)
    )
}

# Helper to compute -log10(p) max
get_neglog10p <- function(df) {
  pcol <- if ("p_val_adj" %in% colnames(df)) "p_val_adj" else "p_val"
  df   <- df[!is.na(df[[pcol]]) & !is.na(df$avg_log2FC), ]
  -log10(df[[pcol]])
}

# Shared Y-axis across both panels
ymax_shared <- max(
  c(get_neglog10p(male_df), get_neglog10p(female_df)),
  na.rm = TRUE
)
ylim_shared <- c(0, ymax_shared * 1.05)

p_male   <- make_ttr_volcano(
  male_df,
  title     = "Male Glutamatergic Cluster 16",
  ttr_color = "red3",
  lfc_cut   = 1,
  p_cut     = 0.05,
  xlim      = c(-2, 1.2),
  ylim      = ylim_shared
)

p_female <- make_ttr_volcano(
  female_df,
  title     = "Female Glutamatergic Cluster 16",
  ttr_color = "dodgerblue3",
  lfc_cut   = 1,
  p_cut     = 0.05,
  xlim      = c(-2, 1.2),
  ylim      = ylim_shared
)

# Save side-by-side with shared Y-axis
pdf(
  file.path(OUT_DIR, "FigureTTRVolcanoMaleFemaleHQ_sideBySide_sameY.pdf"),
  width = 14, height = 6, useDingbats = FALSE
)
print(p_male + p_female)
dev.off()
setEPS()
postscript(
  file.path(OUT_DIR, "FigureTTRVolcanoMaleFemaleHQ_sideBySide_sameY.eps"),
  width = 14, height = 6
)
print(p_male + p_female)
dev.off()

cat("TTR volcano plots saved.\n\n")


################################################################################
# 10. SUMMARY
################################################################################

cat("========================================\n")
cat("SCRIPT 05 COMPLETE — DEG SUMMARY\n")
cat("========================================\n")
cat("Total clusters tested  :", length(DEresults), "\n")
cat("Total genes tested     :", nrow(combined_de), "\n")
cat("Significant (padj<0.05):", nrow(combined_de_filtered), "\n")
cat("  Upregulated (HC>CTRL):",
    sum(combined_de_filtered$avg_log2FC > 0, na.rm = TRUE), "\n")
cat("  Downregulated        :",
    sum(combined_de_filtered$avg_log2FC < 0, na.rm = TRUE), "\n")
cat("\nOutput files:\n")
cat("  AllClustersHeatCallvsControlDEG.csv\n")
cat("  FilteredAllClustersHeatCallvsControlDEG.csv\n")
cat("  Cluster16MaleDEG.csv\n")
cat("  Cluster16FemaleDEG.csv\n")
cat("  Cluster_<N>_DEG.csv  (one per cluster, in OUT_DIR)\n")
cat("\nFigures:\n")
cat("  FigureDEGDotMatrixHeatCallvsControl.pdf/.eps\n")
cat("  FigureDEGDotMatrixWithSeparators.pdf/.eps\n")
cat("  FigureTTRVolcanoMaleFemaleHQ_sideBySide_sameY.pdf/.eps\n\n")


################################################################################
# 11. SAVE OBJECT (unchanged — passed through for downstream scripts)
################################################################################

saveRDS(
  combined_wnn,
  file = file.path(OUT_DIR, "05_combined_wnn_DEG.rds")
)

cat("Object saved: 05_combined_wnn_DEG.rds\n")
cat("Next: Run 06_Differential_Chromatin_Accessibility.R\n\n")

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "session_info_05.txt")
)
