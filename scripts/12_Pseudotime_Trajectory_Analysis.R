# ==============================================================================
# 12_Pseudotime_Trajectory_Analysis.R
# ------------------------------------------------------------------------------
# Slingshot pseudotime trajectory inference + tradeSeq GAM fitting +
# per-cluster developmental shift analysis + smooth expression plots (Fig 5)
# ==============================================================================
# Depends on:  29May2025DAconditioncellcluster.RData (combined_wnn)
# Outputs:     PseudotimeTrajectory_AllCells.csv
#              tradeSeq_AssociationTest_Control.csv
#              tradeSeq_AssociationTest_HeatCall.csv
#              PerClusterPseudotimeShift.csv
#              Fig5A_TrajectoryUMAP.pdf
#              Fig5B_SmoothExpressionPlots.pdf
# ==============================================================================


# ── 0. LIBRARIES ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(slingshot)
  library(tradeSeq)
  library(SingleCellExperiment)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(scales)
  library(viridis)
})
set.seed(1234)
theme_set(theme_cowplot(font_size = 12))


# ── 1. CONFIGURATION ──────────────────────────────────────────────────────────
RDATA_PATH    <- "29May2025DAconditioncellcluster.RData"
OUTDIR        <- "pseudotime_output"
FIGURES       <- file.path(OUTDIR, "figures")
TABLES        <- file.path(OUTDIR, "tables")

# Slingshot / tradeSeq parameters
START_CLUSTER <- "0"      # progenitor cluster (trajectory origin)
N_KNOTS       <- 6        # tradeSeq GAM knots per condition
ASSOC_PVAL    <- 1e-10    # threshold for pseudotime-dynamic genes
N_CORES       <- 4        # parallel cores for tradeSeq

# Genes for smooth expression plots (Figure 5B)
SMOOTH_GENES  <- c("ASCL1", "DLL1", "NOTCH1", "MAML3", "NFIA", "FABP7")

# Conditions
COND1 <- "Control"
COND2 <- "Heat Call"

dir.create(OUTDIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES,  recursive = TRUE, showWarnings = FALSE)

cat("Script: 12_Pseudotime_Trajectory_Analysis.R\n")
cat("OUTDIR:", OUTDIR, "\n")


# ── 2. LOAD DATA ──────────────────────────────────────────────────────────────
cat("\n=== STEP 1: LOADING DATA ===\n")

load(RDATA_PATH)
stopifnot(exists("combined_wnn"))
stopifnot("Condition" %in% colnames(combined_wnn@meta.data))
stopifnot("cluster_annotation" %in% colnames(combined_wnn@meta.data))

cat("Total cells:", ncol(combined_wnn), "\n")
cat("WNN clusters:", length(unique(combined_wnn$seurat_clusters)), "\n")
print(table(combined_wnn$Condition))

# Extract WNN UMAP embedding (used as input to Slingshot)
wnn_umap <- Embeddings(combined_wnn, reduction = "wnn.umap")
cat("WNN UMAP dims:", paste(dim(wnn_umap), collapse = " x "), "\n")


# ── 3. SLINGSHOT TRAJECTORY — BY CONDITION ────────────────────────────────────
cat("\n=== STEP 2: SLINGSHOT TRAJECTORY INFERENCE ===\n")

run_slingshot <- function(seurat_obj, condition_label, start_clust) {
  cells     <- colnames(seurat_obj)[seurat_obj$Condition == condition_label]
  sub_obj   <- subset(seurat_obj, cells = cells)
  umap_sub  <- Embeddings(sub_obj, reduction = "wnn.umap")
  clust_sub <- as.character(sub_obj$seurat_clusters)

  sce <- SingleCellExperiment(
    assays      = list(counts = GetAssayData(sub_obj, assay = "RNA",
                                             slot = "counts")),
    reducedDims = list(UMAP = umap_sub)
  )
  colLabels(sce) <- clust_sub

  cat(sprintf("  Running Slingshot for %s (%d cells)...\n",
              condition_label, ncol(sce)))
  sds <- slingshot(sce,
                   clusterLabels = colLabels(sce),
                   reducedDim    = "UMAP",
                   start.clus    = start_clust)

  pt <- slingPseudotime(sds, na = FALSE)
  # Use principal curve pseudotime (first lineage if multiple)
  pt_vec <- rowMeans(pt, na.rm = TRUE)

  data.frame(
    cell      = colnames(sce),
    pseudotime = pt_vec,
    condition = condition_label,
    cluster   = clust_sub,
    stringsAsFactors = FALSE
  )
}

pt_ctrl <- run_slingshot(combined_wnn, COND1, START_CLUSTER)
pt_heat <- run_slingshot(combined_wnn, COND2, START_CLUSTER)

pseudotime_all <- bind_rows(pt_ctrl, pt_heat)
write.csv(pseudotime_all,
          file.path(TABLES, "PseudotimeTrajectory_AllCells.csv"),
          row.names = FALSE)

cat(sprintf("Pseudotime range (%s): %.1f – %.1f\n",
            COND1, min(pt_ctrl$pseudotime), max(pt_ctrl$pseudotime)))
cat(sprintf("Pseudotime range (%s): %.1f – %.1f\n",
            COND2, min(pt_heat$pseudotime), max(pt_heat$pseudotime)))


# ── 4. PER-CLUSTER PSEUDOTIME SHIFT ANALYSIS ──────────────────────────────────
cat("\n=== STEP 3: PER-CLUSTER PSEUDOTIME SHIFT (Wilcoxon) ===\n")

# Normalize pseudotime within each condition to 0–1 scale for comparability
pseudotime_all <- pseudotime_all %>%
  group_by(condition) %>%
  mutate(pseudotime_norm = (pseudotime - min(pseudotime)) /
                           (max(pseudotime) - min(pseudotime))) %>%
  ungroup()

clusters_to_test <- sort(unique(pseudotime_all$cluster))

cluster_results <- map_dfr(clusters_to_test, function(cl) {
  ctrl_pt <- pseudotime_all %>%
    filter(cluster == cl, condition == COND1) %>%
    pull(pseudotime_norm)
  heat_pt <- pseudotime_all %>%
    filter(cluster == cl, condition == COND2) %>%
    pull(pseudotime_norm)

  if (length(ctrl_pt) < 5 || length(heat_pt) < 5) {
    return(tibble(cluster = cl, n_ctrl = length(ctrl_pt),
                  n_heat = length(heat_pt),
                  median_ctrl = NA_real_, median_heat = NA_real_,
                  delta_pseudotime = NA_real_, pval = NA_real_,
                  pval_adj = NA_real_))
  }

  wt <- wilcox.test(heat_pt, ctrl_pt, alternative = "two.sided")
  tibble(
    cluster          = cl,
    n_ctrl           = length(ctrl_pt),
    n_heat           = length(heat_pt),
    median_ctrl      = median(ctrl_pt),
    median_heat      = median(heat_pt),
    delta_pseudotime = median(heat_pt) - median(ctrl_pt),
    pval             = wt$p.value,
    pval_adj         = NA_real_   # filled below
  )
})

cluster_results$pval_adj <- p.adjust(cluster_results$pval, method = "BH")

cluster_results <- cluster_results %>%
  left_join(
    combined_wnn@meta.data %>%
      select(seurat_clusters, cluster_annotation) %>%
      distinct() %>%
      rename(cluster = seurat_clusters),
    by = "cluster"
  ) %>%
  arrange(pval_adj)

n_sig_clusters <- sum(cluster_results$pval_adj < 0.05, na.rm = TRUE)
cat(sprintf("Clusters with significant pseudotime shift (FDR < 0.05): %d\n",
            n_sig_clusters))
cat("Top 10 shifted clusters:\n")
print(head(cluster_results %>%
             select(cluster, cluster_annotation, n_ctrl, n_heat,
                    delta_pseudotime, pval, pval_adj), 10))

write.csv(cluster_results,
          file.path(TABLES, "PerClusterPseudotimeShift.csv"),
          row.names = FALSE)


# ── 5. FIGURE 5A: TRAJECTORY UMAP ────────────────────────────────────────────
cat("\n=== STEP 4: FIGURE 5A — TRAJECTORY UMAP ===\n")

# Add pseudotime and annotation to combined metadata for plotting
meta_pt <- combined_wnn@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  left_join(pseudotime_all %>% select(cell, pseudotime_norm), by = "cell") %>%
  tibble::column_to_rownames("cell")

# UMAP coordinates
umap_df <- as.data.frame(wnn_umap) %>%
  tibble::rownames_to_column("cell") %>%
  rename(UMAP1 = 2, UMAP2 = 3) %>%
  left_join(meta_pt %>% tibble::rownames_to_column("cell"), by = "cell")

p5a_ctrl <- ggplot(umap_df %>% filter(Condition == COND1),
                   aes(x = UMAP1, y = UMAP2, color = pseudotime_norm)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_viridis(option = "C", name = "Pseudotime\n(normalized)",
                      na.value = "grey80") +
  labs(title = COND1) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p5a_heat <- ggplot(umap_df %>% filter(Condition == COND2),
                   aes(x = UMAP1, y = UMAP2, color = pseudotime_norm)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_viridis(option = "C", name = "Pseudotime\n(normalized)",
                      na.value = "grey80") +
  labs(title = COND2) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p5a <- p5a_ctrl | p5a_heat
ggsave(file.path(FIGURES, "Fig5A_TrajectoryUMAP.pdf"),
       p5a, width = 12, height = 5)
ggsave(file.path(FIGURES, "Fig5A_TrajectoryUMAP.png"),
       p5a, width = 12, height = 5, dpi = 300)
cat("Saved Fig5A_TrajectoryUMAP\n")


# ── 6. TRADSEQ — GAM FITTING BY CONDITION ─────────────────────────────────────
cat("\n=== STEP 5: tradeSeq GAM FITTING ===\n")

fit_tradseq <- function(seurat_obj, pt_df, condition_label, n_knots, n_cores) {
  cells   <- pt_df$cell[pt_df$condition == condition_label]
  sub_obj <- subset(seurat_obj, cells = cells)
  counts  <- GetAssayData(sub_obj, assay = "RNA", slot = "counts")

  pt_vec  <- pt_df %>%
    filter(condition == condition_label) %>%
    arrange(match(cell, colnames(sub_obj))) %>%
    pull(pseudotime_norm)

  cat(sprintf("  Fitting GAM for %s (%d cells, %d knots)...\n",
              condition_label, ncol(sub_obj), n_knots))

  # Gene filter: expressed in ≥ 1% of cells
  min_cells <- max(10, round(0.01 * ncol(sub_obj)))
  keep_genes <- rowSums(counts > 0) >= min_cells
  cat(sprintf("  Genes passing expression filter: %d\n", sum(keep_genes)))

  sce <- fitGAM(counts   = counts[keep_genes, ],
                pseudotime = matrix(pt_vec, ncol = 1),
                cellWeights = matrix(1, nrow = ncol(sub_obj), ncol = 1),
                nknots   = n_knots,
                parallel = n_cores > 1,
                BPPARAM  = BiocParallel::MulticoreParam(n_cores,
                                                        progressbar = TRUE))
  sce
}

sce_ctrl <- fit_tradseq(combined_wnn, pseudotime_all, COND1, N_KNOTS, N_CORES)
sce_heat <- fit_tradseq(combined_wnn, pseudotime_all, COND2, N_KNOTS, N_CORES)

saveRDS(sce_ctrl, file.path(OUTDIR, "tradeSeq_SCE_Control.rds"))
saveRDS(sce_heat, file.path(OUTDIR, "tradeSeq_SCE_HeatCall.rds"))
cat("tradeSeq SCE objects saved.\n")


# ── 7. ASSOCIATION TEST ───────────────────────────────────────────────────────
cat("\n=== STEP 6: ASSOCIATION TEST (pseudotime-dynamic genes) ===\n")

assoc_ctrl <- associationTest(sce_ctrl)
assoc_heat <- associationTest(sce_heat)

assoc_ctrl$gene      <- rownames(assoc_ctrl)
assoc_heat$gene      <- rownames(assoc_heat)
assoc_ctrl$pval_adj  <- p.adjust(assoc_ctrl$pvalue, method = "BH")
assoc_heat$pval_adj  <- p.adjust(assoc_heat$pvalue, method = "BH")

write.csv(assoc_ctrl %>% arrange(pval_adj),
          file.path(TABLES, "tradeSeq_AssociationTest_Control.csv"),
          row.names = FALSE)
write.csv(assoc_heat %>% arrange(pval_adj),
          file.path(TABLES, "tradeSeq_AssociationTest_HeatCall.csv"),
          row.names = FALSE)

sig_ctrl <- sum(assoc_ctrl$pval_adj < ASSOC_PVAL, na.rm = TRUE)
sig_heat <- sum(assoc_heat$pval_adj < ASSOC_PVAL, na.rm = TRUE)
shared   <- length(intersect(
  assoc_ctrl$gene[assoc_ctrl$pval_adj < ASSOC_PVAL],
  assoc_heat$gene[assoc_heat$pval_adj < ASSOC_PVAL]
))

cat(sprintf("Pseudotime-dynamic genes (%s):     %d\n", COND1, sig_ctrl))
cat(sprintf("Pseudotime-dynamic genes (%s): %d\n", COND2, sig_heat))
cat(sprintf("Shared between both conditions:     %d\n", shared))
cat(sprintf("Control-specific:                   %d\n", sig_ctrl - shared))
cat(sprintf("Heat call-specific:                 %d\n", sig_heat - shared))

# Save dynamic gene counts summary
dyn_summary <- tibble(
  Metric = c(
    paste("Pseudotime-dynamic genes:", COND1),
    paste("Pseudotime-dynamic genes:", COND2),
    "Shared between conditions",
    paste("Condition-specific:", COND1),
    paste("Condition-specific:", COND2)
  ),
  Count = c(sig_ctrl, sig_heat, shared,
            sig_ctrl - shared, sig_heat - shared)
)
write.csv(dyn_summary,
          file.path(TABLES, "DynamicGeneCounts_Summary.csv"),
          row.names = FALSE)


# ── 8. FIGURE 5B: SMOOTH EXPRESSION PLOTS ────────────────────────────────────
cat("\n=== STEP 7: FIGURE 5B — SMOOTH EXPRESSION PLOTS ===\n")

# Build smooth expression curves using predictSmooth for each gene + condition
get_smooth_df <- function(sce, gene, condition_label, n_points = 200) {
  if (!gene %in% rownames(sce)) return(NULL)
  sm <- tryCatch(
    predictSmooth(sce, gene = gene, nPoints = n_points),
    error = function(e) NULL
  )
  if (is.null(sm)) return(NULL)
  data.frame(
    pseudotime = seq(0, 1, length.out = n_points),
    expression = as.numeric(sm),
    gene       = gene,
    condition  = condition_label,
    stringsAsFactors = FALSE
  )
}

smooth_list <- list()
for (gene in SMOOTH_GENES) {
  df_ctrl <- get_smooth_df(sce_ctrl, gene, COND1)
  df_heat <- get_smooth_df(sce_heat, gene, COND2)
  if (!is.null(df_ctrl) && !is.null(df_heat))
    smooth_list[[gene]] <- bind_rows(df_ctrl, df_heat)
}

smooth_all <- bind_rows(smooth_list)

if (nrow(smooth_all) > 0) {
  p5b <- ggplot(smooth_all,
                aes(x = pseudotime, y = expression,
                    color = condition, linetype = condition)) +
    geom_line(linewidth = 1.0) +
    facet_wrap(~ gene, scales = "free_y", ncol = 3) +
    scale_color_manual(
      values = setNames(c("#0072B2", "#D55E00"), c(COND1, COND2))
    ) +
    scale_linetype_manual(
      values = setNames(c("solid", "dashed"), c(COND1, COND2))
    ) +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(
      x        = "Normalized pseudotime",
      y        = "Smoothed expression",
      color    = "Condition",
      linetype = "Condition",
      title    = "Key differentiation regulators along astrocyte trajectory"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold"),
      legend.position  = "bottom"
    )

  ggsave(file.path(FIGURES, "Fig5B_SmoothExpressionPlots.pdf"),
         p5b, width = 11, height = 7)
  ggsave(file.path(FIGURES, "Fig5B_SmoothExpressionPlots.png"),
         p5b, width = 11, height = 7, dpi = 300)
  cat("Saved Fig5B_SmoothExpressionPlots\n")
} else {
  cat("WARNING: No smooth expression data generated — check gene names.\n")
}


# ── 9. ASTROCYTE-SPECIFIC PSEUDOTIME PLOTS ───────────────────────────────────
cat("\n=== STEP 8: ASTROCYTE SUBPOPULATION PSEUDOTIME ===\n")

astro_pt <- pseudotime_all %>%
  filter(cell %in% colnames(combined_wnn)[
    combined_wnn$cluster_annotation == "Astrocyte"])

p_astro_violin <- ggplot(astro_pt,
                          aes(x = condition, y = pseudotime_norm,
                              fill = condition)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.alpha = 0, fill = "white") +
  scale_fill_manual(values = c("#0072B2", "#D55E00"),
                    aesthetics = c("fill", "color")) +
  labs(x = NULL, y = "Normalized pseudotime",
       title = "Astrocyte pseudotime distribution by condition") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(FIGURES, "AstroPseudotimeViolin.pdf"),
       p_astro_violin, width = 5, height = 5)

# Most accelerated astrocyte clusters
top_accel <- cluster_results %>%
  filter(cluster_annotation == "Astrocyte",
         pval_adj < 0.05,
         delta_pseudotime > 0) %>%
  arrange(pval_adj) %>%
  head(5)
cat("Top accelerated astrocyte clusters (HC > Control):\n")
print(top_accel %>%
        select(cluster, cluster_annotation, delta_pseudotime, pval, pval_adj))


# ── 10. FINAL SUMMARY ────────────────────────────────────────────────────────
cat("\n========== 12_Pseudotime_Trajectory_Analysis.R COMPLETE ==========\n")
cat(sprintf("Pseudotime range (%s):         %.1f – %.1f\n",
            COND1, min(pt_ctrl$pseudotime), max(pt_ctrl$pseudotime)))
cat(sprintf("Pseudotime range (%s):     %.1f – %.1f\n",
            COND2, min(pt_heat$pseudotime), max(pt_heat$pseudotime)))
cat(sprintf("Significant shifted clusters:  %d\n", n_sig_clusters))
cat(sprintf("Dynamic genes (%s):       %d\n", COND1, sig_ctrl))
cat(sprintf("Dynamic genes (%s):   %d\n", COND2, sig_heat))
cat("\nKey outputs:\n")
cat("  - figures/Fig5A_TrajectoryUMAP.pdf\n")
cat("  - figures/Fig5B_SmoothExpressionPlots.pdf\n")
cat("  - figures/AstroPseudotimeViolin.pdf\n")
cat("  - tables/PseudotimeTrajectory_AllCells.csv\n")
cat("  - tables/PerClusterPseudotimeShift.csv          (Supplemental Table 11)\n")
cat("  - tables/tradeSeq_AssociationTest_Control.csv   (Supplementary File 1)\n")
cat("  - tables/tradeSeq_AssociationTest_HeatCall.csv  (Supplementary File 1)\n")
cat("  - tables/DynamicGeneCounts_Summary.csv\n")
cat("\nPipeline complete. All 12 scripts done.\n")
