################################################################################
# Script:  02_ATAC_QC_and_Processing.R
# Project: Heat_Call_Single_Nucleus_Multiome
# Author:  Prakrit Subba
# Date:    2026
#
# Paper:   Subba et al. 2026 — Multiome Profiling Reveals Astrocyte and
#          Neuroendocrine Targets of Prenatal Acoustic Programming in
#          Zebra Finch Embryos
#
# Description:
#   Processes the snATAC-seq modality from 10X Chromium Single Cell Multiome
#   output for 8 embryonic zebra finch hypothalamus samples. Performs:
#     (1) Cross-sample peak set unification via GRanges reduce()
#     (2) Fragment object and feature matrix creation
#     (3) Chromatin assay construction and sample merging
#     (4) QC metric computation (nucleosome signal, TSS enrichment,
#         fragment counts, reads-in-peaks fraction)
#     (5) QC filtering
#     (6) TF-IDF normalization, top feature selection, and LSI dimensionality
#         reduction (SVD)
#
# Inputs:
#   - <BASE_DIR>/<SAMPLE_ID>_combined_analysis/outs/atac_peaks.bed
#   - <BASE_DIR>/<SAMPLE_ID>_combined_analysis/outs/per_barcode_metrics.csv
#   - <BASE_DIR>/<SAMPLE_ID>_combined_analysis/outs/atac_fragments.tsv.gz
#     (and companion .tsv.gz.tbi index file)
#
# Outputs:
#   - combined_filt : QC-filtered, TF-IDF normalized, LSI-reduced Signac object
#   - QC plots (pre- and post-filter)
#   - 02_combined_atac_filt.rds : Saved Signac/Seurat object for Script 03
#
# Manuscript Reference:
#   Methods — snATAC-seq Processing and QC
#   Figure 1B (ATAC UMAP, carried forward to Script 03)
#   Results — 70,720 of 76,730 nuclei retained post-QC
#
# Requires:  01_RNA_QC_and_Processing.R to have been run (parallel, not
#            sequential — both feed into Script 03)
# Next Script: 03_WNN_Integration_and_Clustering.R
################################################################################


################################################################################
# 0. CONFIGURATION — Update these paths before running
################################################################################

# Root directory containing per-sample CellRanger ARC output folders
# Expected structure: <BASE_DIR>/S1_combined_analysis/outs/
BASE_DIR <- "/path/to/cellranger_arc_outputs"

# Output directory for QC plots and saved objects
OUT_DIR <- "results/02_ATAC_QC"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Sample IDs matching CellRanger ARC output folder prefixes
SAMPLES <- paste0("S", 1:8)

# QC thresholds (as reported in Methods)
MIN_ATAC_COUNTS      <- 100
MAX_ATAC_COUNTS      <- 100000
MIN_PCT_READS_PEAKS  <- 15
MAX_NUCLEOSOME_SIG   <- 4
MIN_TSS_ENRICHMENT   <- 2

# Fragment pre-filter: minimum ATAC fragments per barcode before object creation
MIN_FRAGMENTS        <- 500

# Peak width filter: remove peaks that are too small (likely artifacts)
# or too large (likely structural variants)
MIN_PEAK_WIDTH       <- 20
MAX_PEAK_WIDTH       <- 10000

# TF-IDF feature selection cutoff
TFIDF_MIN_CUTOFF     <- 20


################################################################################
# 1. LOAD LIBRARIES
################################################################################

library(Signac)
library(Seurat)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(patchwork)


################################################################################
# 2. EXPERIMENTAL METADATA
#
# Assigns Condition (Heat Call / Control) and Sex to each sample.
# Matches the design in 01_RNA_QC_and_Processing.R.
################################################################################

metadata <- data.frame(
  Sample    = SAMPLES,
  Condition = c("Heat Call", "Heat Call", "Control",   "Control",
                "Heat Call", "Heat Call", "Control",   "Control"),
  Sex       = c("Male",  "Female", "Male",  "Female",
                "Male",  "Female", "Male",  "Female"),
  stringsAsFactors = FALSE
)

cat("Experimental design:\n")
print(metadata)
cat("\n")


################################################################################
# 3. CREATE UNIFIED PEAK SET
#
# ATAC peak calls from CellRanger ARC are per-sample. To enable cross-sample
# comparison, peaks are read from each sample's BED file, converted to
# GRanges, and merged into a single non-overlapping consensus peak set using
# GenomicRanges::reduce().
#
# Peaks outside 20–10,000 bp are removed:
#   - < 20 bp   : likely sequencing artifacts
#   - > 10,000 bp : likely structural variant signals, not open chromatin
################################################################################

cat("=== Creating unified consensus peak set ===\n")

gr.list <- lapply(SAMPLES, function(s) {
  peaks_path <- file.path(
    BASE_DIR, paste0(s, "_combined_analysis"), "outs", "atac_peaks.bed"
  )
  peaks <- read.table(
    peaks_path,
    col.names = c("chr", "start", "end")
  )
  makeGRangesFromDataFrame(peaks)
})

# Merge overlapping peaks across all samples into consensus peak set
combined.peaks <- reduce(do.call(c, gr.list))

# Filter by peak width
peakwidths    <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < MAX_PEAK_WIDTH &
                                 peakwidths > MIN_PEAK_WIDTH]

cat("Total peaks in consensus set:", length(combined.peaks), "\n")
cat("Peak width range:",
    min(width(combined.peaks)), "–", max(width(combined.peaks)), "bp\n\n")


################################################################################
# 4. LOAD AND FILTER PER-BARCODE METRICS
#
# CellRanger ARC outputs per_barcode_metrics.csv with ATAC QC stats per cell
# barcode. Barcodes with fewer than 500 ATAC fragments are removed before
# object creation to pre-filter low-quality cells.
################################################################################

cat("=== Loading per-barcode metrics ===\n")

md.list <- lapply(SAMPLES, function(s) {
  metrics_path <- file.path(
    BASE_DIR, paste0(s, "_combined_analysis"), "outs", "per_barcode_metrics.csv"
  )
  md <- read.table(
    metrics_path,
    stringsAsFactors = FALSE,
    sep              = ",",
    header           = TRUE,
    row.names        = 1
  )
  # Pre-filter: keep only barcodes with sufficient ATAC fragments
  md[md$atac_fragments > MIN_FRAGMENTS, ]
})

# Add sample-level metadata to each barcode metrics table
for (i in seq_along(SAMPLES)) {
  md.list[[i]]$dataset   <- SAMPLES[i]
  md.list[[i]]$Condition <- metadata$Condition[metadata$Sample == SAMPLES[i]]
  md.list[[i]]$Sex       <- metadata$Sex[metadata$Sample == SAMPLES[i]]
}

# Report barcode counts per sample after pre-filter
cat("Barcodes retained per sample after fragment pre-filter:\n")
for (i in seq_along(SAMPLES)) {
  cat("  ", SAMPLES[i], ":", nrow(md.list[[i]]), "barcodes\n")
}
cat("\n")


################################################################################
# 5. CREATE FRAGMENT OBJECTS
#
# CreateFragmentObject() links each sample's fragment file (.tsv.gz) to its
# filtered barcode list. The companion index (.tsv.gz.tbi) must exist in the
# same directory.
################################################################################

cat("=== Creating fragment objects ===\n")

frags.list <- mapply(
  function(s, md) {
    frag_path <- file.path(
      BASE_DIR, paste0(s, "_combined_analysis"), "outs", "atac_fragments.tsv.gz"
    )
    CreateFragmentObject(
      path  = frag_path,
      cells = rownames(md)
    )
  },
  SAMPLES, md.list,
  SIMPLIFY = FALSE
)

cat("Fragment objects created for all", length(SAMPLES), "samples.\n\n")


################################################################################
# 6. BUILD FEATURE MATRICES
#
# FeatureMatrix() counts ATAC fragment overlaps at each consensus peak for
# each barcode. This is the equivalent of the gene expression count matrix
# for the ATAC modality.
################################################################################

cat("=== Building feature matrices (this may take several minutes) ===\n")

counts.list <- mapply(
  function(frags, md) {
    FeatureMatrix(
      fragments = frags,
      features  = combined.peaks,
      cells     = rownames(md)
    )
  },
  frags.list, md.list,
  SIMPLIFY = FALSE
)

cat("Feature matrices built.\n\n")


################################################################################
# 7. CREATE PER-SAMPLE CHROMATIN ASSAY SEURAT OBJECTS
#
# For each sample:
#   1. CreateChromatinAssay() wraps the count matrix and fragment object
#      into a Signac-compatible ATAC assay
#   2. CreateSeuratObject() creates the Seurat container with the per-barcode
#      metrics from CellRanger ARC as initial metadata
################################################################################

cat("=== Creating per-sample Seurat objects (ATAC) ===\n")

seurat.list <- mapply(
  function(s, counts, frags, md) {
    assay <- CreateChromatinAssay(
      counts    = counts,
      fragments = frags
    )
    obj <- CreateSeuratObject(
      counts    = assay,
      assay     = "ATAC",
      meta.data = md
    )
    cat("  Sample", s, "— Cells:", ncol(obj), "\n")
    obj
  },
  SAMPLES, counts.list, frags.list, md.list,
  SIMPLIFY = FALSE
)

cat("\n")


################################################################################
# 8. MERGE ALL SAMPLES
#
# Cell barcodes are prefixed with sample IDs (e.g., "S1_ATCG...") to prevent
# barcode collision across samples.
################################################################################

cat("=== Merging all ATAC samples ===\n")

combined <- merge(
  seurat.list[[1]],
  y           = seurat.list[-1],
  add.cell.ids = SAMPLES
)

cat("Merged ATAC object — Total cells:", ncol(combined), "\n\n")


################################################################################
# 9. COMPUTE ATAC QC METRICS
#
# Three key QC metrics are computed:
#
#   NucleosomeSignal()  : Ratio of mono-nucleosomal to sub-nucleosomal
#                         fragments. High signal (>4) indicates poor
#                         transposase access (low-quality nuclei).
#
#   TSSEnrichment()     : Ratio of fragments at transcription start sites
#                         vs. flanking regions. Low enrichment (<2) indicates
#                         poor chromatin accessibility signal.
#
#   pct_reads_in_peaks  : Fraction of ATAC reads overlapping called peaks.
#                         Low values (<15%) indicate high background noise.
################################################################################

cat("=== Computing ATAC QC metrics ===\n")
cat("  Running NucleosomeSignal()...\n")
combined <- NucleosomeSignal(object = combined)

cat("  Running TSSEnrichment()...\n")
combined <- TSSEnrichment(object = combined)

# Percent reads in peaks (computed from CellRanger-reported per-barcode metrics)
combined$pct_reads_in_peaks <- combined$atac_peak_region_fragments /
                               combined$atac_fragments * 100

cat("QC metrics computed.\n\n")

# QC metric summaries
cat("QC metric distributions (all cells, pre-filter):\n")
cat("  nCount_ATAC — median:",
    round(median(combined$nCount_ATAC), 0),
    "| range:", round(min(combined$nCount_ATAC), 0),
    "–", round(max(combined$nCount_ATAC), 0), "\n")
cat("  TSS enrichment — median:",
    round(median(combined$TSS.enrichment), 2), "\n")
cat("  Nucleosome signal — median:",
    round(median(combined$nucleosome_signal), 2), "\n")
cat("  Pct reads in peaks — median:",
    round(median(combined$pct_reads_in_peaks), 1), "%\n\n")


################################################################################
# 10. QC VISUALIZATION (PRE-FILTER)
#
# Standard Signac QC plots:
#   - VlnPlot: per-cell distributions of all four QC metrics
#   - FragmentHistogram: nucleosome banding pattern (should show
#     clear mono- and di-nucleosomal bands)
#   - TSSPlot: TSS enrichment profile (should show sharp peak at TSS)
################################################################################

cat("=== Generating pre-filter ATAC QC plots ===\n")

# Violin plots of all four QC metrics
p_vln_pre <- VlnPlot(
  combined,
  features = c("nCount_ATAC", "TSS.enrichment",
               "nucleosome_signal", "pct_reads_in_peaks"),
  ncol     = 4,
  pt.size  = 0
) &
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    plot.title   = element_text(size = 11, face = "bold")
  )

ggsave(
  filename = file.path(OUT_DIR, "ATAC_QC_VlnPlot_prefilter.pdf"),
  plot     = p_vln_pre,
  width    = 16, height = 5
)

# Nucleosome banding pattern
p_nuc <- FragmentHistogram(
  object = combined,
  group.by = "dataset"
)
ggsave(
  filename = file.path(OUT_DIR, "ATAC_QC_NucleosomeHistogram_prefilter.pdf"),
  plot     = p_nuc,
  width    = 10, height = 5
)

# TSS enrichment profile
p_tss <- TSSPlot(combined, group.by = "dataset") + NoLegend()
ggsave(
  filename = file.path(OUT_DIR, "ATAC_QC_TSSPlot_prefilter.pdf"),
  plot     = p_tss,
  width    = 8, height = 5
)

cat("Pre-filter QC plots saved to:", OUT_DIR, "\n\n")


################################################################################
# 11. QC FILTERING
#
# Applies thresholds to remove low-quality nuclei. Thresholds as reported
# in Methods:
#
#   nCount_ATAC        : 100 – 100,000  (removes empty/multiplet droplets)
#   pct_reads_in_peaks : > 15%          (removes high-background nuclei)
#   nucleosome_signal  : < 4            (removes degraded chromatin)
#   TSS.enrichment     : > 2            (removes low-quality ATAC signal)
################################################################################

cat("=== Applying ATAC QC filters ===\n")
cat("  Pre-filter cells:", ncol(combined), "\n")

combined_filt <- subset(
  x      = combined,
  subset = nCount_ATAC       > MIN_ATAC_COUNTS     &
           nCount_ATAC       < MAX_ATAC_COUNTS     &
           pct_reads_in_peaks > MIN_PCT_READS_PEAKS &
           nucleosome_signal  < MAX_NUCLEOSOME_SIG  &
           TSS.enrichment     > MIN_TSS_ENRICHMENT
)

cat("  Post-filter cells:", ncol(combined_filt), "\n")
cat("  Cells removed    :", ncol(combined) - ncol(combined_filt), "\n")
cat("  Retention rate   :",
    round(ncol(combined_filt) / ncol(combined) * 100, 1), "%\n\n")

# Per-sample cell counts after filtering
cat("Post-filter cell counts per sample:\n")
print(table(combined_filt$dataset))
cat("\nPost-filter cell counts per condition:\n")
print(table(combined_filt$Condition))
cat("\n")


################################################################################
# 12. POST-FILTER QC VISUALIZATION
################################################################################

cat("=== Generating post-filter ATAC QC plots ===\n")

p_vln_post <- VlnPlot(
  combined_filt,
  features = c("nCount_ATAC", "TSS.enrichment",
               "nucleosome_signal", "pct_reads_in_peaks"),
  ncol     = 4,
  pt.size  = 0
) &
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    plot.title   = element_text(size = 11, face = "bold")
  )

ggsave(
  filename = file.path(OUT_DIR, "ATAC_QC_VlnPlot_postfilter.pdf"),
  plot     = p_vln_post,
  width    = 16, height = 5
)

cat("Post-filter QC plots saved.\n\n")


################################################################################
# 13. TF-IDF NORMALIZATION
#
# Term Frequency – Inverse Document Frequency (TF-IDF) normalization is the
# standard normalization for ATAC count matrices. It down-weights peaks
# accessible in many cells (high document frequency) and up-weights peaks
# that are rare and cell-type-specific. Applied to the QC-filtered object.
################################################################################

cat("=== Running TF-IDF normalization ===\n")

combined_filt <- RunTFIDF(combined_filt)

cat("TF-IDF normalization complete.\n\n")


################################################################################
# 14. TOP FEATURE SELECTION
#
# Selects the most informative peaks for LSI dimensionality reduction by
# retaining only peaks present in a minimum number of cells.
#   min.cutoff = 20 : peak must be accessible in at least 20 cells
################################################################################

cat("=== Selecting top features (min.cutoff = 20) ===\n")

combined_filt <- FindTopFeatures(combined_filt, min.cutoff = TFIDF_MIN_CUTOFF)

n_var_features <- length(VariableFeatures(combined_filt))
cat("Variable features selected:", n_var_features, "\n\n")


################################################################################
# 15. LSI DIMENSIONALITY REDUCTION (SVD)
#
# Singular Value Decomposition (SVD) on the TF-IDF matrix produces Latent
# Semantic Indexing (LSI) components. LSI is the ATAC-seq equivalent of PCA.
#
# NOTE: LSI component 1 is typically correlated with sequencing depth and is
# excluded from downstream analyses (dims = 2:50 used in Script 03).
################################################################################

cat("=== Running SVD (LSI) ===\n")

combined_filt <- RunSVD(combined_filt)

cat("SVD complete.\n")
cat("LSI components computed:", length(combined_filt@reductions$lsi), "\n")
cat("NOTE: Component 1 is correlated with sequencing depth and will be\n")
cat("      excluded from UMAP/integration in Script 03 (dims = 2:50).\n\n")


################################################################################
# 16. LSI DEPTH CORRELATION CHECK
#
# Confirms that LSI component 1 correlates with sequencing depth (expected)
# and that subsequent components capture biological variation.
################################################################################

cat("=== Checking LSI component depth correlation ===\n")

p_lsi_depth <- DepthCor(combined_filt)
ggsave(
  filename = file.path(OUT_DIR, "ATAC_LSI_DepthCorrelation.pdf"),
  plot     = p_lsi_depth,
  width    = 6, height = 5
)

cat("LSI depth correlation plot saved.\n\n")


################################################################################
# 17. SUMMARY STATISTICS
################################################################################

cat("========================================\n")
cat("SCRIPT 02 COMPLETE — ATAC QC SUMMARY\n")
cat("========================================\n")
cat("Samples processed         :", length(SAMPLES), "\n")
cat("Consensus peaks           :", length(combined.peaks), "\n")
cat("Pre-filter nuclei         :", ncol(combined), "\n")
cat("Post-filter nuclei        :", ncol(combined_filt), "\n")
cat("Nuclei retained           :",
    round(ncol(combined_filt) / ncol(combined) * 100, 1), "%\n")
cat("Variable features (LSI)   :", length(VariableFeatures(combined_filt)), "\n")
cat("\nQC thresholds applied:\n")
cat("  nCount_ATAC         :", MIN_ATAC_COUNTS, "–", MAX_ATAC_COUNTS, "\n")
cat("  pct_reads_in_peaks  : >", MIN_PCT_READS_PEAKS, "%\n")
cat("  nucleosome_signal   : <", MAX_NUCLEOSOME_SIG, "\n")
cat("  TSS.enrichment      : >", MIN_TSS_ENRICHMENT, "\n")
cat("  Pre-filter fragments: >", MIN_FRAGMENTS, "\n")
cat("  Peak width          :", MIN_PEAK_WIDTH, "–", MAX_PEAK_WIDTH, "bp\n")
cat("\nCondition breakdown:\n")
print(table(combined_filt$Condition, combined_filt$Sex))
cat("\n")


################################################################################
# 18. SAVE OUTPUT
#
# Saves the filtered, TF-IDF normalized, LSI-reduced Signac/Seurat object
# for use in:
#   03_WNN_Integration_and_Clustering.R
#     - Requires: combined_filt with lsi reduction (dims 2:50)
#     - Also requires: combined_filt[["lsi"]] for IntegrateEmbeddings()
################################################################################

cat("=== Saving output object ===\n")

saveRDS(
  combined_filt,
  file = file.path(OUT_DIR, "02_combined_atac_filt.rds")
)

cat("Saved: 02_combined_atac_filt.rds\n")
cat("Next:  Run 03_WNN_Integration_and_Clustering.R\n")
cat("       (requires both 01_merged_seurat_filt.rds and",
    "02_combined_atac_filt.rds)\n\n")


