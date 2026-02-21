################################################################################
# Script:  07_Peak_Annotation.R
# Project: Heat_Call_Single_Nucleus_Multiome
# Author:  Prakrit Subba
# Date:    2026
#
# Paper:   Subba et al. 2026 — Multiome Profiling Reveals Astrocyte and
#          Neuroendocrine Targets of Prenatal Acoustic Programming in
#          Zebra Finch Embryos
#
# Description:
#   Annotates significant DA peaks from Script 06 to nearest genes using
#   multiple complementary approaches:
#     (1) ClosestFeature() — Signac nearest-gene assignment from GTF
#     (2) StringToGRanges → export.bed → bedtools intersect pipeline
#         (for direct gene body overlaps)
#     (3) Merge both annotation sources into a single final table
#     (4) Coverage plots for top HC-gained (top 20) and all HC-lost peaks
#
# Inputs:
#   - results/06_DA/06_combined_wnn_DA.rds       (with Annotation + RegionStats)
#   - results/06_DA/DAHeatCallvsControl.csv       (filtered DA peaks)
#   - results/06_DA/DAAstroDAHCgained.csv
#   - results/06_DA/DAAstroDAHClost.csv
#   - results/06_DA/da_cell_cond_peak_gene_overlaps.csv
#       (produced externally by bedtools intersect on da_cell_cond_peaks.bed)
#
# Outputs:
#   Annotation tables:
#     da_cell_cond_peaks.bed                      (BED for bedtools; underscores)
#     DAHeatCallvsControlAnnotated.csv            (final merged annotation)
#     DAAstroDAgainedclosestGene.csv              (HC-gained nearest genes)
#     DAAstroDAlostclosestGene.csv                (HC-lost nearest genes)
#
#   Coverage plots:
#     FIGURES/AstroTop20HCgainedCoveragePlots.pdf (top 20 gained, w/ extend)
#     FIGURES/AstroHClostCoveragePlots.pdf        (all lost peaks)
#
# NOTE on file naming:
#   - da_cell_cond_peaks.bed           : underscores (intermediate file)
#   - da_cell_cond_peak_gene_overlaps  : underscores (bedtools output)
#   - DAHeatCallvsControlAnnotated.csv : NO underscores between words
#   - DAAstroDAgainedclosestGene.csv   : NO underscores between words
#   - DAAstroDAlostclosestGene.csv     : NO underscores between words
#
# NOTE on peak naming:
#   Signac stores peaks as "chr-start-end" (e.g., NC_044211.2-12345-67890).
#   ClosestFeature requires the format "chr_start-end" — gsub replaces the
#   FIRST hyphen (chr-start separator) with underscore.
#   peak_to_gr() converts peak strings back to GRanges for CoveragePlot,
#   converting any remaining dashes in chromosome names to underscores
#   (CRITICAL: rownames use "-", GRanges use "_").
#
# Bedtools command (run on HPC before Section 7):
#   bedtools intersect \
#     -a results/07_Peak_Annotation/da_cell_cond_peaks.bed \
#     -b data/gene_annotation.bed \
#     -wa -wb > results/07_Peak_Annotation/da_cell_cond_peak_gene_overlaps.csv
#
# Previous Script: 06_Differential_Chromatin_Accessibility.R
# Next Script:     08_ChromVAR_TF_Activity.R
################################################################################


################################################################################
# 0. CONFIGURATION
################################################################################

WNN_RDS         <- "results/06_DA/06_combined_wnn_DA.rds"
DA_RESULTS_CSV  <- "results/06_DA/DAHeatCallvsControl.csv"
DA_GAINED_CSV   <- "results/06_DA/DAAstroDAHCgained.csv"
DA_LOST_CSV     <- "results/06_DA/DAAstroDAHClost.csv"

OUT_DIR         <- "results/07_Peak_Annotation"
dir.create(OUT_DIR,                    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "FIGURES"), showWarnings = FALSE)

# Bedtools output (produced externally; script will use it if present)
BEDTOOLS_OVERLAPS <- file.path(OUT_DIR, "da_cell_cond_peak_gene_overlaps.csv")

# Significance threshold (matches Script 06)
PADJ_CUTOFF <- 0.05

set.seed(1234)


################################################################################
# 1. LOAD LIBRARIES
################################################################################

library(Seurat)
library(Signac)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome)
library(Biostrings)
library(GenomeInfoDb)
library(Rsamtools)
library(ggrepel)
library(dplyr)
library(ggplot2)


################################################################################
# 2. LOAD OBJECTS AND DA RESULTS
################################################################################

cat("=== Loading objects and DA results ===\n")

combined_wnn <- readRDS(WNN_RDS)
DefaultAssay(combined_wnn) <- "ATAC"

cat("combined_wnn — cells:", ncol(combined_wnn),
    "| ATAC peaks:", nrow(combined_wnn@assays$ATAC), "\n")

# Validate that Annotation was added in Script 06
if (is.null(Annotation(combined_wnn))) {
  stop(
    "Annotation(combined_wnn) is NULL.\n",
    "Run 06_Differential_Chromatin_Accessibility.R first to add GTF annotation."
  )
}
cat("Annotation present:", length(Annotation(combined_wnn)), "features\n\n")

# Load significant DA peaks from Script 06
da_cell_cond           <- read.csv(DA_RESULTS_CSV)
colnames(da_cell_cond) <- c("peak", "pval", "avg_log2FC", "pct.1", "pct.2", "pval_adj")

cat("DA peaks loaded:", nrow(da_cell_cond), "\n\n")

# Load astrocyte HC-gained and HC-lost peaks
dagained <- read.csv(DA_GAINED_CSV)
dalost   <- read.csv(DA_LOST_CSV)

cat("Astrocyte HC-gained peaks:", nrow(dagained), "\n")
cat("Astrocyte HC-lost peaks  :", nrow(dalost),   "\n\n")


################################################################################
# 3. PEAK NAME FORMAT FIX FOR CLOSESTFEATURE
#
# ClosestFeature() requires peak names in the format "chr_start-end"
# (underscore between chromosome name and start coordinate).
# Signac rownames use "chr-start-end" (hyphen throughout).
#
# gsub replaces only the FIRST hyphen in each peak name with underscore:
#   NC_044211.2-12345-67890  →  NC_044211.2_12345-67890
################################################################################

cat("=== Adjusting peak name format for ClosestFeature ===\n")
cat("  Before: ", head(da_cell_cond$peak, 2), "\n")

da_cell_cond$peak <- gsub("^(\\S+)-", "\\1_", da_cell_cond$peak)

cat("  After : ", head(da_cell_cond$peak, 2), "\n")

# Filter to padj < 0.05 (should already be filtered, but enforce here)
da_cell_cond <- da_cell_cond[da_cell_cond$pval_adj < PADJ_CUTOFF, ]

cat("Peaks after format fix and filter:", nrow(da_cell_cond), "\n\n")


################################################################################
# 4. PEAK-TO-GRANGES HELPER FUNCTION
#
# Converts a Signac peak string to a single GRanges object.
# Handles both 3-part (standard) and 4-part (chromosome name contains hyphen)
# peak formats from the zebra finch RefSeq genome.
#
# CRITICAL: chromosome names use "-" in peak rownames but "_" in GRanges.
#           gsub("-", "_", chr) performs this conversion.
#
# Examples:
#   "NC_044211.2-12345-67890"        → 3-part → chr = "NC_044211.2"
#   "NC-044211.2-12345-67890"        → 4-part → chr = "NC-044211.2" → "NC_044211.2"
################################################################################

peak_to_gr <- function(pk) {
  parts <- strsplit(pk, "-")[[1]]
  if (length(parts) == 4) {
    # Chromosome name contains a hyphen (e.g., "NC-044211.2")
    chr   <- paste0(parts[1], "-", parts[2])
    start <- as.integer(parts[3])
    end   <- as.integer(parts[4])
  } else if (length(parts) == 3) {
    chr   <- parts[1]
    start <- as.integer(parts[2])
    end   <- as.integer(parts[3])
  } else {
    stop("Bad peak format: ", pk)
  }
  # CRITICAL: rownames use "-", GRanges use "_"
  chr <- gsub("-", "_", chr)
  GRanges(
    seqnames = chr,
    ranges   = IRanges(start = start, end = end)
  )
}


################################################################################
# 5. CLOSESTFEATURE ANNOTATION (GENOME-WIDE DA PEAKS)
#
# Maps each significant DA peak to its nearest gene using the GTF annotation
# added to combined_wnn in Script 06. Returns gene_id, gene_name, distance,
# type (promoter/distal), and strand for each peak.
################################################################################

cat("=== Running ClosestFeature on significant DA peaks ===\n")

annot_da_cell_cond <- ClosestFeature(
  combined_wnn,
  regions = da_cell_cond$peak
)

cat("Peaks annotated by ClosestFeature:", nrow(annot_da_cell_cond), "\n\n")


################################################################################
# 6. EXPORT DA PEAKS TO BED FOR BEDTOOLS ANNOTATION
#
# Converts significant DA peaks to GRanges and exports as BED file for
# external gene-body overlap annotation using bedtools intersect.
# Run the bedtools command on HPC (see header) before Section 7.
################################################################################

cat("=== Exporting DA peaks to BED for bedtools ===\n")

peaks_gr <- StringToGRanges(da_cell_cond$peak)

export.bed(
  peaks_gr,
  con = file.path(OUT_DIR, "da_cell_cond_peaks.bed")
)

cat("Saved: da_cell_cond_peaks.bed (", length(peaks_gr), "peaks)\n\n")


################################################################################
# 7. MERGE BEDTOOLS OVERLAPS WITH CLOSESTFEATURE ANNOTATION
#
# Reads the bedtools gene-body overlap output and merges it with the
# ClosestFeature results to capture both:
#   (1) Peaks directly overlapping gene bodies (bedtools)
#   (2) Peaks nearest to a gene (ClosestFeature)
#
# Peak ID is reconstructed from BED columns:
#   paste0(chrom, "-", as.integer(start) + 1, "-", end)
#   (BED is 0-based; adding 1 to start converts to 1-based peak ID)
#
# Final table: da_cell_cond_gene + matching_peaks → final_df
# Removes duplicates and rows with no gene annotation.
################################################################################

if (file.exists(BEDTOOLS_OVERLAPS)) {

  cat("=== Merging bedtools peak-gene overlaps ===\n")

  peak_gene_overlaps <- read.csv(file = BEDTOOLS_OVERLAPS)

  # Reconstruct peak ID from BED columns (BED is 0-based → add 1 to start)
  peak_gene_overlaps <- peak_gene_overlaps %>%
    mutate(peak = paste0(chrom, "-", as.integer(start) + 1, "-", end)) %>%
    select(peak, gene_id) %>%
    distinct()

  cat("Unique peak-gene pairs from bedtools:", nrow(peak_gene_overlaps), "\n")

  # Merge bedtools gene_id with DA results
  da_cell_cond_gene <- merge(
    da_cell_cond,
    peak_gene_overlaps,
    by    = "peak",
    all.x = TRUE
  )

  # Pull relevant columns from ClosestFeature for secondary merge
  annot_subset            <- annot_da_cell_cond[, c("query_region", "gene_id")]
  colnames(annot_subset)[1] <- "peak"

  # Merge: update gene_id from ClosestFeature where available
  matching_peaks         <- merge(da_cell_cond_gene, annot_subset, by = "peak")
  matching_peaks$gene_id <- matching_peaks$gene_id.y

  # Combine all annotated results; remove duplicates and unannotated peaks
  final_df <- rbind(da_cell_cond_gene, matching_peaks)
  final_df <- unique(final_df)
  final_df <- final_df[!is.na(final_df$gene_id), ]

  cat("Final annotated DA peaks:", nrow(final_df), "\n")

} else {

  cat("NOTE: da_cell_cond_peak_gene_overlaps.csv not found.\n")
  cat("      Proceeding with ClosestFeature annotation only.\n")
  cat("      Run bedtools intersect on HPC, then re-run Section 7.\n\n")

  # Fallback: use ClosestFeature results only
  annot_subset            <- annot_da_cell_cond[, c("query_region", "gene_id")]
  colnames(annot_subset)[1] <- "peak"
  final_df                <- merge(da_cell_cond, annot_subset, by = "peak", all.x = TRUE)
  final_df                <- final_df[!is.na(final_df$gene_id), ]
  cat("ClosestFeature-only annotated peaks:", nrow(final_df), "\n")
}

write.csv(
  final_df,
  file      = file.path(OUT_DIR, "DAHeatCallvsControlAnnotated.csv"),
  row.names = FALSE
)

cat("Saved: DAHeatCallvsControlAnnotated.csv\n\n")


################################################################################
# 8. CLOSESTFEATURE ANNOTATION OF ASTROCYTE DA PEAKS
#
# Annotates astrocyte-specific HC-gained and HC-lost peaks to nearest genes.
# Uses the peak_to_gr() helper to convert peak strings to GRanges.
#
# Output files (NO underscores between words — exact names required by
# astrocyte rewiring scripts):
#   DAAstroDAgainedclosestGene.csv
#   DAAstroDAlostclosestGene.csv
################################################################################

cat("=== Annotating astrocyte DA peaks to nearest genes ===\n")

# Subset to astrocyte object for ClosestFeature context
astroobj           <- subset(combined_wnn, subset = cluster_annotation == "Astrocyte")
DefaultAssay(astroobj) <- "ATAC"

cat("Astrocyte object — cells:", ncol(astroobj), "\n\n")

# HC-gained peaks
if (nrow(dagained) > 0) {
  gr_gained        <- do.call(c, lapply(dagained$peak, peak_to_gr))
  names(gr_gained) <- dagained$peak

  closest_gained <- ClosestFeature(astroobj, regions = gr_gained)

  write.csv(
    closest_gained,
    file.path(OUT_DIR, "DAAstroDAgainedclosestGene.csv"),
    row.names = FALSE
  )
  cat("HC-gained annotated:", nrow(closest_gained), "peaks\n")
  cat("Saved: DAAstroDAgainedclosestGene.csv\n\n")
}

# HC-lost peaks
if (nrow(dalost) > 0) {
  gr_lost        <- do.call(c, lapply(dalost$peak, peak_to_gr))
  names(gr_lost) <- dalost$peak

  closest_lost <- ClosestFeature(astroobj, regions = gr_lost)

  write.csv(
    closest_lost,
    file.path(OUT_DIR, "DAAstroDAlostclosestGene.csv"),
    row.names = FALSE
  )
  cat("HC-lost annotated:", nrow(closest_lost), "peaks\n")
  cat("Saved: DAAstroDAlostclosestGene.csv\n\n")
}


################################################################################
# 9. COVERAGE PLOTS — HC-GAINED PEAKS (TOP 20)
#
# CoveragePlot parameters:
#   split.by          = "Condition"
#   peaks             = TRUE
#   annotation        = FALSE
#   extend.upstream   = 2000   (2 kb extension for gained peaks)
#   extend.downstream = 2000
#   scale_fill_manual : Control = "#0072B2", Heat Call = "#D55E00"
#
# Title format: paste0(i, ". ", pk, " log2FC: ", round(avg_log2FC, 2))
# Output: FIGURES/AstroTop20HCgainedCoveragePlots.pdf (14 × 4 inches)
################################################################################

cat("=== Coverage plots — HC-gained (top 20) ===\n")

Idents(astroobj) <- "Condition"

if (nrow(dagained) > 0) {

  pdf(
    file.path(OUT_DIR, "FIGURES", "AstroTop20HCgainedCoveragePlots.pdf"),
    width = 14, height = 4
  )

  for (i in 1:min(20, nrow(dagained))) {
    pk <- dagained$peak[i]
    gr <- peak_to_gr(pk)

    p <- try(
      CoveragePlot(
        astroobj,
        region            = gr,
        split.by          = "Condition",
        peaks             = TRUE,
        annotation        = FALSE,
        extend.upstream   = 2000,
        extend.downstream = 2000
      ) +
        scale_fill_manual(
          values = c("Control" = "#0072B2", "Heat Call" = "#D55E00")
        ) +
        ggtitle(
          paste0(i, ". ", pk, " log2FC: ", round(dagained$avg_log2FC[i], 2))
        ),
      silent = TRUE
    )

    if (!inherits(p, "try-error")) print(p)
  }

  dev.off()
  cat("Plotted top 20 HC-gained peaks.\n")
  cat("Saved: FIGURES/AstroTop20HCgainedCoveragePlots.pdf\n\n")

} else {
  cat("No HC-gained peaks — skipping coverage plots.\n\n")
}


################################################################################
# 10. COVERAGE PLOTS — HC-LOST PEAKS (ALL)
#
# CoveragePlot parameters (same as above minus extension):
#   split.by   = "Condition"
#   peaks      = TRUE
#   annotation = FALSE
#   (no extend.upstream / extend.downstream for lost peaks)
#   scale_fill_manual: Control = "#0072B2", Heat Call = "#D55E00"
#
# Title: peak name only (pk)
# Output: FIGURES/AstroHClostCoveragePlots.pdf (14 × 4 inches)
################################################################################

cat("=== Coverage plots — HC-lost (all", nrow(dalost), "peaks) ===\n")

if (nrow(dalost) > 0) {

  pdf(
    file.path(OUT_DIR, "FIGURES", "AstroHClostCoveragePlots.pdf"),
    width = 14, height = 4
  )

  for (i in 1:nrow(dalost)) {
    pk <- dalost$peak[i]
    gr <- peak_to_gr(pk)

    p <- try(
      CoveragePlot(
        astroobj,
        region     = gr,
        split.by   = "Condition",
        peaks      = TRUE,
        annotation = FALSE
      ) +
        scale_fill_manual(
          values = c("Control" = "#0072B2", "Heat Call" = "#D55E00")
        ) +
        ggtitle(pk),
      silent = TRUE
    )

    if (!inherits(p, "try-error")) print(p)
  }

  dev.off()
  cat("Plotted", nrow(dalost), "HC-lost peaks.\n")
  cat("Saved: FIGURES/AstroHClostCoveragePlots.pdf\n\n")

} else {
  cat("No HC-lost peaks — skipping coverage plots.\n\n")
}


################################################################################
# 11. SUMMARY
################################################################################

cat("========================================\n")
cat("SCRIPT 07 COMPLETE — ANNOTATION SUMMARY\n")
cat("========================================\n")
cat("Genome-wide DA peaks annotated  :", nrow(da_cell_cond), "\n")
cat("Final annotated table rows      :", nrow(final_df), "\n")
cat("Astrocyte HC-gained annotated   :",
    if (nrow(dagained) > 0) nrow(closest_gained) else 0, "\n")
cat("Astrocyte HC-lost annotated     :",
    if (nrow(dalost)   > 0) nrow(closest_lost)   else 0, "\n")
cat("\nOutput files:\n")
cat("  da_cell_cond_peaks.bed\n")
cat("  DAHeatCallvsControlAnnotated.csv\n")
cat("  DAAstroDAgainedclosestGene.csv\n")
cat("  DAAstroDAlostclosestGene.csv\n")
cat("  FIGURES/AstroTop20HCgainedCoveragePlots.pdf\n")
cat("  FIGURES/AstroHClostCoveragePlots.pdf\n\n")
cat("Next: Run 08_ChromVAR_TF_Activity.R\n\n")


