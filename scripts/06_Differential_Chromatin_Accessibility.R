################################################################################
# Script:  06_Differential_Chromatin_Accessibility.R
# Project: Heat_Call_Single_Nucleus_Multiome
# Author:  Prakrit Subba
# Date:    2026
#
# Paper:   Subba et al. 2026 — Multiome Profiling Reveals Astrocyte and
#          Neuroendocrine Targets of Prenatal Acoustic Programming in
#          Zebra Finch Embryos
#
# Description:
#   Genome-wide differential chromatin accessibility (DA) analysis across all
#   49 WNN clusters, plus astrocyte-specific DA analysis. Also adds genome
#   annotation and GC content to the Signac object (required for ChromVAR
#   in Script 07 and LinkPeaks in the astrocyte rewiring scripts).
#
#   Analyses:
#     (1) GTF gene annotation → Annotation(combined_wnn)
#     (2) RegionStats GC content → fasta reference
#     (3) All-cluster DA: Heat Call vs Control (MAST, per WNN cluster)
#     (4) Combine, filter, and annotate DA peaks (ClosestFeature + bedtools)
#     (5) Astrocyte-specific DA: Heat Call vs Control
#
# Inputs:
#   - results/05_DEG/05_combined_wnn_DEG.rds
#   - data/genomic.gtf                                  (zebra finch GTF)
#   - data/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna  (zebra finch genome)
#
# Outputs (exact file names required by downstream scripts):
#
#   Genome-wide DA:
#     Cluster_<N>_DA.csv                          (one per cluster)
#     AllClustersHeatCallvsControlDA.csv           (combined raw)
#     DAHeatCallvsControl.csv                      (filtered padj<0.05;
#                                                   required by Scripts 09+)
#     da_cell_cond_peaks.bed                       (for bedtools, underscores)
#     da_cell_cond_peak_gene_overlaps.csv          (bedtools output, underscores)
#     DAHeatCallvsControlAnnotated.csv             (final annotated table)
#
#   Astrocyte-specific DA:
#     DAAstrocyteDAMASTall.csv                     (all peaks tested)
#     DAAstroDAsig.csv                             (FDR<0.05, |LFC|>0.25)
#     DAAstroDAHCgained.csv                        (HC-gained peaks)
#     DAAstroDAHClost.csv                          (HC-lost peaks)
#     DAAstroDAgainedclosestGene.csv               (HC-gained annotated)
#     DAAstroDAlostclosestGene.csv                 (HC-lost annotated)
#
#   Object:
#     06_combined_wnn_DA.rds                       (with Annotation + RegionStats)
#
# NOTE on file naming:
#   - Cluster_<N>_DA.csv       : underscores around cluster number
#   - AllClusters...DA.csv     : NO underscores between words
#   - DAHeatCallvsControl.csv  : NO underscores (exact name for downstream)
#   - da_cell_cond_*.csv/.bed  : underscores (intermediate files)
#   - DAAstro*.csv             : NO underscores between prefix and words
#
# Column names in DAHeatCallvsControl.csv (set by downstream read code):
#   peak, pval, avg_log2FC, pct.1, pct.2, pval_adj
#
# Previous Script: 05_Differential_Gene_Expression.R
# Next Script:     07_ChromVAR_TF_Activity.R
################################################################################


################################################################################
# 0. CONFIGURATION
################################################################################

WNN_RDS  <- "results/05_DEG/05_combined_wnn_DEG.rds"

# Reference genome files (zebra finch GCF_003957565.2)
GTF_PATH  <- "data/genomic.gtf"
FASTA_PATH <- "data/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna"

OUT_DIR   <- "results/06_DA"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# DA thresholds
DA_MIN_PCT        <- 0.1
DA_TEST           <- "MAST"
PADJ_CUTOFF       <- 0.05     # for DAHeatCallvsControl.csv (genome-wide)
FDR_THRESHOLD     <- 0.05     # for astrocyte-specific DA
LFC_THRESHOLD     <- 0.25     # for astrocyte-specific DA (scATAC standard)

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
# 2. LOAD ANNOTATED WNN OBJECT
################################################################################

cat("=== Loading annotated WNN object ===\n")

combined_wnn <- readRDS(WNN_RDS)

cat("Cells      :", ncol(combined_wnn), "\n")
cat("ATAC peaks :", nrow(combined_wnn@assays$ATAC), "\n")
cat("Conditions :", paste(unique(combined_wnn$Condition), collapse = ", "), "\n\n")

stopifnot(
  "ATAC"             %in% names(combined_wnn@assays),
  "wsnn_res.0.6"     %in% colnames(combined_wnn@meta.data),
  "Condition"        %in% colnames(combined_wnn@meta.data),
  "cluster_annotation" %in% colnames(combined_wnn@meta.data)
)


################################################################################
# 3. ADD GENE ANNOTATION FROM GTF
#
# Imports the zebra finch GTF and filters to feature types compatible with
# Signac's Annotation slot. Required for:
#   - ClosestFeature() (peak-to-gene annotation)
#   - LinkPeaks()      (peak-gene correlation; astrocyte rewiring scripts)
#   - TSSEnrichment()  (if re-run)
#
# Key steps:
#   - Rename "transcript_id" column → "tx_id" (Signac expected name)
#   - Add "gene_name" from "gene_id" if absent (zebra finch GTF may lack it)
#   - Filter to chromosomes present in the ATAC assay (coarse pruning)
################################################################################

cat("=== Adding GTF gene annotation ===\n")
cat("  Importing GTF:", GTF_PATH, "\n")

gtf        <- import(GTF_PATH)
gene.coords <- gtf[gtf$type %in%
                   c("gene", "transcript", "exon", "CDS",
                     "start_codon", "stop_codon")]

# Rename transcript_id → tx_id (required by Signac)
names(mcols(gene.coords))[
  names(mcols(gene.coords)) == "transcript_id"
] <- "tx_id"

# Ensure gene_name column exists (zebra finch GTF uses gene_id only)
if (!"gene_name" %in% names(mcols(gene.coords))) {
  gene.coords$gene_name <- gene.coords$gene_id
  cat("  NOTE: gene_name absent in GTF — using gene_id as gene_name\n")
}

# Filter to chromosomes present in the ATAC assay (leave seqnames as RefSeq)
valid_chroms <- seqlevels(granges(combined_wnn@assays$ATAC))
gene.coords  <- keepSeqlevels(
  gene.coords,
  valid_chroms,
  pruning.mode = "coarse"
)

# Add annotation to Signac object
Annotation(combined_wnn) <- gene.coords

cat("  Features in annotation:", length(gene.coords), "\n")
cat("  Chromosomes retained  :", length(valid_chroms), "\n")
cat("  Annotation added to combined_wnn.\n\n")


################################################################################
# 4. COMPUTE GC CONTENT WITH REGIONSTATS
#
# RegionStats() computes GC content for each peak using the reference genome.
# Required by RunChromVAR() (Script 07) for bias correction.
# The FASTA index (.fai) is created automatically if absent.
################################################################################

cat("=== Running RegionStats (GC content) ===\n")
cat("  FASTA:", FASTA_PATH, "\n")

fasta <- FaFile(FASTA_PATH)

# Index if not already done
indexFa(fasta)

combined_wnn <- RegionStats(combined_wnn, genome = fasta)

cat("  RegionStats complete.\n\n")


################################################################################
# 5. SET UP IDENTITIES FOR DA
#
# Identical to DEG setup (Script 05) but on the ATAC assay.
# cluster_condition values: "0_Heat Call", "0_Control", "1_Heat Call", ...
################################################################################

cat("=== Setting up cluster_condition identities (ATAC) ===\n")

DefaultAssay(combined_wnn) <- "ATAC"

combined_wnn$cluster_condition <- paste(
  combined_wnn$wsnn_res.0.6,
  combined_wnn$Condition,
  sep = "_"
)

Idents(combined_wnn) <- "cluster_condition"

cat("Example identities:\n")
print(head(sort(unique(combined_wnn$cluster_condition))))
cat("\n")


################################################################################
# 6. ALL-CLUSTER DA: HEAT CALL VS CONTROL
#
# Parameters (manuscript Methods):
#   ident.1  = paste0(cluster, "_Heat Call")
#   ident.2  = paste0(cluster, "_Control")
#   assay    = "ATAC"
#   test.use = "MAST"
#   min.pct  = 0.1
#
# Adds $Cluster and $Peak columns before saving each per-cluster CSV.
# Per-cluster output: Cluster_<N>_DA.csv  (underscores around cluster number)
################################################################################

cat("=== Running all-cluster DA (Heat Call vs Control) ===\n")
cat("(MAST test on ATAC — run on HPC)\n\n")

clusters  <- levels(combined_wnn$wsnn_res.0.6)
DAresults <- list()

for (cluster in clusters) {

  ident_hc   <- paste0(cluster, "_Heat Call")
  ident_ctrl <- paste0(cluster, "_Control")

  if (!ident_hc   %in% Idents(combined_wnn) ||
      !ident_ctrl %in% Idents(combined_wnn)) {
    cat("  Cluster", cluster, ": skipped (identity not found)\n")
    next
  }

  cat("  Cluster", cluster, "...")

  DAresults[[cluster]] <- FindMarkers(
    combined_wnn,
    ident.1  = ident_hc,
    ident.2  = ident_ctrl,
    assay    = "ATAC",
    test.use = DA_TEST,
    min.pct  = DA_MIN_PCT
  )

  if (nrow(DAresults[[cluster]]) > 0) {
    DAresults[[cluster]]$Cluster <- cluster
    DAresults[[cluster]]$Peak    <- rownames(DAresults[[cluster]])

    write.csv(
      DAresults[[cluster]],
      file = file.path(OUT_DIR, paste0("Cluster_", cluster, "_DA.csv"))
    )
    cat(" done —", nrow(DAresults[[cluster]]), "peaks tested\n")
  } else {
    cat(" no results\n")
  }
}

cat("\nAll-cluster DA complete.\n\n")


################################################################################
# 7. COMBINE AND SAVE ALL-CLUSTER DA
#
# Combined raw output: AllClustersHeatCallvsControlDA.csv
# (clean name for the date-stamped original 29May2025DAconditioncellcluster.csv)
################################################################################

cat("=== Combining all DA results ===\n")

combined_da <- do.call(rbind, DAresults)

write.csv(
  combined_da,
  file      = file.path(OUT_DIR, "AllClustersHeatCallvsControlDA.csv"),
  row.names = FALSE
)

cat("Total peaks tested across all clusters:", nrow(combined_da), "\n")
cat("Saved: AllClustersHeatCallvsControlDA.csv\n\n")


################################################################################
# 8. FILTER TO SIGNIFICANT DA AND SAVE DAHeatCallvsControl.csv
#
# Filters to p_val_adj < 0.05. Saves as DAHeatCallvsControl.csv with column
# order exactly as expected by downstream scripts (Scripts 09+):
#   peak, pval, avg_log2FC, pct.1, pct.2, pval_adj
#
# IMPORTANT: This exact file name (DAHeatCallvsControl.csv) and column order
# are required by all downstream astrocyte analysis scripts.
################################################################################

cat("=== Filtering to significant DA (p_val_adj <", PADJ_CUTOFF, ") ===\n")

da_cell_cond <- combined_da %>%
  filter(p_val_adj < PADJ_CUTOFF) %>%
  select(
    peak       = Peak,
    pval       = p_val,
    avg_log2FC,
    pct.1,
    pct.2,
    pval_adj   = p_val_adj
  )

write.csv(
  da_cell_cond,
  file      = file.path(OUT_DIR, "DAHeatCallvsControl.csv"),
  row.names = FALSE
)

cat("Significant DA peaks :", nrow(da_cell_cond), "\n")
cat("  HC-gained (LFC>0)  :", sum(da_cell_cond$avg_log2FC > 0), "\n")
cat("  HC-lost   (LFC<0)  :", sum(da_cell_cond$avg_log2FC < 0), "\n")
cat("Saved: DAHeatCallvsControl.csv\n\n")

# Summary by cluster
cat("DA peaks per cluster (top 10):\n")
print(
  combined_da %>%
    filter(p_val_adj < PADJ_CUTOFF) %>%
    count(Cluster, name = "n_DA_peaks") %>%
    arrange(desc(n_DA_peaks)) %>%
    head(10)
)
cat("\n")


################################################################################
# 9. PEAK NAME FORMAT FIX FOR CLOSESTFEATURE
#
# Signac stores peaks as "chr-start-end" but some chromosome names contain
# hyphens (e.g., NC_044211.2-12345-67890). ClosestFeature requires the
# chromosome-start separator to be an underscore, not a hyphen.
# The gsub replaces the FIRST hyphen in each peak name with underscore:
#   NC_044211.2-12345-67890 → NC_044211.2_12345-67890
################################################################################

cat("=== Adjusting peak name format for ClosestFeature ===\n")

da_cell_cond$peak <- gsub("^(\\S+)-", "\\1_", da_cell_cond$peak)

cat("Example adjusted peak names:\n")
print(head(da_cell_cond$peak))
cat("\n")


################################################################################
# 10. CLOSESTFEATURE ANNOTATION
#
# Annotates each significant DA peak to its nearest gene using the GTF
# annotation added in Section 3. Returns gene_id, gene_name, distance,
# type (promoter/distal), and strand.
################################################################################

cat("=== Running ClosestFeature annotation ===\n")

annot_da_cell_cond <- ClosestFeature(
  combined_wnn,
  regions = da_cell_cond$peak
)

cat("Peaks annotated:", nrow(annot_da_cell_cond), "\n\n")


################################################################################
# 11. EXPORT BED AND BEDTOOLS PEAK-GENE OVERLAP
#
# Exports significant DA peaks as BED for external bedtools intersect with
# a gene body annotation. The bedtools output (da_cell_cond_peak_gene_overlaps.csv)
# is then read back and merged with the DA results.
#
# NOTE: The bedtools intersect command should be run on the HPC between
#       Sections 11 and 12:
#
#   bedtools intersect \
#     -a da_cell_cond_peaks.bed \
#     -b gene_annotation.bed \
#     -wa -wb > da_cell_cond_peak_gene_overlaps.csv
################################################################################

cat("=== Exporting DA peaks to BED for bedtools annotation ===\n")

peaks_gr <- StringToGRanges(da_cell_cond$peak)

export.bed(
  peaks_gr,
  con = file.path(OUT_DIR, "da_cell_cond_peaks.bed")
)

cat("BED file saved: da_cell_cond_peaks.bed\n")
cat("NOTE: Run bedtools intersect on HPC before Section 12.\n\n")


################################################################################
# 12. MERGE BEDTOOLS PEAK-GENE OVERLAPS WITH DA RESULTS
#
# Reads bedtools output and merges with ClosestFeature results to produce
# a comprehensive DA peak annotation table.
################################################################################

peak_gene_file <- file.path(OUT_DIR, "da_cell_cond_peak_gene_overlaps.csv")

if (file.exists(peak_gene_file)) {

  cat("=== Merging bedtools peak-gene overlaps ===\n")

  peak_gene_overlaps <- read.csv(peak_gene_file)

  # Reconstruct peak identifier from BED columns (chr, start, end)
  peak_gene_overlaps <- peak_gene_overlaps %>%
    mutate(peak = paste0(chrom, "-", as.integer(start) + 1, "-", end)) %>%
    select(peak, gene_id) %>%
    distinct()

  # Merge with DA results
  da_cell_cond_gene <- merge(
    da_cell_cond,
    peak_gene_overlaps,
    by   = "peak",
    all.x = TRUE
  )

  # Pull closest gene annotations from ClosestFeature
  annot_subset <- annot_da_cell_cond[, c("query_region", "gene_id")]
  colnames(annot_subset)[1] <- "peak"

  # Update gene_id from ClosestFeature where available
  matching_peaks <- merge(
    da_cell_cond_gene,
    annot_subset,
    by = "peak"
  )
  matching_peaks$gene_id <- matching_peaks$gene_id.y

  # Combine all annotated peaks
  final_df <- rbind(da_cell_cond_gene, matching_peaks)
  final_df <- unique(final_df)
  final_df <- final_df[!is.na(final_df$gene_id), ]

  write.csv(
    final_df,
    file      = file.path(OUT_DIR, "DAHeatCallvsControlAnnotated.csv"),
    row.names = FALSE
  )

  cat("Final annotated DA peaks:", nrow(final_df), "\n")
  cat("Saved: DAHeatCallvsControlAnnotated.csv\n\n")

} else {
  cat("NOTE: da_cell_cond_peak_gene_overlaps.csv not found.\n")
  cat("      Run bedtools intersect first, then re-run this section.\n")
  cat("      Proceeding with ClosestFeature annotation only.\n\n")

  write.csv(
    annot_da_cell_cond,
    file      = file.path(OUT_DIR, "DAHeatCallvsControlAnnotated.csv"),
    row.names = FALSE
  )
}


################################################################################
# 13. ASTROCYTE-SPECIFIC DA ANALYSIS
#
# Subsets combined_wnn to Astrocyte cluster only, then runs FindMarkers
# comparing Heat Call vs Control directly on the astrocyte subset.
#
# This is separate from the all-cluster loop above — identities here are
# set to Condition directly (not cluster_condition), since all cells in
# the subset are already astrocytes.
#
# Parameters identical to genome-wide analysis:
#   test.use = "MAST"
#   min.pct  = 0.1
#   assay    = "ATAC"
#
# Filtering thresholds (manuscript Methods):
#   FDR_THRESHOLD = 0.05
#   LFC_THRESHOLD = 0.25  (scATAC standard)
#
# Output file names (NO underscores between words — exact names required
# by astrocyte rewiring scripts):
#   DAAstrocyteDAMASTall.csv
#   DAAstroDAsig.csv
#   DAAstroDAHCgained.csv
#   DAAstroDAHClost.csv
#   DAAstroDAgainedclosestGene.csv
#   DAAstroDAlostclosestGene.csv
################################################################################

cat("=== Astrocyte-specific DA analysis ===\n")

# Subset to astrocytes
astroobj <- subset(combined_wnn, subset = cluster_annotation == "Astrocyte")
DefaultAssay(astroobj) <- "ATAC"
Idents(astroobj)       <- "Condition"

cat("Astrocyte cells:", ncol(astroobj), "\n")
print(table(astroobj$Condition))
cat("\n")

cat("  Running DA (MAST)...\n")

DAresults_astro <- FindMarkers(
  astroobj,
  ident.1  = "Heat Call",
  ident.2  = "Control",
  assay    = "ATAC",
  test.use = DA_TEST,
  min.pct  = DA_MIN_PCT
)

DAresults_astro$peak <- rownames(DAresults_astro)

write.csv(
  DAresults_astro,
  file      = file.path(OUT_DIR, "DAAstrocyteDAMASTall.csv"),
  row.names = FALSE
)

cat("  Peaks tested:", nrow(DAresults_astro), "\n")
cat("  Saved: DAAstrocyteDAMASTall.csv\n\n")

# Filter: FDR < 0.05, |LFC| > 0.25
dasig <- DAresults_astro %>%
  filter(p_val_adj < FDR_THRESHOLD, abs(avg_log2FC) > LFC_THRESHOLD)

dagained <- dasig %>% filter(avg_log2FC > 0) %>% arrange(p_val_adj)
dalost   <- dasig %>% filter(avg_log2FC < 0) %>% arrange(p_val_adj)

write.csv(dasig,    file.path(OUT_DIR, "DAAstroDAsig.csv"),      row.names = FALSE)
write.csv(dagained, file.path(OUT_DIR, "DAAstroDAHCgained.csv"), row.names = FALSE)
write.csv(dalost,   file.path(OUT_DIR, "DAAstroDAHClost.csv"),   row.names = FALSE)

cat("=== Astrocyte DA results ===\n")
cat("  Total DA peaks (FDR<", FDR_THRESHOLD, ", |LFC|>", LFC_THRESHOLD, "):", nrow(dasig), "\n")
cat("  HC-gained:", nrow(dagained), "\n")
cat("  HC-lost  :", nrow(dalost),   "\n\n")


################################################################################
# 14. CLOSESTFEATURE ANNOTATION OF ASTROCYTE DA PEAKS
#
# Annotates HC-gained and HC-lost astrocyte peaks to nearest genes.
# Peak-to-gene conversion uses a helper to handle the chromosome-hyphen format.
################################################################################

cat("=== Annotating astrocyte DA peaks to nearest genes ===\n")

# Helper: convert Signac peak string → single GRanges object
peak_to_gr <- function(peak) {
  p <- strsplit(peak, "-")[[1]]
  if (length(p) == 4) {
    # NC-XXX-start-end format (chromosome name has hyphen)
    chr   <- paste(p[1], p[2], sep = "-")
    start <- as.integer(p[3])
    end   <- as.integer(p[4])
  } else if (length(p) == 3) {
    chr   <- p[1]
    start <- as.integer(p[2])
    end   <- as.integer(p[3])
  } else {
    stop("Unexpected peak format: ", peak)
  }
  chr <- gsub("_", "-", chr)   # restore chromosome name hyphens
  GRanges(seqnames = chr, ranges = IRanges(start = start, end = end),
          peak_id  = peak)
}

# Annotate HC-gained peaks
if (nrow(dagained) > 0) {
  gr_gained      <- do.call(c, lapply(dagained$peak, peak_to_gr))
  names(gr_gained) <- dagained$peak

  closest_gained <- ClosestFeature(astroobj, regions = gr_gained)

  write.csv(
    closest_gained,
    file.path(OUT_DIR, "DAAstroDAgainedclosestGene.csv"),
    row.names = FALSE
  )
  cat("  HC-gained annotated:", nrow(closest_gained), "peaks\n")
}

# Annotate HC-lost peaks
if (nrow(dalost) > 0) {
  gr_lost        <- do.call(c, lapply(dalost$peak, peak_to_gr))
  names(gr_lost) <- dalost$peak

  closest_lost   <- ClosestFeature(astroobj, regions = gr_lost)

  write.csv(
    closest_lost,
    file.path(OUT_DIR, "DAAstroDAlostclosestGene.csv"),
    row.names = FALSE
  )
  cat("  HC-lost annotated  :", nrow(closest_lost), "peaks\n")
}

cat("\n")


################################################################################
# 15. ASTROCYTE DA VISUALIZATION
#
# Volcano plot: HC-gained (D55E00) vs HC-lost (0072B2) vs NS (grey)
# MA plot: mean accessibility vs log2FC
################################################################################

cat("=== Generating astrocyte DA figures ===\n")

dir.create(file.path(OUT_DIR, "FIGURES"), showWarnings = FALSE)

DAresults_astro <- DAresults_astro %>%
  mutate(direction = case_when(
    p_val_adj < FDR_THRESHOLD &  avg_log2FC >  LFC_THRESHOLD ~ "HC-gained",
    p_val_adj < FDR_THRESHOLD &  avg_log2FC < -LFC_THRESHOLD ~ "HC-lost",
    TRUE                                                       ~ "NS"
  ))

# Volcano plot
p_volcano <- ggplot(
  DAresults_astro,
  aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)
) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(
    values = c("HC-gained" = "#D55E00", "HC-lost" = "#0072B2", "NS" = "grey"),
    breaks = c("HC-gained", "HC-lost", "NS")
  ) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), lty = 2) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD), lty = 2) +
  labs(
    title = paste0(
      "Astrocyte DA — ", nrow(dagained), " gained, ", nrow(dalost), " lost"
    ),
    x = "log2 Fold Change (Heat Call / Control)",
    y = "-log10(FDR)"
  ) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.85))

ggsave(
  file.path(OUT_DIR, "FIGURES", "AstroDAVolcano.pdf"),
  p_volcano,
  width = 8, height = 6
)

# MA plot
DAresults_astro$baseMean <- (DAresults_astro$pct.1 + DAresults_astro$pct.2) / 2

p_ma <- ggplot(
  DAresults_astro,
  aes(x = baseMean, y = avg_log2FC, color = direction)
) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(
    values = c("HC-gained" = "#D55E00", "HC-lost" = "#0072B2", "NS" = "grey")
  ) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             lty = 2, color = "red", alpha = 0.5) +
  labs(
    title = "Astrocyte DA — MA Plot",
    x     = "Mean Cells Accessible",
    y     = "log2FC"
  ) +
  theme_classic()

ggsave(
  file.path(OUT_DIR, "FIGURES", "AstroDAMAplot.pdf"),
  p_ma,
  width = 8, height = 6
)

cat("Figures saved to FIGURES/\n\n")


################################################################################
# 16. SUMMARY
################################################################################

cat("========================================\n")
cat("SCRIPT 06 COMPLETE — DA SUMMARY\n")
cat("========================================\n")
cat("Genome-wide DA:\n")
cat("  Clusters tested        :", length(DAresults), "\n")
cat("  Total peaks tested     :", nrow(combined_da), "\n")
cat("  Significant (padj<0.05):", nrow(da_cell_cond), "\n")
cat("    HC-gained            :", sum(da_cell_cond$avg_log2FC > 0), "\n")
cat("    HC-lost              :", sum(da_cell_cond$avg_log2FC < 0), "\n")
cat("\nAstrocyte-specific DA:\n")
cat("  Total peaks tested     :", nrow(DAresults_astro), "\n")
cat(sprintf("  Significant (FDR<%.2f, |LFC|>%.2f): %d\n",
            FDR_THRESHOLD, LFC_THRESHOLD, nrow(dasig)))
cat("    HC-gained            :", nrow(dagained), "\n")
cat("    HC-lost              :", nrow(dalost),   "\n")
cat("\nOutput files:\n")
cat("  AllClustersHeatCallvsControlDA.csv\n")
cat("  DAHeatCallvsControl.csv             (required by Scripts 09+)\n")
cat("  Cluster_<N>_DA.csv                  (one per cluster)\n")
cat("  da_cell_cond_peaks.bed\n")
cat("  DAHeatCallvsControlAnnotated.csv\n")
cat("  DAAstrocyteDAMASTall.csv\n")
cat("  DAAstroDAsig.csv\n")
cat("  DAAstroDAHCgained.csv\n")
cat("  DAAstroDAHClost.csv\n")
cat("  DAAstroDAgainedclosestGene.csv\n")
cat("  DAAstroDAlostclosestGene.csv\n")
cat("  FIGURES/AstroDAVolcano.pdf\n")
cat("  FIGURES/AstroDAMAplot.pdf\n\n")


################################################################################
# 17. SAVE OBJECT
#
# combined_wnn now has:
#   - Annotation(combined_wnn) : GTF gene coordinates
#   - GC content columns       : from RegionStats (required by ChromVAR)
# Required by Scripts 07 (ChromVAR) and all astrocyte rewiring scripts.
################################################################################

saveRDS(
  combined_wnn,
  file = file.path(OUT_DIR, "06_combined_wnn_DA.rds")
)

cat("Saved: 06_combined_wnn_DA.rds\n")
cat("Next:  Run 07_ChromVAR_TF_Activity.R\n\n")

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT_DIR, "session_info_06.txt")
)
