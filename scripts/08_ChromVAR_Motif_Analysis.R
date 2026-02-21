################################################################################
# Script:  08_ChromVAR_Motif_Analysis.R
# Project: Heat_Call_Single_Nucleus_Multiome
# Author:  Prakrit Subba
# Date:    2026
#
# Paper:   Subba et al. 2026 — Multiome Profiling Reveals Astrocyte and
#          Neuroendocrine Targets of Prenatal Acoustic Programming in
#          Zebra Finch Embryos
#
# Description:
#   Transcription factor motif analysis using JASPAR 2020 + ChromVAR. Runs
#   both genome-wide (all cell types) and astrocyte-specific analyses.
#
#   Steps:
#     (1) Build JASPAR 2020 CORE vertebrate motif matrix
#     (2) Add motif annotations to combined_wnn → RunChromVAR (genome-wide)
#     (3) Differential TF activity: Heat Call vs Control across all cells
#         → FigureChromVARTopMotifsLogos.pdf
#         → TableChromVARDAHeatCallvsControl.csv
#     (4) Astrocyte-specific ChromVAR → differential TF activity
#         → AstroTFactivityDAchromVAR.csv
#     (5) Motif enrichment on astrocyte HC-gained DA peaks (FindMotifs)
#         → AstroDAgainedmotifenrichment.csv
#         → DAgainedmotifenrichment.pdf (bar + scatter)
#     (6) TF family summaries: HSF, NFI, SOX, Notch, bHLH
#
# Inputs:
#   - results/06_DA/06_combined_wnn_DA.rds   (has Annotation + RegionStats)
#   - results/06_DA/DAAstroDAHCgained.csv
#   - data/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna
#
# Outputs:
#   Object:
#     motif.combined_wnn.rds              (DOT notation — exact name)
#
#   Genome-wide ChromVAR:
#     TableChromVARDAHeatCallvsControl.csv (row.names = TRUE — exact)
#     FigureChromVARTopMotifsLogos.pdf
#
#   Astrocyte ChromVAR:
#     AstroTFactivityDAchromVAR.csv       (row.names = FALSE)
#
#   Motif enrichment (FindMotifs):
#     AstroDAgainedmotifenrichment.csv
#     DAgainedmotifenrichment.pdf
#
# NOTE on naming:
#   - motif.combined_wnn.rds              : DOT notation (not underscore)
#   - TableChromVARDAHeatCallvsControl.csv : NO underscores between words
#   - AstroTFactivityDAchromVAR.csv       : NO underscores between words
#   - AstroDAgainedmotifenrichment.csv    : NO underscores between words
#   - DAgainedmotifenrichment.pdf         : NO underscores between words
#
# NOTE on row.names:
#   - TableChromVARDAHeatCallvsControl.csv : row.names = TRUE  (motif IDs as rows)
#   - AstroTFactivityDAchromVAR.csv        : row.names = FALSE (motif_id column added)
#   - AstroDAgainedmotifenrichment.csv     : row.names = FALSE (motif_id column added)
#
# NOTE on motif matrix rownames:
#   rownames(motif.matrix) <- gsub("_", "-", rownames(motif.matrix))
#   UNDERSCORE → HYPHEN (not the reverse). Ensures peak names in the
#   motif matrix match Signac format (which uses hyphens).
#
# Previous Script: 07_Peak_Annotation.R
# Next Script:     09_Astrocyte_Rewiring.R
################################################################################


################################################################################
# 0. CONFIGURATION
################################################################################

WNN_RDS       <- "results/06_DA/06_combined_wnn_DA.rds"
DA_GAINED_CSV <- "results/06_DA/DAAstroDAHCgained.csv"
FASTA_PATH    <- "data/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna"

OUT_DIR    <- "results/08_ChromVAR"
DIR_FIG    <- file.path(OUT_DIR, "FIGURES")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG,  showWarnings = FALSE)

# Significance threshold
FDR_THRESHOLD <- 0.05

set.seed(1234)


################################################################################
# 1. LOAD LIBRARIES
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(chromVAR)
  library(JASPAR2020)
  library(TFBSTools)
  library(motifmatchr)
  library(Rsamtools)
  library(BSgenome)
  library(BiocParallel)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})


################################################################################
# 2. LOAD OBJECTS
################################################################################

cat("=== Loading objects ===\n")

combined_wnn <- readRDS(WNN_RDS)
DefaultAssay(combined_wnn) <- "ATAC"

cat("combined_wnn — cells:", ncol(combined_wnn),
    "| ATAC peaks:", nrow(combined_wnn@assays$ATAC), "\n")

# Validate Annotation exists (added in Script 06)
if (is.null(Annotation(combined_wnn))) {
  stop(
    "Annotation(combined_wnn) is NULL. ",
    "Run 06_Differential_Chromatin_Accessibility.R first."
  )
}

# Load astrocyte HC-gained DA peaks
dagained <- read.csv(DA_GAINED_CSV)
cat("Astrocyte HC-gained peaks:", nrow(dagained), "\n\n")

# Reference genome FASTA
fasta <- FaFile(FASTA_PATH)
indexFa(fasta)    # creates .fai index if absent


################################################################################
# 3. JASPAR 2020 MOTIF DATABASE
#
# Retrieves all CORE vertebrate position frequency matrices from JASPAR 2020.
# Parameters (exact from original code):
#   collection   = "CORE"
#   taxgroup     = "vertebrates"
#   all_versions = FALSE   (latest version of each motif only)
################################################################################

cat("=== Loading JASPAR 2020 CORE vertebrate motifs ===\n")

pfm <- getMatrixSet(
  x    = JASPAR2020,
  opts = list(
    collection   = "CORE",
    taxgroup     = "vertebrates",
    all_versions = FALSE
  )
)

cat("Motifs loaded:", length(pfm), "\n\n")


################################################################################
# 4. BUILD MOTIF MATRIX AND ADD TO COMBINED_WNN
#
# CreateMotifMatrix() scans each ATAC peak for each motif using the FASTA.
# CRITICAL rownames step:
#   rownames(motif.matrix) <- gsub("_", "-", rownames(motif.matrix))
#   Replaces UNDERSCORES with HYPHENS in peak rownames so they match
#   the Signac ATAC assay peak names (which use hyphens).
#
# RegionStats() recomputed here if needed (idempotent — safe to re-run).
################################################################################

cat("=== Building motif matrix for combined_wnn ===\n")
cat("  (CreateMotifMatrix — computationally intensive)\n\n")

# Ensure GC content stats are present (RegionStats from Script 06)
combined_wnn <- RegionStats(combined_wnn, genome = fasta)

motif.matrix <- CreateMotifMatrix(
  features = granges(combined_wnn@assays$ATAC),
  pwm      = pfm,
  genome   = fasta
)

# CRITICAL: replace underscores with hyphens in peak rownames
rownames(motif.matrix) <- gsub("_", "-", rownames(motif.matrix))

motif <- CreateMotifObject(data = motif.matrix, pwm = pfm)
Motifs(combined_wnn) <- motif

cat("Motif matrix dimensions:", dim(motif.matrix), "\n")
cat("Motifs added to combined_wnn.\n\n")


################################################################################
# 5. RUN CHROMVAR (GENOME-WIDE)
#
# Computes per-cell TF deviation z-scores across all cells and peaks.
# registerSerialParam() sets BiocParallel to single-core (avoids HPC fork
# issues). Adds a "chromvar" assay to the object.
#
# Object saved as: motif.combined_wnn.rds   (DOT notation — exact filename)
################################################################################

cat("=== Running ChromVAR (genome-wide) ===\n")
cat("  (Run on HPC — computationally intensive)\n\n")

registerSerialParam()

combined_wnn <- RunChromVAR(
  object = combined_wnn,
  genome = fasta,
  assay  = "ATAC"
)

# Save immediately with DOT notation filename (not underscore)
saveRDS(
  combined_wnn,
  file = file.path(OUT_DIR, "motif.combined_wnn.rds")
)

cat("ChromVAR complete.\n")
cat("Saved: motif.combined_wnn.rds\n\n")

# Convenience alias matching original code variable name
motif.combined_wnn <- combined_wnn


################################################################################
# 6. DIFFERENTIAL TF ACTIVITY — GENOME-WIDE (Heat Call vs Control)
#
# Parameters (exact from original code):
#   ident.1  = "Heat Call"
#   ident.2  = "Control"
#   only.pos = TRUE         (motifs with HIGHER activity in Heat Call only)
#   mean.fxn = rowMeans     (required for chromVAR deviation z-scores)
#   fc.name  = "avgdiff"    (Signac-recommended column name for chromVAR)
#
# Output: differential.activity
#   Column "avgdiff" (not "avg_log2FC") — chromVAR deviation difference.
#
# Saved as: TableChromVARDAHeatCallvsControl.csv
#   row.names = TRUE (motif IDs are preserved as rownames — exact original)
################################################################################

cat("=== Differential TF activity — genome-wide ===\n")

DefaultAssay(motif.combined_wnn) <- "chromvar"
Idents(motif.combined_wnn)       <- "Condition"

differential.activity <- FindMarkers(
  object   = motif.combined_wnn,
  ident.1  = "Heat Call",
  ident.2  = "Control",
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name  = "avgdiff"
)

cat("Significant TF motifs (padj <", FDR_THRESHOLD, "):",
    sum(differential.activity$p_val_adj < FDR_THRESHOLD, na.rm = TRUE), "\n")

# NOTE: row.names = TRUE — motif IDs are kept as rownames (exact original)
write.csv(
  differential.activity,
  file      = file.path(OUT_DIR, "TableChromVARDAHeatCallvsControl.csv"),
  row.names = TRUE
)

cat("Saved: TableChromVARDAHeatCallvsControl.csv  (row.names = TRUE)\n\n")

# Quick preview of top motifs
cat("Top 10 differentially active TF motifs:\n")
print(head(differential.activity, 10))
cat("\n")


################################################################################
# 7. MOTIF LOGO PLOT — TOP 12 TF MOTIFS (Figure 3)
#
# MotifPlot() requires the ATAC assay (not chromvar) for PFM retrieval.
# Plots PWM logos for the top 12 motifs ranked by avgdiff.
# Output: FigureChromVARTopMotifsLogos.pdf (12 × 8 inches; NO underscores)
################################################################################

cat("=== Generating motif logo plot (top 12) ===\n")

top_motifs_logos <- head(rownames(differential.activity), 12)

p_logos <- MotifPlot(
  object = motif.combined_wnn,
  motifs = top_motifs_logos,
  assay  = "ATAC"
) +
  patchwork::plot_annotation(
    title = "Top chromVAR Motifs: Heat Call vs Control"
  )

pdf(
  file.path(OUT_DIR, "FigureChromVARTopMotifsLogos.pdf"),
  width = 12, height = 8
)
print(p_logos)
dev.off()

cat("Saved: FigureChromVARTopMotifsLogos.pdf\n\n")


################################################################################
# 8. ASTROCYTE-SPECIFIC CHROMVAR
#
# Subsets to Astrocyte cluster, adds motifs + RegionStats if absent,
# then runs ChromVAR on the astrocyte subset.
#
# Differential TF activity parameters differ from genome-wide:
#   only.pos        = FALSE   (returns BOTH higher and lower activity)
#   mean.fxn        = rowMeans
#   fc.name         = "avgdiff"
#   logfc.threshold = 0       (return all tested motifs regardless of FC)
#
# Variable: difftfactivity  (not differential.activity)
# Column added: difftfactivity$motif_id <- rownames(difftfactivity)
# Saved as: AstroTFactivityDAchromVAR.csv  (row.names = FALSE)
################################################################################

cat("=== Astrocyte-specific ChromVAR ===\n")

# Subset to astrocytes
astroobj           <- subset(combined_wnn, subset = cluster_annotation == "Astrocyte")
DefaultAssay(astroobj) <- "ATAC"

cat("Astrocyte cells:", ncol(astroobj), "\n")
print(table(astroobj$Condition))
cat("\n")

# Add motifs and RegionStats to astrocyte subset if absent
if (is.null(Motifs(astroobj))) {
  cat("  Adding motifs to astrocyte object...\n")

  astroobj <- RegionStats(astroobj, genome = fasta)

  motif.matrix.astro <- CreateMotifMatrix(
    features = granges(astroobj@assays$ATAC),
    pwm      = pfm,
    genome   = fasta
  )

  # CRITICAL: underscore → hyphen (same as genome-wide step above)
  rownames(motif.matrix.astro) <- gsub("_", "-", rownames(motif.matrix.astro))

  motif.astro <- CreateMotifObject(data = motif.matrix.astro, pwm = pfm)
  Motifs(astroobj) <- motif.astro

  cat("  Motifs added to astrocyte object.\n\n")
} else {
  cat("  Motifs already present in astrocyte object.\n\n")
}

# Run ChromVAR on astrocytes
cat("  Running ChromVAR on astrocyte subset...\n")

registerSerialParam()

astroobj <- RunChromVAR(
  object = astroobj,
  genome = fasta,
  assay  = "ATAC"
)

cat("  ChromVAR complete for astrocytes.\n\n")

# Differential TF activity in astrocytes
DefaultAssay(astroobj) <- "chromvar"
Idents(astroobj)       <- "Condition"

difftfactivity         <- FindMarkers(
  object          = astroobj,
  ident.1         = "Heat Call",
  ident.2         = "Control",
  only.pos        = FALSE,      # Both directions (increased + decreased)
  mean.fxn        = rowMeans,
  fc.name         = "avgdiff",
  logfc.threshold = 0           # Return all motifs regardless of FC
)

difftfactivity$motif_id <- rownames(difftfactivity)

cat("Astrocyte TF motifs tested:", nrow(difftfactivity), "\n")
cat("Significant (padj <", FDR_THRESHOLD, "):",
    sum(difftfactivity$p_val_adj < FDR_THRESHOLD, na.rm = TRUE), "\n\n")

write.csv(
  difftfactivity,
  file      = file.path(OUT_DIR, "AstroTFactivityDAchromVAR.csv"),
  row.names = FALSE
)

cat("Saved: AstroTFactivityDAchromVAR.csv  (row.names = FALSE)\n\n")

# Reset to ATAC assay for downstream FindMotifs
DefaultAssay(astroobj) <- "ATAC"
Idents(astroobj)       <- "Condition"


################################################################################
# 9. FINDMOTIFS ON ASTROCYTE HC-GAINED DA PEAKS
#
# Identifies which TF motifs are significantly enriched in HC-gained DA peaks
# compared to the background set of all ATAC peaks in astrocytes.
#
# TF names annotated via: sapply(rownames(enriched_motifs), function(x) pfm[[x]]@name)
# Results sorted by pvalue ascending.
# Saved as: AstroDAgainedmotifenrichment.csv  (row.names = FALSE)
################################################################################

cat("=== FindMotifs on astrocyte HC-gained peaks ===\n")
cat("HC-gained peaks to test:", nrow(dagained), "\n")

if (nrow(dagained) > 0) {

  enriched_motifs <- FindMotifs(
    object   = astroobj,
    features = dagained$peak
  )

  enriched_motifs$motif_id <- rownames(enriched_motifs)

  # Annotate TF names from JASPAR PFM object
  enriched_motifs$tf_name <- sapply(
    rownames(enriched_motifs),
    function(x) {
      motif_info <- pfm[[x]]
      if (!is.null(motif_info)) motif_info@name else NA
    }
  )

  # Sort by p-value ascending
  enriched_motifs <- enriched_motifs %>% arrange(pvalue)

  write.csv(
    enriched_motifs,
    file      = file.path(OUT_DIR, "AstroDAgainedmotifenrichment.csv"),
    row.names = FALSE
  )

  cat("Motifs tested:", nrow(enriched_motifs), "\n")
  cat("Saved: AstroDAgainedmotifenrichment.csv\n\n")

  # Preview top 30
  cat("Top 30 enriched motifs:\n")
  print(
    enriched_motifs %>%
      select(tf_name, fold.enrichment, pvalue, percent.observed, percent.background) %>%
      head(30)
  )
  cat("\n")

} else {
  cat("No HC-gained peaks — skipping FindMotifs.\n\n")
}


################################################################################
# 10. TF FAMILY SUMMARIES
#
# Prints family-specific subsets of enriched_motifs for five biologically
# relevant TF families. Grep patterns (case-insensitive):
#   HSF/HSFA    : "HSF|HSFA|heat"     (Heat shock factors)
#   NFI         : "NFI"               (Nuclear Factor I — gliogenesis)
#   SOX         : "SOX"               (SOX family — neural)
#   Notch       : "RBPJ|HES|HEY"     (Notch signaling pathway)
#   bHLH        : "ASCL|NEUROG|HAND|TCF12"  (bHLH — stress/neural)
################################################################################

if (nrow(dagained) > 0 && exists("enriched_motifs")) {

  cat("=== TF Family Summaries ===\n\n")

  # Heat shock factors
  hsf_motifs <- enriched_motifs %>%
    filter(grepl("HSF|HSFA|heat", tf_name, ignore.case = TRUE))
  cat("HEAT SHOCK FACTORS (HSF/HSFA):", nrow(hsf_motifs), "\n")
  if (nrow(hsf_motifs) > 0)
    print(hsf_motifs %>% select(tf_name, fold.enrichment, pvalue))

  cat("\n")

  # Nuclear Factor I (gliogenesis)
  nfi_motifs <- enriched_motifs %>%
    filter(grepl("NFI", tf_name, ignore.case = TRUE))
  cat("NFI FAMILY — Gliogenesis:", nrow(nfi_motifs), "\n")
  if (nrow(nfi_motifs) > 0)
    print(nfi_motifs %>% select(tf_name, fold.enrichment, pvalue))

  cat("\n")

  # SOX family (neural)
  sox_motifs <- enriched_motifs %>%
    filter(grepl("SOX", tf_name, ignore.case = TRUE))
  cat("SOX FAMILY — Neural:", nrow(sox_motifs), "\n")
  if (nrow(sox_motifs) > 0)
    print(sox_motifs %>% select(tf_name, fold.enrichment, pvalue))

  cat("\n")

  # Notch pathway (RBPJ/HES/HEY)
  notch_motifs <- enriched_motifs %>%
    filter(grepl("RBPJ|HES|HEY", tf_name, ignore.case = TRUE))
  cat("NOTCH PATHWAY (RBPJ/HES/HEY):", nrow(notch_motifs), "\n")
  if (nrow(notch_motifs) > 0)
    print(notch_motifs %>% select(tf_name, fold.enrichment, pvalue))

  cat("\n")

  # bHLH family (stress/neural)
  bhlh_motifs <- enriched_motifs %>%
    filter(grepl("ASCL|NEUROG|HAND|TCF12", tf_name, ignore.case = TRUE))
  cat("bHLH FAMILY (ASCL/NEUROG/HAND/TCF12):", nrow(bhlh_motifs), "\n")
  if (nrow(bhlh_motifs) > 0)
    print(bhlh_motifs %>% select(tf_name, fold.enrichment, pvalue))

  cat("\n")
}


################################################################################
# 11. MOTIF ENRICHMENT FIGURES — DAgainedmotifenrichment.pdf
#
# Two panels saved to one PDF (12 × 8 inches):
#
# Panel 1 (p1): Horizontal bar chart of top 30 motifs
#   x    : TF name (reordered by -log10(pvalue))
#   y    : -log10(pvalue)
#   fill : fold.enrichment (gradient: lightblue → darkred)
#   Reference line: -log10(0.05)
#
# Panel 2 (p2): Scatter — fold enrichment vs significance
#   x     : fold.enrichment
#   y     : -log10(pvalue)
#   size  : percent.observed
#   color : -log10(pvalue) (gradient: blue → red)
#   Labels: top 20 TF names only (others NA — silently dropped by ggrepel)
#
# Label assignment (exact from original code):
#   enriched_motifs$label <- NA
#   enriched_motifs$label[1:min(20, nrow(enriched_motifs))] <-
#       enriched_motifs$tf_name[1:min(20, nrow(enriched_motifs))]
################################################################################

if (nrow(dagained) > 0 && exists("enriched_motifs")) {

  cat("=== Generating motif enrichment figures ===\n")

  pdf(
    file.path(DIR_FIG, "DAgainedmotifenrichment.pdf"),
    width = 12, height = 8
  )

  # Panel 1: top 30 bar chart
  top_motifs_bar <- head(enriched_motifs, 30)

  p1 <- ggplot(
    top_motifs_bar,
    aes(x = reorder(tf_name, -log10(pvalue)), y = -log10(pvalue),
        fill = fold.enrichment)
  ) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "lightblue", high = "darkred") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(
      title = paste0(
        "Top 30 Motifs Enriched in DA-Gained Peaks (n=", nrow(dagained), ")"
      ),
      x    = "Transcription Factor",
      y    = "-log10(p-value)",
      fill = "Fold Enrichment"
    ) +
    theme_minimal(base_size = 11)

  print(p1)

  # Panel 2: scatter — top 20 labeled, rest NA (silently dropped by ggrepel)
  enriched_motifs$label <- NA
  enriched_motifs$label[1:min(20, nrow(enriched_motifs))] <-
    enriched_motifs$tf_name[1:min(20, nrow(enriched_motifs))]

  p2 <- ggplot(
    enriched_motifs,
    aes(x = fold.enrichment, y = -log10(pvalue))
  ) +
    geom_point(
      aes(size = percent.observed, color = -log10(pvalue)),
      alpha = 0.6
    ) +
    geom_text_repel(
      aes(label = label),
      size         = 3,
      max.overlaps = 30
    ) +
    scale_color_gradient(low = "blue", high = "red") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(
      title = paste0(
        "Motif Enrichment in ", nrow(dagained), " HC-Gained Astrocyte Peaks"
      ),
      x    = "Fold Enrichment",
      y    = "-log10(p-value)",
      size = "Peaks with Motif"
    ) +
    theme_minimal()

  print(p2)

  dev.off()

  cat("Saved: DAgainedmotifenrichment.pdf\n\n")

} else {
  cat("No HC-gained peaks — skipping motif enrichment figures.\n\n")
}


################################################################################
# 12. SUMMARY
################################################################################

cat("========================================\n")
cat("SCRIPT 08 COMPLETE — CHROMVAR SUMMARY\n")
cat("========================================\n")
cat("Genome-wide ChromVAR:\n")
cat("  TF motifs tested              :", nrow(differential.activity), "\n")
cat("  Significant (only.pos=TRUE)   :",
    sum(differential.activity$p_val_adj < FDR_THRESHOLD, na.rm = TRUE), "\n")
cat("\nAstrocyte-specific ChromVAR:\n")
cat("  TF motifs tested              :", nrow(difftfactivity), "\n")
cat("  Significant (both directions) :",
    sum(difftfactivity$p_val_adj < FDR_THRESHOLD, na.rm = TRUE), "\n")
if (exists("enriched_motifs")) {
cat("\nFindMotifs (HC-gained peaks):\n")
cat("  HC-gained peaks tested        :", nrow(dagained), "\n")
cat("  Motifs enriched               :", nrow(enriched_motifs), "\n")
}
cat("\nOutput files:\n")
cat("  motif.combined_wnn.rds               (DOT notation)\n")
cat("  TableChromVARDAHeatCallvsControl.csv (row.names=TRUE)\n")
cat("  FigureChromVARTopMotifsLogos.pdf\n")
cat("  AstroTFactivityDAchromVAR.csv        (row.names=FALSE)\n")
cat("  AstroDAgainedmotifenrichment.csv\n")
cat("  FIGURES/DAgainedmotifenrichment.pdf\n\n")
cat("Next: Run 09_Astrocyte_Rewiring.R\n\n")


