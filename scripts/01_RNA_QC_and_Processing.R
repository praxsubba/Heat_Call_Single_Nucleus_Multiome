################################################################################
# Script:  01_RNA_QC_and_Processing.R
# Project: Heat_Call_Single_Nucleus_Multiome
# Author:  Prakrit Subba
# Date:    2026
#
# Paper:   Subba et al. 2026 — Multiome Profiling Reveals Astrocyte and
#          Neuroendocrine Targets of Prenatal Acoustic Programming in
#          Zebra Finch Embryos
#
# Description:
#   Loads 10X Chromium Single Cell Multiome (Gene Expression modality) output
#   for 8 embryonic zebra finch hypothalamus samples. Performs per-sample
#   Seurat object creation, QC filtering, merging, and SCTransform
#   normalization. Mitochondrial genes are identified from a reference GTF
#   and excluded from the analysis.
#
# Inputs:
#   - <BASE_DIR>/<SAMPLE_ID>_combined_analysis/outs/filtered_feature_bc_matrix.h5
#   - MT_GTF_PATH : GTF file listing mitochondrial chromosome genes
#
# Outputs:
#   - merged_seurat_filt        : QC-filtered, SCTransform-normalized Seurat object
#   - QC_VlnPlot_prefilter.pdf  : Pre-filter QC violin plots
#   - QC_VlnPlot_postfilter.pdf : Post-filter QC violin plots
#   - 01_merged_seurat_filt.rds : Saved Seurat object for downstream scripts
#
# Manuscript Reference:
#   Methods — Single-Nucleus RNA-seq Processing and QC
#   Figure 1B (RNA UMAP, carried forward to Script 03)
#
# Next Script: 02_ATAC_QC_and_Processing.R
################################################################################


################################################################################
# 0. CONFIGURATION — Update these paths before running
################################################################################

# Root directory containing per-sample CellRanger ARC output folders
# Expected structure: <BASE_DIR>/S1_combined_analysis/outs/
BASE_DIR <- "/path/to/cellranger_arc_outputs"

# GTF file for mitochondrial chromosome genes (species-specific)
# Used to identify and remove mt genes before downstream analysis
MT_GTF_PATH <- "/path/to/reference/mt_chromosome.gtf"

# Output directory for QC plots and saved objects
OUT_DIR <- "results/01_RNA_QC"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Sample IDs matching CellRanger output folder prefixes
SAMPLE_IDS <- paste0("S", 1:8)

# QC thresholds (as reported in Methods)
MIN_FEATURES  <- 200
MAX_FEATURES  <- 3500
MIN_COUNTS    <- 500
MAX_COUNTS    <- 25000
MAX_MT_PCT    <- 5


################################################################################
# 1. LOAD LIBRARIES
################################################################################

library(Seurat)
library(SeuratObject)
library(stringr)
library(dplyr)
library(ggplot2)
library(patchwork)


################################################################################
# 2. LOAD MITOCHONDRIAL GENE LIST
#
# Mitochondrial genes are identified from a species-specific GTF and removed
# before downstream analysis to avoid confounding from mt-gene expression,
# which can indicate low-quality nuclei or cytoplasmic contamination in
# snRNA-seq.
################################################################################

cat("=== Loading mitochondrial gene list ===\n")

mt_genes     <- read.delim(MT_GTF_PATH, header = FALSE, comment.char = "#")
mt_gene_ids  <- unique(
  str_extract(mt_genes$V9, "(?<=gene_id\\s)[^\";]+")
) %>%
  na.omit() %>%
  as.character()

cat("Mitochondrial genes identified:", length(mt_gene_ids), "\n\n")


################################################################################
# 3. CREATE PER-SAMPLE SEURAT OBJECTS
#
# Reads filtered_feature_bc_matrix.h5 from each sample's CellRanger ARC
# output directory. For multiome outputs, extracts the "Gene Expression"
# modality. Mitochondrial genes are removed from the feature set, and
# percent.mt is calculated for use as an additional QC metric.
#
# Filters at object creation:
#   min.cells = 3    : gene must be detected in ≥3 cells to be retained
#   min.features = 200 : cell must have ≥200 detected genes
################################################################################

cat("=== Creating per-sample Seurat objects ===\n")

seurat_objects <- list()

for (id in SAMPLE_IDS) {

  data_dir <- file.path(BASE_DIR, paste0(id, "_combined_analysis"), "outs")
  h5_path  <- file.path(data_dir, "filtered_feature_bc_matrix.h5")

  # Read 10X multiome h5; extract Gene Expression slot if multiome format
  data_raw <- Read10X_h5(h5_path)
  if (is.list(data_raw) && "Gene Expression" %in% names(data_raw)) {
    data <- data_raw$`Gene Expression`
  } else {
    data <- data_raw
  }

  # Create Seurat object with initial gene/cell filters
  seurat_obj <- CreateSeuratObject(
    counts      = data,
    project     = id,
    min.cells   = 3,
    min.features = 200
  )

  # Remove mitochondrial genes from the feature set
  # (GTF uses underscores; Seurat rownames use hyphens)
  mt_gene_ids_seurat <- mt_gene_ids %>% str_replace_all("_", "-")
  seurat_obj <- subset(
    seurat_obj,
    features = setdiff(rownames(seurat_obj), mt_gene_ids_seurat)
  )

  # Calculate percent mitochondrial reads as QC metric
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  seurat_objects[[id]] <- seurat_obj
  cat("  Sample", id, "— Cells:", ncol(seurat_obj),
      "| Genes:", nrow(seurat_obj), "\n")
}

cat("\n")


################################################################################
# 4. ADD EXPERIMENTAL METADATA
#
# Assigns Condition (Heat Call / Control) and Sex (Male / Female) to each
# sample based on the experimental design:
#   S1: Heat Call, Male    |  S2: Heat Call, Female
#   S3: Control,   Male    |  S4: Control,   Female
#   S5: Heat Call, Male    |  S6: Heat Call, Female
#   S7: Control,   Male    |  S8: Control,   Female
################################################################################

cat("=== Adding experimental metadata ===\n")

metadata <- data.frame(
  Sample    = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"),
  Condition = c("Heat Call", "Heat Call", "Control",   "Control",
                "Heat Call", "Heat Call", "Control",   "Control"),
  Sex       = c("Male",  "Female", "Male",  "Female",
                "Male",  "Female", "Male",  "Female"),
  stringsAsFactors = FALSE
)

for (i in seq_along(seurat_objects)) {
  seurat_objects[[i]]$Condition <- metadata$Condition[i]
  seurat_objects[[i]]$Sex       <- metadata$Sex[i]
}

cat("Metadata assigned. Sample breakdown:\n")
print(metadata)
cat("\n")


################################################################################
# 5. MERGE ALL SAMPLES INTO A SINGLE SEURAT OBJECT
#
# Cell barcodes are prefixed with sample IDs (e.g., "S1_ATCG...") to avoid
# barcode collisions across samples.
################################################################################

cat("=== Merging all samples ===\n")

merged_seurat <- merge(
  seurat_objects[[1]],
  y           = seurat_objects[-1],
  add.cell.ids = SAMPLE_IDS,
  project     = "HeatCall_Multiome"
)

cat("Merged object — Total cells:", ncol(merged_seurat),
    "| Total genes:", nrow(merged_seurat), "\n\n")


################################################################################
# 6. QC VISUALIZATION (PRE-FILTER)
#
# Violin plots of nFeature_RNA, nCount_RNA, and percent.mt across all cells
# before applying thresholds. Saved for reporting and QC documentation.
################################################################################

cat("=== Generating pre-filter QC plots ===\n")

p_prefilter <- VlnPlot(
  merged_seurat,
  features  = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol      = 3,
  pt.size   = 0,
  split.by  = "Condition"
) &
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
    plot.title   = element_text(size = 12, face = "bold")
  )

ggsave(
  filename = file.path(OUT_DIR, "QC_VlnPlot_prefilter.pdf"),
  plot     = p_prefilter,
  width    = 14, height = 6
)

# Feature scatter: nCount vs nFeature and nCount vs percent.mt
p_scatter1 <- FeatureScatter(merged_seurat, "nCount_RNA", "nFeature_RNA")
p_scatter2 <- FeatureScatter(merged_seurat, "nCount_RNA", "percent.mt")

ggsave(
  filename = file.path(OUT_DIR, "QC_FeatureScatter_prefilter.pdf"),
  plot     = p_scatter1 + p_scatter2,
  width    = 12, height = 5
)

cat("Pre-filter QC plots saved to:", OUT_DIR, "\n\n")


################################################################################
# 7. QC FILTERING
#
# Applies thresholds to remove:
#   - Low-quality nuclei (too few genes/counts — likely empty droplets)
#   - Potential multiplets (too many genes/counts)
#   - Mitochondrially contaminated nuclei (percent.mt > 5%)
#
# Final thresholds (as reported in Methods):
#   nFeature_RNA : 200 – 3,500
#   nCount_RNA   : 500 – 25,000
#   percent.mt   : < 5%
################################################################################

cat("=== Applying QC filters ===\n")
cat("  Pre-filter cells :", ncol(merged_seurat), "\n")

merged_seurat_filt <- subset(
  merged_seurat,
  subset = nFeature_RNA > MIN_FEATURES &
           nFeature_RNA < MAX_FEATURES &
           nCount_RNA   > MIN_COUNTS   &
           nCount_RNA   < MAX_COUNTS   &
           percent.mt   < MAX_PCT_MT
)

cat("  Post-filter cells:", ncol(merged_seurat_filt), "\n")
cat("  Cells removed    :", ncol(merged_seurat) - ncol(merged_seurat_filt), "\n\n")

# Per-sample cell counts after filtering
cat("Post-filter cell counts per sample:\n")
print(table(merged_seurat_filt$orig.ident))
cat("\n")

# Per-condition cell counts
cat("Post-filter cell counts per condition:\n")
print(table(merged_seurat_filt$Condition))
cat("\n")


################################################################################
# 8. POST-FILTER QC VISUALIZATION
################################################################################

cat("=== Generating post-filter QC plots ===\n")

p_postfilter <- VlnPlot(
  merged_seurat_filt,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol     = 3,
  pt.size  = 0,
  split.by = "Condition"
) &
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title  = element_text(size = 12, face = "bold")
  )

ggsave(
  filename = file.path(OUT_DIR, "QC_VlnPlot_postfilter.pdf"),
  plot     = p_postfilter,
  width    = 14, height = 6
)

cat("Post-filter QC plots saved.\n\n")


################################################################################
# 9. NORMALIZATION — SCTransform
#
# SCTransform is applied to the merged filtered object. This regularized
# negative binomial regression removes the influence of sequencing depth on
# gene expression, producing Pearson residuals stored in the "SCT" assay.
# The resulting object is used as input for PCA and WNN clustering in
# Script 03_WNN_Integration_and_Clustering.R.
################################################################################

cat("=== Running SCTransform normalization ===\n")
cat("(This may take several minutes depending on cell count)\n\n")

merged_seurat_filt <- SCTransform(
  merged_seurat_filt,
  verbose = TRUE
)

cat("\nSCTransform complete.\n")
cat("Default assay:", DefaultAssay(merged_seurat_filt), "\n")
cat("Assays in object:", paste(names(merged_seurat_filt@assays), collapse = ", "), "\n\n")


################################################################################
# 10. SUMMARY STATISTICS
################################################################################

cat("========================================\n")
cat("SCRIPT 01 COMPLETE — QC SUMMARY\n")
cat("========================================\n")
cat("Samples processed       :", length(SAMPLE_IDS), "\n")
cat("Pre-filter nuclei       :", ncol(merged_seurat), "\n")
cat("Post-filter nuclei      :", ncol(merged_seurat_filt), "\n")
cat("Nuclei retained         :",
    round(ncol(merged_seurat_filt) / ncol(merged_seurat) * 100, 1), "%\n")
cat("Genes in SCT assay      :", nrow(merged_seurat_filt[["SCT"]]), "\n")
cat("QC thresholds applied:\n")
cat("  nFeature_RNA :", MIN_FEATURES, "–", MAX_FEATURES, "\n")
cat("  nCount_RNA   :", MIN_COUNTS,   "–", MAX_COUNTS,   "\n")
cat("  percent.mt   : <", MAX_MT_PCT, "%\n")
cat("\nCondition breakdown:\n")
print(table(merged_seurat_filt$Condition, merged_seurat_filt$Sex))
cat("\n")


################################################################################
# 11. SAVE OUTPUT
#
# Saves the filtered, SCTransform-normalized Seurat object for use in:
#   02_ATAC_QC_and_Processing.R  (ATAC side processed in parallel)
#   03_WNN_Integration_and_Clustering.R (requires this object)
################################################################################

cat("=== Saving output object ===\n")

saveRDS(
  merged_seurat_filt,
  file = file.path(OUT_DIR, "01_merged_seurat_filt.rds")
)

cat("Saved: 01_merged_seurat_filt.rds\n")
cat("Next:  Run 02_ATAC_QC_and_Processing.R\n\n")



