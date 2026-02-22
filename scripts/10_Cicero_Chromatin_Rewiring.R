################################################################################
# Script:  10_Cicero_Chromatin_Rewiring.R
# Project: Heat_Call_Single_Nucleus_Multiome
# Author:  Prakrit Subba
# Date:    2026
#
# Produces:
#   Manuscript Figure 2 (gene integration), Figure 5 (locus plots via
#   seuratwithlinks.rds), Supplemental Tables (rewired connections,
#   promoter-anchored gene ranking, DA-Cicero overlap, LinkPeaks).
#
# Key thresholds (manuscript Methods):
#   ACCTHRESHOLD    = 0.01   (peak accessible if >1% cells open)
#   ACTIVETHRESHOLD = 0.05   (link active if both endpoints >5% accessible)
#   COACCESSCUTOFF  = 0.15   (Cicero co-accessibility cutoff)
#
# Inputs:
#   - results/06_DA/06_combined_wnn_DA.rds
#   - <CICERO_DIR>/baselineciceroconnections.rds
#   - <CICERO_DIR>/baselineciceroccans.rds
#   - results/06_DA/DAAstroDAHCgained.csv
#   - results/06_DA/DAAstroDAgainedclosestGene.csv
#   - results/05_DEG/FilteredAllClustersHeatCallvsControlDEG.csv
#   - data/GCF_003957565.2_bTaeGut1.4.pri_genomic.gtf
#
# Outputs (exact names used by downstream / manuscript):
#   rewiredhc.rds / .csv
#   rewiredctrl.rds / .csv
#   connsclusterannotated.rds / .csv
#   ccanscluster.rds / .csv
#   peakaccessibility.rds
#   connsclustersignacformat.rds
#   ccansclustersignacformat.rds
#   seuratwithlinks.rds               ← Figure 5 CoveragePlot source
#   HCpromotergenerank.csv
#   HCpromotergenerankwithRNA.csv
#   rewiredhcpromoteranchored.csv
#   CicerorewiringDAoverlap.csv
#   peaksHCrewiredANDDAgained.txt
#   Overlap57peaksDAstats.csv
#   Overlap57peaksRNAlinks.csv
#   Overlap57peaksFULLannotation.csv
#   Overlap57peaksgenesummary.csv
#   AllDAgainedLinkPeakslinks.csv
#   AstroDAgeneswithstats.csv
#   AstroDATop50genesmanuscript.csv
#   FIGURES/Fig2GeneIntegration.pdf   ← Figure 2
#
# Previous Script: 09_hdWGCNA_Coexpression_Networks.R
# Next Script:     11_Figure5_Locus_Plots.R
################################################################################


################################################################################
# 0. CONFIGURATION
################################################################################

WNN_RDS       <- "results/06_DA/06_combined_wnn_DA.rds"

CICERO_DIR    <- "data2/george_lab/psubba/cicero_runs/zfinch_cicero_2025-12-13_15:44:35"
CONNS_RDS     <- file.path(CICERO_DIR, "baselineciceroconnections.rds")
CCANS_RDS     <- file.path(CICERO_DIR, "baselineciceroccans.rds")

GTF_PATH      <- "data/GCF_003957565.2_bTaeGut1.4.pri_genomic.gtf"

DA_GAINED_CSV       <- "results/06_DA/DAAstroDAHCgained.csv"
DA_LOST_CSV         <- "results/06_DA/DAAstroDAHClost.csv"
CLOSEST_GAINED_CSV  <- "results/07_Peak_Annotation/DAAstroDAgainedclosestGene.csv"
DEG_CSV             <- "results/05_DEG/FilteredAllClustersHeatCallvsControlDEG.csv"

OUT_DIR  <- "results/10_Cicero"
DIR_FIG  <- file.path(OUT_DIR, "FIGURES")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG, showWarnings = FALSE)

ACCTHRESHOLD    <- 0.01
ACTIVETHRESHOLD <- 0.05
COACCESSCUTOFF  <- 0.15
CLUSTERID       <- "Astrocyte"

set.seed(1234)


################################################################################
# 1. LIBRARIES
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(GenomicRanges)
  library(rtracklayer)
  library(Matrix)
})


################################################################################
# 2. LOAD OBJECTS
################################################################################

combined_wnn <- readRDS(WNN_RDS)
DefaultAssay(combined_wnn) <- "ATAC"

connsbaseline <- readRDS(CONNS_RDS)
ccansbaseline <- readRDS(CCANS_RDS)

# Ensure character (not factor) peak names
for (col in c("Peak1", "Peak2")) {
  if (is.factor(connsbaseline[[col]]))
    connsbaseline[[col]] <- as.character(connsbaseline[[col]])
}
if (is.factor(ccansbaseline$Peak))
  ccansbaseline$Peak <- as.character(ccansbaseline$Peak)

dagained       <- read.csv(DA_GAINED_CSV)
dalost         <- read.csv(DA_LOST_CSV)
closest_gained <- read.csv(CLOSEST_GAINED_CSV)


################################################################################
# 3. SUBSET ASTROCYTES AND COMPUTE PEAK ACCESSIBILITY
#
# Per-peak fraction of cells with non-zero counts in each condition.
# Used to define accessible (>ACCTHRESHOLD) and active (>ACTIVETHRESHOLD) links.
################################################################################

stopifnot(
  "cluster_annotation" %in% colnames(combined_wnn@meta.data),
  "Condition"          %in% colnames(combined_wnn@meta.data)
)

cells_clust <- colnames(combined_wnn)[combined_wnn$cluster_annotation == CLUSTERID]
cells_hc    <- cells_clust[combined_wnn$Condition[cells_clust] == "Heat Call"]
cells_ctrl  <- cells_clust[combined_wnn$Condition[cells_clust] == "Control"]

atac_counts <- GetAssayData(combined_wnn, assay = "ATAC", layer = "counts")
atac_clust  <- atac_counts[, cells_clust]

peaks_hc   <- Matrix::rowMeans(atac_clust[, cells_hc]   > 0)
peaks_ctrl <- Matrix::rowMeans(atac_clust[, cells_ctrl] > 0)

peak_accessibility <- data.frame(
  peak     = rownames(atac_clust),
  acchc    = peaks_hc,
  accctrl  = peaks_ctrl,
  deltaacc = peaks_hc - peaks_ctrl,
  stringsAsFactors = FALSE
)

saveRDS(peak_accessibility,
        file.path(OUT_DIR, "peakaccessibility.rds"))


################################################################################
# 4. FILTER TO ASTROCYTE-ACCESSIBLE CONNECTIONS
################################################################################

accessible_peaks <- peak_accessibility$peak[
  peak_accessibility$acchc   > ACCTHRESHOLD |
  peak_accessibility$accctrl > ACCTHRESHOLD
]

connsclust <- connsbaseline %>%
  filter(Peak1 %in% accessible_peaks, Peak2 %in% accessible_peaks)


################################################################################
# 5. ANNOTATE CONNECTIONS WITH CONDITION-SPECIFIC ACCESSIBILITY
################################################################################

connsclust <- connsclust %>%
  left_join(peak_accessibility %>% select(peak, acchc, accctrl, deltaacc),
            by = c("Peak1" = "peak")) %>%
  rename(Peak1acchc = acchc, Peak1accctrl = accctrl, Peak1delta = deltaacc) %>%
  left_join(peak_accessibility %>% select(peak, acchc, accctrl, deltaacc),
            by = c("Peak2" = "peak")) %>%
  rename(Peak2acchc = acchc, Peak2accctrl = accctrl, Peak2delta = deltaacc)

connsclust$meanacchc   <- (connsclust$Peak1acchc   + connsclust$Peak2acchc)   / 2
connsclust$meanaccctrl <- (connsclust$Peak1accctrl + connsclust$Peak2accctrl) / 2
connsclust$deltameanacc <- connsclust$meanacchc - connsclust$meanaccctrl

connsclust$activehc   <- connsclust$Peak1acchc   > ACTIVETHRESHOLD &
                         connsclust$Peak2acchc   > ACTIVETHRESHOLD
connsclust$activectrl <- connsclust$Peak1accctrl > ACTIVETHRESHOLD &
                         connsclust$Peak2accctrl > ACTIVETHRESHOLD


################################################################################
# 6. DEFINE HC-SPECIFIC AND CTRL-SPECIFIC REWIRED CONNECTIONS
################################################################################

rewiredhc   <- connsclust %>%
  filter(activehc == TRUE, activectrl == FALSE, coaccess > COACCESSCUTOFF)

rewiredctrl <- connsclust %>%
  filter(activectrl == TRUE, activehc == FALSE, coaccess > COACCESSCUTOFF)

ccansclust <- ccansbaseline %>%
  filter(Peak %in% accessible_peaks) %>%
  left_join(peak_accessibility %>% select(peak, deltaacc, acchc, accctrl),
            by = c("Peak" = "peak"))


################################################################################
# 7. SAVE OUTPUTS (ORIGINAL 4-PART PEAK FORMAT)
################################################################################

saveRDS(connsclust,        file.path(OUT_DIR, "connsclusterannotated.rds"))
saveRDS(rewiredhc,         file.path(OUT_DIR, "rewiredhc.rds"))
saveRDS(rewiredctrl,       file.path(OUT_DIR, "rewiredctrl.rds"))
saveRDS(ccansclust,        file.path(OUT_DIR, "ccanscluster.rds"))

write.csv(connsclust,  file.path(OUT_DIR, "connsclusterannotated.csv"),  row.names = FALSE)
write.csv(rewiredhc,   file.path(OUT_DIR, "rewiredhc.csv"),              row.names = FALSE)
write.csv(rewiredctrl, file.path(OUT_DIR, "rewiredctrl.csv"),            row.names = FALSE)
write.csv(ccansclust,  file.path(OUT_DIR, "ccanscluster.csv"),           row.names = FALSE)


################################################################################
# 8. CONVERT PEAK FORMAT FOR SIGNAC
#
# Cicero uses 4-part peak names for NCBI chromosomes: NC-044211.2-start-end
# Signac requires 3-part: NC044211.2-start-end (hyphen → underscore in chr)
################################################################################

convert_peak_format <- function(peaks) {
  parts_list <- strsplit(peaks, "-")
  sapply(parts_list, function(p) {
    if      (length(p) == 4) paste0(p[1], "_", p[2], "-", p[3], "-", p[4])
    else if (length(p) == 3) paste(p, collapse = "-")
    else                     paste(p, collapse = "-")
  }, USE.NAMES = FALSE)
}

conns_signac       <- connsclust
conns_signac$Peak1 <- convert_peak_format(connsclust$Peak1)
conns_signac$Peak2 <- convert_peak_format(connsclust$Peak2)

ccans_signac       <- ccansclust
ccans_signac$Peak  <- convert_peak_format(ccansclust$Peak)

# Verify 3-part conversion
test_split <- strsplit(conns_signac$Peak1[1], "-")[[1]]
stopifnot("Peak conversion to 3-part failed" = length(test_split) == 3)

saveRDS(conns_signac, file.path(OUT_DIR, "connsclustersignacformat.rds"))
saveRDS(ccans_signac, file.path(OUT_DIR, "ccansclustersignacformat.rds"))


################################################################################
# 9. ADD LINKS TO SEURAT → seuratwithlinks.rds
#
# Required for CoveragePlot(links = TRUE) in Figure 5 locus plots.
################################################################################

# Update ATAC assay rownames to Signac 3-part format
old_peaks <- rownames(combined_wnn@assays$ATAC)
new_peaks <- convert_peak_format(old_peaks)

atac_assay <- combined_wnn@assays$ATAC
rownames(atac_assay@counts) <- new_peaks
if (!is.null(atac_assay@data))          rownames(atac_assay@data)          <- new_peaks
if (!is.null(atac_assay@meta.features)) rownames(atac_assay@meta.features) <- new_peaks
combined_wnn@assays$ATAC <- atac_assay

links <- SignacConnectionsToLinks(
  conns = conns_signac,
  ccans = ccans_signac
)
Links(combined_wnn@assays$ATAC) <- links

saveRDS(combined_wnn, file.path(OUT_DIR, "seuratwithlinks.rds"))


################################################################################
# 10. PROMOTER ANCHORING + RNA INTEGRATION
#
# Identifies HC-rewired connections where at least one endpoint overlaps a
# gene promoter (±2 kb). Ranks genes by number of HC-specific promoter loops
# and joins with RNA DEG results.
#
# Outputs:
#   rewiredhcpromoteranchored.csv      (promoter-anchored HC links)
#   HCpromotergenerank.csv             (gene ranking by n HC loops)
#   HCpromotergenerankwithRNA.csv      (joined with RNA DE)
################################################################################

gtf   <- import(GTF_PATH)
genes <- gtf[gtf$type == "gene"]

# Use best available gene name column
name_col <- if ("gene_name"  %in% names(mcols(genes))) "gene_name"  else
            if ("gene"       %in% names(mcols(genes))) "gene"       else
            if ("Name"       %in% names(mcols(genes))) "Name"       else "gene_id"
mcols(genes)$gene_name <- mcols(genes)[[name_col]]

# Build promoter windows (2 kb upstream, 2 kb downstream)
# Fix chromosome naming: NC_044211.2 → NC-044211.2 to match Cicero peaks
promoters_gr              <- GenomicRanges::promoters(genes, upstream = 2000, downstream = 2000)
seqlevels(genes)          <- gsub("_", "-", seqlevels(genes))
seqlevels(promoters_gr)   <- seqlevels(genes)

# Helper: peak string → GRanges (handles 4-part NCBI chr names)
peak_to_gr <- function(peaks_chr_start_end) {
  parts_list <- strsplit(peaks_chr_start_end, "-")
  chr   <- character(length(parts_list))
  start <- integer(length(parts_list))
  end   <- integer(length(parts_list))
  for (i in seq_along(parts_list)) {
    p <- parts_list[[i]]
    if (length(p) == 4) {
      chr[i]   <- paste0(p[1], "-", p[2])
      start[i] <- as.integer(p[3])
      end[i]   <- as.integer(p[4])
    } else if (length(p) == 3) {
      chr[i]   <- p[1]
      start[i] <- as.integer(p[2])
      end[i]   <- as.integer(p[3])
    }
  }
  GRanges(seqnames = chr, ranges = IRanges(start = start, end = end),
          peak_id  = peaks_chr_start_end)
}

# Convert rewiredhc Peak1 and Peak2 to GRanges and find promoter overlaps
p1_gr  <- peak_to_gr(rewiredhc$Peak1)
p2_gr  <- peak_to_gr(rewiredhc$Peak2)

rewiredhc$Peak1gene <- NA_character_
rewiredhc$Peak2gene <- NA_character_

p1_hit <- findOverlaps(p1_gr, promoters_gr, ignore.strand = TRUE)
p2_hit <- findOverlaps(p2_gr, promoters_gr, ignore.strand = TRUE)

rewiredhc$Peak1gene[queryHits(p1_hit)] <-
  mcols(promoters_gr)$gene_name[subjectHits(p1_hit)]
rewiredhc$Peak2gene[queryHits(p2_hit)] <-
  mcols(promoters_gr)$gene_name[subjectHits(p2_hit)]

# Keep links where at least one endpoint overlaps a promoter
rewiredhc_prom <- rewiredhc %>%
  filter(!is.na(Peak1gene) | !is.na(Peak2gene)) %>%
  mutate(promoter_gene = dplyr::coalesce(Peak1gene, Peak2gene))

write.csv(rewiredhc_prom,
          file.path(OUT_DIR, "rewiredhcpromoteranchored.csv"),
          row.names = FALSE)

# Rank genes by number of HC-specific promoter loops
hc_gene_rank <- rewiredhc_prom %>%
  group_by(promoter_gene) %>%
  summarise(
    nhclinks     = n(),
    meancoaccess = mean(coaccess,    na.rm = TRUE),
    meanPeak1acchc = mean(Peak1acchc, na.rm = TRUE),
    meanPeak2acchc = mean(Peak2acchc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(nhclinks), desc(meancoaccess))

write.csv(hc_gene_rank,
          file.path(OUT_DIR, "HCpromotergenerank.csv"),
          row.names = FALSE)

# Join with RNA DEG results
deg_all <- read.csv(DEG_CSV)

rna_da <- FindMarkers(
  object          = subset(combined_wnn, subset = cluster_annotation == "Astrocyte"),
  ident.1         = "Heat Call",
  ident.2         = "Control",
  assay           = "SCT",
  logfc.threshold = 0,
  min.pct         = 0.1
)
rna_da$gene <- rownames(rna_da)

hc_gene_rank_rna <- hc_gene_rank %>%
  left_join(rna_da, by = c("promoter_gene" = "gene"))

write.csv(hc_gene_rank_rna,
          file.path(OUT_DIR, "HCpromotergenerankwithRNA.csv"),
          row.names = FALSE)


################################################################################
# 11. CICERO–DA OVERLAP
#
# Counts overlap between HC-rewired peaks and DA-gained/lost peaks.
# Saves overlap table and the specific dual-hit peak list.
# Outputs:
#   CicerorewiringDAoverlap.csv
#   peaksHCrewiredANDDAgained.txt
################################################################################

rewiredhc_peaks  <- unique(c(rewiredhc$Peak1,    rewiredhc$Peak2))
rewiredct_peaks  <- unique(c(rewiredctrl$Peak1,  rewiredctrl$Peak2))

overlap_tbl <- tibble(
  comparison = c("HCrewired-DAgained", "HCrewired-DAlost",
                 "Ctrlrewired-DAgained", "Ctrlrewired-DAlost"),
  nrewired   = c(length(rewiredhc_peaks),  length(rewiredhc_peaks),
                 length(rewiredct_peaks),  length(rewiredct_peaks)),
  nda        = c(nrow(dagained), nrow(dalost), nrow(dagained), nrow(dalost)),
  noverlap   = c(sum(rewiredhc_peaks  %in% dagained$peak),
                 sum(rewiredhc_peaks  %in% dalost$peak),
                 sum(rewiredct_peaks  %in% dagained$peak),
                 sum(rewiredct_peaks  %in% dalost$peak))
) %>%
  mutate(
    pctrewired = round(100 * noverlap / nrewired, 1),
    pctda      = round(100 * noverlap / pmax(nda, 1), 1)
  )

write.csv(overlap_tbl,
          file.path(OUT_DIR, "CicerorewiringDAoverlap.csv"),
          row.names = FALSE)

# Dual-hit peaks: HC-rewired AND DA-gained
hcboth <- intersect(rewiredhc_peaks, dagained$peak)
writeLines(hcboth, file.path(OUT_DIR, "peaksHCrewiredANDDAgained.txt"))


################################################################################
# 12. LINKPEAKS ON DUAL-HIT PEAKS (HC-REWIRED ∩ DA-GAINED)
#
# Computes peak-gene expression correlations for the dual-hit peaks using
# Signac's LinkPeaks (Pearson, distance = 500 kb, min.cells = 10).
# Genes searched are those within 500 kb of dual-hit peaks.
#
# Outputs:
#   Overlap57peaksDAstats.csv
#   Overlap57peaksRNAlinks.csv
#   Overlap57peaksFULLannotation.csv
#   Overlap57peaksgenesummary.csv
#   AllDAgainedLinkPeakslinks.csv
#   AstroDAgeneswithstats.csv
#   AstroDATop50genesmanuscript.csv
################################################################################

if (length(hcboth) > 0) {

  overlap_da <- dagained %>% filter(peak %in% hcboth) %>% arrange(p_val_adj)
  write.csv(overlap_da, file.path(OUT_DIR, "Overlap57peaksDAstats.csv"),
            row.names = FALSE)

  # Subset astrocyte object with SCT assay for LinkPeaks
  astrolink <- subset(combined_wnn, subset = cluster_annotation == "Astrocyte")
  DefaultAssay(astrolink) <- "ATAC"
  stopifnot("SCT" %in% names(astrolink@assays))

  # Find genes within 500 kb of dual-hit peaks
  peak_ranges          <- peak_to_gr(hcboth)
  peak_ranges_extended <- peak_ranges + 500000
  annotation           <- Annotation(astrolink)

  if (is.null(annotation)) stop("No Annotation found in astrolink — run Script 06 first.")

  overlaps    <- findOverlaps(peak_ranges_extended, annotation)
  nearby_genes <- unique(annotation$gene_name[subjectHits(overlaps)])
  nearby_genes <- nearby_genes[!is.na(nearby_genes)]

  # LinkPeaks: dual-hit peaks only
  astrolink <- LinkPeaks(
    object          = astrolink,
    peak.assay      = "ATAC",
    expression.assay = "SCT",
    peak.slot       = "counts",
    expression.slot = "data",
    genes.use       = nearby_genes,
    distance        = 500000,
    min.cells       = 10,
    method          = "pearson"
  )

  links         <- Links(astrolink)
  links_df      <- data.frame(gene  = links$gene,
                               score = links$score,
                               pvalue = links$pvalue,
                               peak  = links$peak)
  overlap_links_clean <- links_df %>%
    filter(peak %in% hcboth) %>%
    arrange(desc(abs(score))) %>%
    select(peak, gene, score, pvalue)

  overlap_full <- overlap_links_clean %>%
    left_join(overlap_da %>% select(peak, avg_log2FC, p_val_adj, pct.1, pct.2),
              by = "peak")

  gene_summary <- overlap_links_clean %>%
    group_by(gene) %>%
    summarise(npeaks      = n(),
              maxabscorr  = max(abs(score)),
              avgcorr     = mean(score),
              minpvalue   = min(pvalue),
              .groups     = "drop") %>%
    arrange(desc(maxabscorr))

  write.csv(overlap_links_clean, file.path(OUT_DIR, "Overlap57peaksRNAlinks.csv"),    row.names = FALSE)
  write.csv(overlap_full,        file.path(OUT_DIR, "Overlap57peaksFULLannotation.csv"), row.names = FALSE)
  write.csv(gene_summary,        file.path(OUT_DIR, "Overlap57peaksgenesummary.csv"), row.names = FALSE)

  # LinkPeaks: all DA-gained peaks (full gene universe)
  peak_gene_da <- closest_gained %>%
    select(peak_id = query_region, gene_name) %>%
    filter(!is.na(gene_name))

  write.csv(peak_gene_da,
            file.path(OUT_DIR, "AstroDAgeneswithstats.csv"),
            row.names = FALSE)

  genes_use_all <- unique(peak_gene_da$gene_name)
  genes_use_all <- genes_use_all[!is.na(genes_use_all)]

  if (length(genes_use_all) > 0) {
    astrolink <- LinkPeaks(
      object           = astrolink,
      peak.assay       = "ATAC",
      expression.assay = "SCT",
      peak.slot        = "counts",
      expression.slot  = "data",
      genes.use        = genes_use_all,
      distance         = 500000,
      min.cells        = 10,
      method           = "pearson"
    )

    links_all    <- Links(astrolink)
    links_all_df <- data.frame(gene   = links_all$gene,
                                score  = links_all$score,
                                pvalue = links_all$pvalue,
                                peak   = links_all$peak)

    write.csv(links_all_df,
              file.path(OUT_DIR, "AllDAgainedLinkPeakslinks.csv"),
              row.names = FALSE)
  }

  # Top 50 DA genes for manuscript table
  top50_manuscript <- peak_gene_da %>%
    left_join(dagained %>% select(peak = peak, avg_log2FC, p_val_adj, pct.1, pct.2),
              by = c("peak_id" = "peak")) %>%
    filter(!is.na(p_val_adj), !grepl("^LOC", gene_name)) %>%
    arrange(p_val_adj) %>%
    head(50)

  write.csv(top50_manuscript,
            file.path(OUT_DIR, "AstroDATop50genesmanuscript.csv"),
            row.names = FALSE)
}


################################################################################
# 13. FIGURE 2 — GENE INTEGRATION PLOT
#
# Scatter plots correlating HC-specific promoter loop count with RNA
# log2FC (p2a) and -log10(p-value) (p2b). Top 10 genes labeled.
# Output: FIGURES/Fig2GeneIntegration.pdf (12 × 5.5 inches)
################################################################################

if (exists("hc_gene_rank_rna") && nrow(hc_gene_rank_rna) > 0) {

  plot_df <- hc_gene_rank_rna %>%
    filter(!is.na(avg_log2FC)) %>%
    mutate(
      neglog10p = -log10(p_val + 1e-300),
      sig       = p_val_adj < 0.05
    )

  top_genes <- plot_df %>%
    arrange(desc(nhclinks), p_val) %>%
    slice_head(n = 10)

  # Pearson correlation test for caption
  wilcox_res <- tryCatch(
    wilcox.test(nhclinks ~ sig, data = plot_df),
    error = function(e) list(p.value = NA)
  )

  p2a <- ggplot(plot_df, aes(x = avg_log2FC, y = nhclinks)) +
    geom_point(aes(color = sig), alpha = 0.7, size = 1.5) +
    scale_color_manual(
      values = c("TRUE" = "firebrick", "FALSE" = "grey40"),
      labels = c("TRUE" = "FDR < 0.05", "FALSE" = "n.s.")
    ) +
    geom_text_repel(data = top_genes, aes(label = promoter_gene),
                    size = 3, max.overlaps = 15, na.rm = TRUE) +
    stat_cor(method = "pearson", label.x.npc = 0.05, label.y.npc = 0.95,
             size = 4) +
    labs(x     = "RNA log2FC (Heat Call / Control)",
         y     = "HC-specific promoter-anchored links",
         color = "RNA significance") +
    theme(legend.position = c(0.8, 0.15),
          legend.background = element_rect(fill = "white", color = "black"))

  p2b <- ggplot(plot_df, aes(x = neglog10p, y = nhclinks)) +
    geom_point(aes(color = sig), alpha = 0.7, size = 1.5) +
    scale_color_manual(
      values = c("TRUE" = "firebrick", "FALSE" = "grey40")
    ) +
    geom_text_repel(data = top_genes, aes(label = promoter_gene),
                    size = 3, max.overlaps = 15, na.rm = TRUE) +
    labs(x       = "-log10(p-value) RNA",
         y       = "HC-specific promoter-anchored links",
         caption = paste0("Wilcoxon high- vs low-link genes p = ",
                          format(wilcox_res$p.value, digits = 3))) +
    theme(legend.position = "none",
          plot.caption     = element_text(size = 9))

  ggsave(
    file.path(DIR_FIG, "Fig2GeneIntegration.pdf"),
    p2a + p2b,
    width  = 12,
    height = 5.5
  )
}


################################################################################
# 14. SUMMARY
################################################################################

cat("SCRIPT 10 COMPLETE\n")
cat("HC-specific rewired connections :", nrow(rewiredhc),   "\n")
cat("CTRL-specific rewired connections:", nrow(rewiredctrl), "\n")
cat("Dual-hit peaks (HC-rewired ∩ DA-gained):", length(hcboth), "\n")
cat("Promoter-anchored HC links      :", nrow(rewiredhc_prom), "\n")
cat("Genes ranked by HC loops        :", nrow(hc_gene_rank),   "\n")
cat("Next: Run 11_Figure5_Locus_Plots.R\n")

