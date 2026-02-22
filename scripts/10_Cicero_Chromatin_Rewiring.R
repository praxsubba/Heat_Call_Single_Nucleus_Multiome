# ==============================================================================
# 10_Cicero_Chromatin_Rewiring.R
# ------------------------------------------------------------------------------
# Promoter Anchoring, RNA Integration & Chromatin–RNA Concordance
# (HC-specific rewired connections → gene ranking → RNA DE validation)
# ==============================================================================
# Depends on:  Script 01 outputs — rewiredhc.csv, rewiredctrl.csv,
#                                   connsclusterannotated.csv
#              29May2025DAconditioncellcluster.RData (Seurat multiome object)
# Outputs:     HCpromotergenerank.csv
#              HCpromotergenerankwithRNA.csv
#              rewiredhcpromoteranchored.csv
#              Figures: SupplFig2_RewiredGenesVsRNA.pdf
#                       Fig2GeneIntegration.pdf
#              ManuscriptTables/Table1GlobalStats.csv
#              ManuscriptTables/Table3TopRewiringGenes.csv
#              GOenrichmentLists/ (8 gene-list files)
# ==============================================================================


# ── 0. LIBRARIES ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(IRanges)
  library(tidyverse)    # dplyr, tidyr, ggplot2, readr, stringr
  library(ggpubr)       # CORRECTION: required for stat_cor() in scatter plot
  library(ggrepel)      # geom_text_repel()
})
set.seed(1234)
options(stringsAsFactors = FALSE)


# ── 1. CONFIGURATION ──────────────────────────────────────────────────────────
RUNDIR     <- "data2georgelabpsubbacicerorunszfinchcicero20251213154435"
OUTDIR     <- file.path(RUNDIR, "rewiringAstrocyteHCvsCTRL")
FIGURES    <- file.path(OUTDIR, "figures")
MANU       <- file.path(OUTDIR, "ManuscriptTables")
GO_DIR     <- file.path(OUTDIR, "GOenrichmentLists")
RDATA_PATH <- "29May2025DAconditioncellcluster.RData"

# CORRECTION: DEG_CSV removed entirely — RNA DE is run fresh below using
# paste0(cluster_annotation, Condition) idents, consistent with Script 05.

# Thresholds
COACCESS_CUTOFF  <- 0.15   # Cicero co-accessibility threshold
PROMOTER_WINDOW  <- 2000   # bp upstream / downstream of TSS

dir.create(OUTDIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)
dir.create(MANU,    recursive = TRUE, showWarnings = FALSE)
dir.create(GO_DIR,  recursive = TRUE, showWarnings = FALSE)

cat("Script: 10_Cicero_Chromatin_Rewiring.R\n")
cat("RUNDIR:", RUNDIR, "\n")
cat("OUTDIR:", OUTDIR, "\n")


# ── 2. LOAD DATA ──────────────────────────────────────────────────────────────
cat("\n=== STEP 1: LOADING DATA ===\n")

load(RDATA_PATH)
stopifnot(exists("combined_wnn"))

hcedges <- read_csv(file.path(OUTDIR, "rewiredhc.csv"),
                    show_col_types = FALSE)
ctedges <- read_csv(file.path(OUTDIR, "rewiredctrl.csv"),
                    show_col_types = FALSE)
conns   <- read_csv(file.path(OUTDIR, "connsclusterannotated.csv"),
                    show_col_types = FALSE)

cat("Loaded:\n")
cat(" -", nrow(conns),   "chromatin connections\n")
cat(" -", nrow(hcedges), "HC-specific edges\n")
cat(" -", nrow(ctedges), "Control-specific edges\n")

# Unique peak names (character vector; used downstream for overlap tests)
rewiredhc       <- unique(c(hcedges$Peak1, hcedges$Peak2))
hc_unique_peaks <- length(rewiredhc)


# ── 3. SUBSET TO ASTROCYTES ───────────────────────────────────────────────────
cat("\n=== STEP 2: SUBSETTING TO ASTROCYTES ===\n")

obj <- subset(combined_wnn, subset = cluster_annotation == "Astrocyte")
DefaultAssay(obj) <- "ATAC"

stopifnot("cluster_annotation" %in% colnames(obj@meta.data))
stopifnot("Condition"           %in% colnames(obj@meta.data))

cat("Astrocyte cells:  ", ncol(obj), "\n")
cat("Control cells:    ", sum(obj$Condition == "Control"), "\n")
cat("Heat Call cells:  ", sum(obj$Condition == "Heat Call"), "\n")


# ── 4. HELPER: peak string → GRanges ─────────────────────────────────────────
# Handles both NC-044211.2-start-end (4 parts) and chr-start-end (3 parts)
peak_to_gr <- function(pk) {
  p <- strsplit(pk, "-")[[1]]
  if (length(p) == 4) {
    chr   <- paste0(p[1], "-", p[2])   # reconstruct: NC-044211.2
    start <- as.integer(p[3])
    end   <- as.integer(p[4])
  } else if (length(p) == 3) {
    chr   <- p[1]
    start <- as.integer(p[2])
    end   <- as.integer(p[3])
  } else {
    stop("Unexpected peak format: ", pk)
  }
  chr <- gsub("-", "_", chr)           # NC-044211.2 → NC_044211.2 (GTF format)
  GRanges(seqnames = chr,
          ranges   = IRanges(start = start, end = end),
          peakid   = pk)
}


# ── 5. BUILD PROMOTER COORDINATES FROM GTF ANNOTATION ────────────────────────
cat("\n=== STEP 3: BUILDING PROMOTER COORDINATES ===\n")

gtf_annotation <- Annotation(obj)
if (is.null(gtf_annotation)) stop("No gene annotation found in Seurat object.")

ann_gene <- gtf_annotation[gtf_annotation$type == "gene"]
ann_gene <- ann_gene[!duplicated(ann_gene$gene_name)]   # one entry per gene

# Strand-aware TSS → ± PROMOTER_WINDOW
tss <- ifelse(as.character(strand(ann_gene)) == "+",
              start(ann_gene), end(ann_gene))

promoters_gr <- GRanges(
  seqnames = seqnames(ann_gene),
  ranges   = IRanges(start = tss - PROMOTER_WINDOW,
                     end   = tss + PROMOTER_WINDOW)
)
mcols(promoters_gr)$gene_name <- ann_gene$gene_name
cat("Total promoters defined:", length(promoters_gr), "\n")


# ── 6. PROMOTER ANCHORING ─────────────────────────────────────────────────────
cat("\n=== STEP 4: PROMOTER ANCHORING (HC-REWIRED CONNECTIONS) ===\n")

p1gr <- do.call(c, lapply(hcedges$Peak1, peak_to_gr))
p2gr <- do.call(c, lapply(hcedges$Peak2, peak_to_gr))

cat("Example peak:   ", hcedges$Peak1[1], "\n")
cat("Converted chr:  ", as.character(seqnames(p1gr)[1]), "\n")

p1hit <- findOverlaps(p1gr, promoters_gr, ignore.strand = TRUE)
p2hit <- findOverlaps(p2gr, promoters_gr, ignore.strand = TRUE)

cat("Peak1 promoter overlaps:", length(p1hit), "\n")
cat("Peak2 promoter overlaps:", length(p2hit), "\n")

# Assign gene names to each endpoint
hcedges$Peak1gene <- NA_character_
hcedges$Peak2gene <- NA_character_
hcedges$Peak1gene[queryHits(p1hit)] <- promoters_gr$gene_name[subjectHits(p1hit)]
hcedges$Peak2gene[queryHits(p2hit)] <- promoters_gr$gene_name[subjectHits(p2hit)]

# FIX: use pivot_longer instead of coalesce().
# coalesce(Peak1gene, Peak2gene) silently drops Peak2gene when Peak1gene is
# already non-NA — meaning promoter–promoter connections are counted only once.
# pivot_longer preserves both endpoints as separate rows, then distinct() removes
# true duplicates (same Peak1–Peak2 pair mapping to the same gene via both ends).
rewiredhc_prom <- hcedges %>%
  filter(!is.na(Peak1gene) | !is.na(Peak2gene)) %>%
  tidyr::pivot_longer(
    cols      = c(Peak1gene, Peak2gene),
    names_to  = "endpoint",
    values_to = "promoter_gene"
  ) %>%
  filter(!is.na(promoter_gene)) %>%
  distinct(Peak1, Peak2, coaccess, promoter_gene, .keep_all = TRUE)

cat("Promoter-anchored HC links:", nrow(rewiredhc_prom), "\n")


# ── 7. RANK GENES BY HC-SPECIFIC PROMOTER LOOPS ──────────────────────────────
cat("\n=== STEP 5: GENE RANKING BY REWIRING ===\n")

hc_gene_rank <- rewiredhc_prom %>%
  group_by(promoter_gene) %>%
  summarise(
    nhclinks      = n(),
    mean_coaccess = mean(coaccess,   na.rm = TRUE),
    mean_p1_hc    = mean(Peak1acchc, na.rm = TRUE),
    mean_p2_hc    = mean(Peak2acchc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(nhclinks), desc(mean_coaccess))

cat("Genes with HC-specific promoter links:", nrow(hc_gene_rank), "\n")
cat("Top 15 genes:\n")
print(head(hc_gene_rank, 15))

write.csv(hc_gene_rank,
          file.path(OUTDIR, "HCpromotergenerank.csv"),
          row.names = FALSE)
write.csv(rewiredhc_prom,
          file.path(OUTDIR, "rewiredhcpromoteranchored.csv"),
          row.names = FALSE)


# ── 8. RNA DIFFERENTIAL EXPRESSION — ASTROCYTES ───────────────────────────────
cat("\n=== STEP 6: RNA DE — ASTROCYTES (HC vs CTRL) ===\n")

# CORRECTION: Do NOT load a pre-computed DEG CSV.
# Create combined cluster+condition idents (matching Script 05 approach) and
# run FindMarkers fresh on combined_wnn so the DE test is fully reproducible
# and always reflects the current object's cell composition.
stopifnot("SCT" %in% names(combined_wnn@assays))

Idents(combined_wnn) <- paste0(combined_wnn$cluster_annotation,
                                combined_wnn$Condition)

ident1 <- "AstrocyteHeat Call"
ident2 <- "AstrocyteControl"

if (!ident1 %in% levels(Idents(combined_wnn)) |
    !ident2 %in% levels(Idents(combined_wnn))) {
  stop("Could not find idents '", ident1, "' and/or '", ident2,
       "' — check cluster_annotation and Condition labels in meta.data.")
}

rna_da <- FindMarkers(
  object          = combined_wnn,
  ident.1         = ident1,
  ident.2         = ident2,
  assay           = "SCT",
  logfc.threshold = 0,
  min.pct         = 0.1
)
rna_da$gene <- rownames(rna_da)

cat("DE genes (Astrocyte HC vs CTRL, any FC, pct ≥ 0.1):", nrow(rna_da), "\n")


# ── 9. JOIN PROMOTER RANK WITH RNA ────────────────────────────────────────────
cat("\n=== STEP 7: CHROMATIN–RNA INTEGRATION ===\n")

hc_gene_rank_rna <- hc_gene_rank %>%
  left_join(
    rna_da %>% select(gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2),
    by = c("promoter_gene" = "gene")
  )

cat("Genes with RNA data:", sum(!is.na(hc_gene_rank_rna$p_val)), "\n")

write.csv(hc_gene_rank_rna,
          file.path(OUTDIR, "HCpromotergenerankwithRNA.csv"),
          row.names = FALSE)


# ── 10. CHROMATIN–RNA CONCORDANCE STATISTICS ─────────────────────────────────
cat("\n=== STEP 8: CONCORDANCE STATISTICS ===\n")

rna_test <- hc_gene_rank_rna %>%
  filter(!is.na(avg_log2FC)) %>%
  mutate(
    high_links   = nhclinks > median(nhclinks, na.rm = TRUE),
    rewire_group = case_when(
      nhclinks >  10 ~ "High (>10)",
      nhclinks >=  5 ~ "Moderate (5-9)",
      TRUE           ~ "Low (1-4)"
    )
  )

cor_result     <- cor.test(rna_test$nhclinks, rna_test$avg_log2FC,
                           method = "pearson")
wilcox_result  <- wilcox.test(avg_log2FC ~ high_links, data = rna_test)
kruskal_result <- kruskal.test(avg_log2FC ~ rewire_group, data = rna_test)

cat(sprintf("Pearson r  (links vs RNA log2FC):     %.3f, p = %.3e\n",
            cor_result$estimate, cor_result$p.value))
cat(sprintf("Wilcoxon p (high- vs low-link genes): %.3e\n",
            wilcox_result$p.value))
cat(sprintf("Kruskal-Wallis p (3-group):           %.3e\n",
            kruskal_result$p.value))
cat("(Weak/null correlation = chromatin PRIMING model — expected)\n")


# ── 11. FIGURE: SCATTER REWIRING vs RNA (SupplFig2) ──────────────────────────
cat("\n=== STEP 9: SUPPLEMENTAL FIGURE 2 — REWIRING vs RNA ===\n")

plot_df <- hc_gene_rank_rna %>%
  filter(!is.na(avg_log2FC)) %>%
  mutate(
    neg_log10p = -log10(p_val + 1e-300),
    sig        = p_val_adj < 0.05
  )

top_genes <- plot_df %>%
  arrange(desc(nhclinks), p_val) %>%
  slice_head(n = 10)

p_scatter <- ggplot(plot_df, aes(x = avg_log2FC, y = nhclinks)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 1.5) +
  scale_color_manual(
    values = c("TRUE" = "firebrick", "FALSE" = "grey40"),
    labels = c("TRUE" = "FDR < 0.05", "FALSE" = "n.s.")
  ) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.8) +
  geom_text_repel(data         = top_genes,
                  aes(label    = promoter_gene),
                  size         = 3,
                  max.overlaps = 15,
                  na.rm        = TRUE) +
  stat_cor(method      = "pearson",     # requires ggpubr
           label.x.npc = 0.05,
           label.y.npc = 0.95,
           size         = 4) +
  labs(
    x        = "RNA log2FC (Heat Call vs Control)",
    y        = "HC-specific promoter-anchored loops",
    color    = "RNA significance",
    title    = sprintf("Promoter rewiring vs RNA expression (%d genes)",
                       nrow(plot_df)),
    subtitle = sprintf("Pearson r = %.3f, p = %.3e",
                       cor_result$estimate, cor_result$p.value)
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position   = c(0.82, 0.15),
        legend.background = element_rect(fill = "white", color = "black"))

ggsave(file.path(FIGURES, "SupplFig2_RewiredGenesVsRNA.pdf"),
       p_scatter, width = 7, height = 5.5)
ggsave(file.path(FIGURES, "SupplFig2_RewiredGenesVsRNA.png"),
       p_scatter, width = 7, height = 5.5, dpi = 300)
cat("Saved SupplFig2_RewiredGenesVsRNA (.pdf + .png)\n")


# ── 12. FIGURE: GENE INTEGRATION (Fig2) ──────────────────────────────────────
cat("\n=== STEP 10: FIGURE 2 — GENE INTEGRATION ===\n")

p2a <- ggplot(plot_df, aes(x = avg_log2FC, y = nhclinks)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 1.5) +
  scale_color_manual(
    values = c("TRUE" = "firebrick", "FALSE" = "grey40"),
    labels = c("TRUE" = "FDR < 0.05", "FALSE" = "n.s.")
  ) +
  geom_text_repel(data         = top_genes,
                  aes(label    = promoter_gene),
                  size         = 3,
                  max.overlaps = 15,
                  na.rm        = TRUE) +
  stat_cor(method      = "pearson",
           label.x.npc = 0.05,
           label.y.npc = 0.95,
           size         = 4) +
  labs(x     = "RNA log2FC (HC vs Control)",
       y     = "HC-specific promoter-anchored loops",
       color = "RNA significance") +
  theme_bw(base_size = 12) +
  theme(legend.position   = c(0.82, 0.15),
        legend.background = element_rect(fill = "white", color = "black"))

p2b <- ggplot(plot_df, aes(x = neg_log10p, y = nhclinks)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 1.5) +
  scale_color_manual(
    values = c("TRUE" = "firebrick", "FALSE" = "grey40")
  ) +
  geom_text_repel(data         = top_genes,
                  aes(label    = promoter_gene),
                  size         = 3,
                  max.overlaps = 15,
                  na.rm        = TRUE) +
  labs(x       = "-log10(p-value, RNA)",
       y       = "HC-specific promoter-anchored loops",
       caption = paste0("Wilcoxon (high- vs low-link genes): p = ",
                        format(wilcox_result$p.value, digits = 3))) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.caption    = element_text(size = 9))

ggsave(file.path(FIGURES, "Fig2GeneIntegration.pdf"),
       p2a + p2b, width = 12, height = 5.5)
cat("Saved Fig2GeneIntegration.pdf\n")


# ── 13. HIGH-REWIRING GENE SUMMARY ───────────────────────────────────────────
cat("\n=== STEP 11: TOP REWIRED GENES ===\n")

top10 <- hc_gene_rank_rna %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(nhclinks), p_val) %>%
  select(promoter_gene, nhclinks, avg_log2FC, p_val, p_val_adj) %>%
  head(10)
cat("TOP 10 REWIRED GENES:\n")
print(top10)

high_rewire <- hc_gene_rank_rna %>%
  filter(nhclinks > 10) %>%
  arrange(desc(nhclinks)) %>%
  select(promoter_gene, nhclinks, avg_log2FC, p_val_adj)
cat(sprintf("\nGenes with >10 HC-specific loops: %d\n", nrow(high_rewire)))
if (nrow(high_rewire) > 0) print(high_rewire)

# Multi-criteria ranking for locus-plot candidates (Fig5)
top_candidates <- hc_gene_rank_rna %>%
  filter(!is.na(avg_log2FC)) %>%
  mutate(
    rank_links    = rank(-nhclinks,        ties.method = "average"),
    rank_pval     = rank( p_val,           ties.method = "average"),
    rank_fc       = rank(-abs(avg_log2FC), ties.method = "average"),
    combined_rank = rank_links + rank_pval + rank_fc,
    category      = case_when(
      nhclinks > 10 & p_val < 0.05 & avg_log2FC > 0 ~ "High-rewiring Upregulated",
      nhclinks > 10 & p_val < 0.05 & avg_log2FC < 0 ~ "High-rewiring Downregulated",
      nhclinks > 10 & p_val >= 0.05                  ~ "High-rewiring Primed",
      nhclinks <= 10 & p_val < 0.05                  ~ "Low-rewiring DE",
      TRUE                                            ~ "Other"
    )
  ) %>%
  arrange(combined_rank) %>%
  slice_head(n = 20)

write.csv(top_candidates,
          file.path(FIGURES, "Fig5LocusCandidatesTop20.csv"),
          row.names = FALSE)


# ── 14. MANUSCRIPT SUMMARY TABLE ─────────────────────────────────────────────
cat("\n=== STEP 12: MANUSCRIPT SUMMARY TABLE ===\n")

global_stats <- tibble(
  Metric = c(
    "Total astrocyte connections (Cicero)",
    "Heat call-specific rewired connections",
    "Control-specific rewired connections",
    "Unique peaks in HC-rewired connections",
    "Genes with >= 1 HC promoter-anchored loop",
    "Genes with >= 5 HC promoter-anchored loops",
    "Genes with > 10 HC promoter-anchored loops",
    "Pearson r (links vs RNA log2FC)",
    "Pearson p-value",
    "Wilcoxon p (high- vs low-link genes)",
    "Kruskal-Wallis p (3-group comparison)"
  ),
  Value = c(
    format(nrow(conns),   big.mark = ","),
    format(nrow(hcedges), big.mark = ","),
    format(nrow(ctedges), big.mark = ","),
    format(hc_unique_peaks, big.mark = ","),
    as.character(sum(hc_gene_rank_rna$nhclinks >= 1,  na.rm = TRUE)),
    as.character(sum(hc_gene_rank_rna$nhclinks >= 5,  na.rm = TRUE)),
    as.character(sum(hc_gene_rank_rna$nhclinks >  10, na.rm = TRUE)),
    as.character(round(cor_result$estimate, 3)),
    format(cor_result$p.value,     digits = 3),
    format(wilcox_result$p.value,  digits = 3),
    format(kruskal_result$p.value, digits = 3)
  )
)
print(global_stats)

write.csv(global_stats,
          file.path(MANU, "Table1GlobalStats.csv"), row.names = FALSE)

top_genes_ranked <- hc_gene_rank_rna %>%
  filter(!is.na(avg_log2FC)) %>%
  mutate(
    rank_links    = rank(-nhclinks,        ties.method = "average"),
    rank_pval     = rank( p_val,           ties.method = "average"),
    rank_fc       = rank(-abs(avg_log2FC), ties.method = "average"),
    combined_rank = rank_links + rank_pval + rank_fc
  ) %>%
  arrange(combined_rank) %>%
  slice_head(n = 50) %>%
  select(promoter_gene, nhclinks, avg_log2FC, p_val, p_val_adj, combined_rank)

write.csv(top_genes_ranked,
          file.path(MANU, "Table3TopRewiringGenes.csv"), row.names = FALSE)


# ── 15. GO ENRICHMENT GENE LISTS ─────────────────────────────────────────────
cat("\n=== STEP 13: GO ENRICHMENT GENE LISTS ===\n")

genes_any_rewiring  <- hc_gene_rank_rna %>% filter(nhclinks >= 1) %>%
  pull(promoter_gene) %>% unique() %>% sort()
genes_high_rewiring <- hc_gene_rank_rna %>% filter(nhclinks >= 5) %>%
  pull(promoter_gene) %>% unique() %>% sort()
genes_very_high     <- hc_gene_rank_rna %>% filter(nhclinks >  10) %>%
  pull(promoter_gene) %>% unique() %>% sort()
genes_rewiring_up   <- hc_gene_rank_rna %>%
  filter(nhclinks >= 5, !is.na(avg_log2FC), avg_log2FC > 0, p_val < 0.05) %>%
  pull(promoter_gene) %>% unique() %>% sort()
genes_rewiring_down <- hc_gene_rank_rna %>%
  filter(nhclinks >= 5, !is.na(avg_log2FC), avg_log2FC < 0, p_val < 0.05) %>%
  pull(promoter_gene) %>% unique() %>% sort()
genes_primed_only   <- hc_gene_rank_rna %>%
  filter(nhclinks >= 5, !is.na(avg_log2FC), p_val >= 0.05) %>%
  pull(promoter_gene) %>% unique() %>% sort()
top50_genes         <- top_genes_ranked$promoter_gene[
  seq_len(min(50, nrow(top_genes_ranked)))]
genes_de_only       <- hc_gene_rank_rna %>%
  filter(!is.na(avg_log2FC), p_val < 0.05) %>%
  pull(promoter_gene) %>% unique() %>% sort()

writeLines(genes_any_rewiring,  file.path(GO_DIR, "1AllRewiringGenes.txt"))
writeLines(genes_high_rewiring, file.path(GO_DIR, "2HighRewiringGenes5plus.txt"))
writeLines(genes_very_high,     file.path(GO_DIR, "3VeryHighRewiringGenes10plus.txt"))
writeLines(genes_rewiring_up,   file.path(GO_DIR, "4RewiringPlusUpregulated.txt"))
writeLines(genes_rewiring_down, file.path(GO_DIR, "5RewiringPlusDownregulated.txt"))
writeLines(genes_primed_only,   file.path(GO_DIR, "6PrimedOnlyNoRNAChange.txt"))
writeLines(top50_genes,         file.path(GO_DIR, "7Top50MultiCriteria.txt"))
writeLines(genes_de_only,       file.path(GO_DIR, "8DEGenesRNAOnly.txt"))

gene_list_summary <- tibble(
  File = c(
    "1AllRewiringGenes.txt",          "2HighRewiringGenes5plus.txt",
    "3VeryHighRewiringGenes10plus.txt","4RewiringPlusUpregulated.txt",
    "5RewiringPlusDownregulated.txt", "6PrimedOnlyNoRNAChange.txt",
    "7Top50MultiCriteria.txt",         "8DEGenesRNAOnly.txt"
  ),
  n_genes = c(
    length(genes_any_rewiring),  length(genes_high_rewiring),
    length(genes_very_high),     length(genes_rewiring_up),
    length(genes_rewiring_down), length(genes_primed_only),
    length(top50_genes),         length(genes_de_only)
  )
)
print(gene_list_summary)
write.csv(gene_list_summary,
          file.path(GO_DIR, "README_GeneLists.csv"), row.names = FALSE)


# ── 16. FINAL SUMMARY ────────────────────────────────────────────────────────
cat("\n========== 10_Cicero_Chromatin_Rewiring.R COMPLETE ==========\n")
cat("Outputs saved to:                 ", normalizePath(OUTDIR), "\n\n")
cat(sprintf("HC-rewired connections:            %d\n",  nrow(hcedges)))
cat(sprintf("Unique HC peaks:                   %d\n",  hc_unique_peaks))
cat(sprintf("Promoter-anchored HC links:        %d\n",  nrow(rewiredhc_prom)))
cat(sprintf("Genes with >= 1 HC loop:           %d\n",  nrow(hc_gene_rank)))
cat(sprintf("Genes with RNA data:               %d\n",  sum(!is.na(hc_gene_rank_rna$p_val))))
cat(sprintf("Genes with > 10 HC loops:          %d\n",  nrow(high_rewire)))
cat(sprintf("Pearson r (rewiring vs RNA):        %.3f\n", cor_result$estimate))
cat(sprintf("Pearson p:                          %.3e\n", cor_result$p.value))
cat(sprintf("Kruskal-Wallis p (3-group):         %.3e\n", kruskal_result$p.value))
cat("\nKey output files:\n")
cat("  - HCpromotergenerankwithRNA.csv    (main integration table)\n")
cat("  - HCpromotergenerank.csv           (gene-level rewiring rank)\n")
cat("  - rewiredhcpromoteranchored.csv    (edge-level promoter-anchored table)\n")
cat("  - figures/SupplFig2_RewiredGenesVsRNA.pdf\n")
cat("  - figures/Fig2GeneIntegration.pdf\n")
cat("  - ManuscriptTables/Table1GlobalStats.csv\n")
cat("  - ManuscriptTables/Table3TopRewiringGenes.csv\n")
cat("  - GOenrichmentLists/ (8 gene-list .txt files)\n")
cat("\nNext step: run 03TFactivityandenrichment.R\n")
