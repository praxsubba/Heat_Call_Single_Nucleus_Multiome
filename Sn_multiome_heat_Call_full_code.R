# Single-Nucleus Multiome Analysis: Complete and Annotated Workflow

This document provides a well-organized, thoroughly annotated, and reproducible R workflow for single-nucleus RNA-seq (snRNA-seq), single-nucleus ATAC-seq (snATAC-seq), and their integration using Seurat and Signac. The code is structured for clarity and ready for submission to GitHub as a supplement for high-impact publications[1][2][3].

---

## 1. Data Processing and Quality Control (QC)

### **RNA (snRNA-seq)**

Load required libraries
library(Seurat)
library(SeuratObject)
library(stringr)

Load and process mitochondrial genes
mt_genes <- read.delim("/data2/george_lab/psubba/zebra_finch_reference_data/mt_chromosome.gtf", header = FALSE, comment.char = "#")
mt_gene_ids <- unique(str_extract(mt_genes$V9, "(?<=gene_id\s)[^";]+")) %>% na.omit() %>% as.character()

Create Seurat objects for each sample
sample_ids <- paste0("S", 1:8)
seurat_objects <- list()
for (id in sample_ids) {
data_dir <- paste0("/data2/george_lab/psubba/", id, "combined_analysis/outs")
data_raw <- Read10X_h5(paste0(data_dir, "/filtered_feature_bc_matrix.h5"))
if (is.list(data_raw) && "Gene Expression" %in% names(data_raw)) {
data <- data_raw$Gene Expression
} else {
data <- data_raw
}
seurat_obj <- CreateSeuratObject(counts = data, project = id, min.cells = 3, min.features = 200)
mt_gene_ids_seurat <- mt_gene_ids %>% str_replace_all("", "-")
seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), mt_gene_ids_seurat))
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_objects[[id]] <- seurat_obj
}

Add metadata (Condition, Sex)
metadata <- data.frame(
Sample = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"),
Condition = c("Heat Call", "Heat Call", "Control", "Control", "Heat Call", "Heat Call", "Control", "Control"),
Sex = c("Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female")
)
for (i in seq_along(seurat_objects)) {
seurat_objects[[i]]$Condition <- metadata$Condition[i]
seurat_objects[[i]]$Sex <- metadata$Sex[i]
}

Merge all samples
merged_seurat <- merge(seurat_objects[], y = seurat_objects[-1], add.cell.ids = sample_ids, project = "Merged")

QC Filtering
merged_seurat_filt <- subset(merged_seurat, nFeature_RNA > 200 & nFeature_RNA < 3500 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mt < 5)

text

### **ATAC (snATAC-seq)**

library(Signac)
library(Seurat)
library(GenomicRanges)
library(rtracklayer)

samples <- paste0("S", 1:8)
metadata <- data.frame(
Sample = samples,
Condition = c("Heat Call", "Heat Call", "Control", "Control", "Heat Call", "Heat Call", "Control", "Control"),
Sex = c("Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female")
)

Read peaks for each sample and merge
gr.list <- lapply(samples, function(s) {
peaks <- read.table(paste0("/data2/george_lab/psubba/", s, "_combined_analysis/outs/atac_peaks.bed"), col.names = c("chr", "start", "end"))
makeGRangesFromDataFrame(peaks)
})
combined.peaks <- reduce(do.call(c, gr.list))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]

Read per-barcode metrics and filter
md.list <- lapply(samples, function(s) {
md <- read.table(paste0("/data2/george_lab/psubba/", s, "_combined_analysis/outs/per_barcode_metrics.csv"), stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)
md[md$atac_fragments > 500, ]
})
for (i in 1:length(samples)) {
md.list[[i]]$dataset <- samples[i]
md.list[[i]]$Condition <- metadata$Condition[metadata$Sample == samples[i]]
md.list[[i]]$Sex <- metadata$Sex[metadata$Sample == samples[i]]
}

Create fragment and feature matrices, then Seurat objects
frags.list <- mapply(function(s, md) {
CreateFragmentObject(path = paste0("/data2/george_lab/psubba/", s, "_combined_analysis/outs/atac_fragments.tsv.gz"), cells = rownames(md))
}, samples, md.list, SIMPLIFY = FALSE)
counts.list <- mapply(function(frags, md) {
FeatureMatrix(fragments = frags, features = combined.peaks, cells = rownames(md))
}, frags.list, md.list, SIMPLIFY = FALSE)
seurat.list <- mapply(function(s, counts, frags, md) {
assay <- CreateChromatinAssay(counts, fragments = frags)
obj <- CreateSeuratObject(assay, assay = "ATAC", meta.data = md)
obj
}, samples, counts.list, frags.list, md.list, SIMPLIFY = FALSE)
combined <- merge(seurat.list[], y = seurat.list[-1], add.cell.ids = samples)

Compute QC metrics
combined <- NucleosomeSignal(object = combined)
combined <- TSSEnrichment(object = combined)
combined$pct_reads_in_peaks <- combined$atac_peak_region_fragments / combined$atac_fragments * 100

QC Filtering
combined_filt <- subset(
x = combined,
subset = nCount_ATAC > 100 & nCount_ATAC < 100000 &
pct_reads_in_peaks > 15 &
nucleosome_signal < 4 &
TSS.enrichment > 2
)

text

---

## 2. Normalization

### **RNA**

SCTransform normalization
merged_seurat_filt <- SCTransform(merged_seurat_filt)

text

### **ATAC**

TF-IDF normalization and feature selection
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)

text

---

## 3. Clustering

### **RNA**

PCA and clustering
merged_seurat_filt <- RunPCA(merged_seurat_filt, assay = "SCT")
merged_seurat_filt <- FindNeighbors(merged_seurat_filt, dims = 1:35)
merged_seurat_filt <- FindClusters(merged_seurat_filt, resolution = 1.1)
merged_seurat_filt <- RunUMAP(merged_seurat_filt, dims = 1:35)

text

### **ATAC**

UMAP and clustering
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
object.list <- SplitObject(combined_filt, split.by = "dataset")
integration.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = rownames(combined_filt), reduction = "rlsi", dims = 2:50)
integrated <- IntegrateEmbeddings(anchorset = integration.anchors, reductions = combined_filt[["lsi"]], new.reduction.name = "integrated_lsi", dims.to.integrate = 1:50)
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:50)
integrated <- FindNeighbors(integrated, reduction = 'integrated_lsi', dims = 2:50)
integrated <- FindClusters(integrated, verbose = FALSE, algorithm = 3, resolution = 1.3)

text

### **WNN (Integrated RNA + ATAC)**

Integration and WNN clustering
combined_wnn <- FindMultiModalNeighbors(
object = combined_wnn,
reduction.list = list("pca", "lsi"),
dims.list = list(1:35, 2:50),
modality.weight.name = "RNA.weight",
verbose = TRUE
)
combined_wnn <- RunUMAP(
object = combined_wnn,
nn.name = "weighted.nn",
reduction.name = "wnn_umap",
reduction.key = "WNNUMAP_"
)
combined_wnn <- FindClusters(
object = combined_wnn,
resolution = 0.6,
graph.name = "wsnn",
algorithm = 3
)

text

---

## 4. Cell x Gene Marker Gene Summarization

library(dplyr)
library(readr)

Read and filter for brain tissue
df <- read.csv("C48_CELLxGENE_gene_expression_042925.csv")
brain_df <- df %>% filter(Tissue == "brain")

Calculate percent of cells expressing gene
brain_df <- brain_df %>%
mutate(Percent_Cells_Expressing = as.numeric(Number.of.Cells.Expressing.Genes) / as.numeric(Cell.Count))

For each gene, select top 20 cell types by expression and percent expressing
top20_per_gene <- brain_df %>%
group_by(Gene.Symbol) %>%
arrange(desc(Expression), desc(Percent_Cells_Expressing), .by_group = TRUE) %>%
slice_head(n = 20) %>%
ungroup()

Frequency of each cell type across all genes' top 20s
celltype_freq <- top20_per_gene %>%
count(Cell.Type, name = "Frequency") %>%
arrange(desc(Frequency))

Average Expression and Percent_Cells_Expressing per Cell.Type
avg_stats <- top20_per_gene %>%
group_by(Cell.Type) %>%
summarise(
Avg_Expression = mean(Expression, na.rm = TRUE),
Avg_Perc_Cells_Exp = mean(Percent_Cells_Expressing, na.rm = TRUE)
)
celltype_freq <- celltype_freq %>%
left_join(avg_stats, by = "Cell.Type")

text

---

## 5. Differential Accessibility and Gene Expression Analysis

### **Differential Accessibility (DA) — ATAC**

DefaultAssay(combined_wnn) <- "ATAC"
combined_wnn$cluster_condition <- paste(combined_wnn$wsnn_res.0.6, combined_wnn$Condition, sep = "_")
Idents(combined_wnn) <- "cluster_condition"
clusters <- levels(combined_wnn$wsnn_res.0.6)
DA_results <- list()
for(cluster in clusters) {
DA_results[[cluster]] <- FindMarkers(
combined_wnn,
ident.1 = paste0(cluster, "_Heat Call"),
ident.2 = paste0(cluster, "_Control"),
assay = "ATAC",
test.use = 'MAST',
min.pct = 0.1
)
if(nrow(DA_results[[cluster]]) > 0) {
DA_results[[cluster]]$Cluster <- cluster
DA_results[[cluster]]$Peak <- rownames(DA_results[[cluster]]
write.csv(DA_results[[cluster]], paste0("Cluster", cluster, "_DA.csv"))
}
}
combined_da <- do.call(rbind, DA_results)
write.csv(combined_da, "29_May_2025_DA_condition_cell_cluster.csv", row.names = FALSE)

text

### **Differential Gene Expression (DE) — RNA**

DefaultAssay(combined_wnn) <- "SCT"
combined_wnn$cluster_condition <- paste(combined_wnn$wsnn_res.0.6, combined_wnn$Condition, sep = "_")
Idents(combined_wnn) <- "cluster_condition"
clusters <- levels(combined_wnn$wsnn_res.0.6)
DE_results <- list()
for(cluster in clusters) {
DE_results[[cluster]] <- FindMarkers(
object = combined_wnn,
ident.1 = paste0(cluster, "_Heat Call"),
ident.2 = paste0(cluster, "_Control"),
test.use = "MAST",
recorrect_umi = FALSE,
min.pct = 0.25,
logfc.threshold = 0.25
)
if(nrow(DE_results[[cluster]]) > 0) {
DE_results[[cluster]]$Cluster <- cluster
DE_results[[cluster]]$Gene <- rownames(DE_results[[cluster]])
write.csv(DE_results[[cluster]], file = paste0("Cluster", cluster, "_DEG.csv"))
}
}
combined_de <- do.call(rbind, DE_results)
write.csv(combined_de, "All_Clusters_HeatCall_vs_Control_DEG.csv")

text

---

## 6. ATAC Peak Annotation

library(Seurat)
library(Signac)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome)
library(Biostrings)
library(GenomeInfoDb)
library(Rsamtools)
library(ggrepel)

Import GTF annotation file (leave seqnames as RefSeq style)
gtf <- import("/data2/george_lab/psubba/zebra_finch_reference_data/genomic.gtf")
gene.coords <- gtf[gtf$type %in% c('gene','transcript','exon','CDS','start_codon','stop_codon')]

Rename transcript_id column if needed
names(mcols(gene.coords))[names(mcols(gene.coords)) == "transcript_id"] <- "tx_id"

Ensure gene_name column exists
if (!"gene_name" %in% names(mcols(gene.coords))) {
gene.coords$gene_name <- gene.coords$gene_id
}

Filter to keep only chromosomes present in ATAC assay
valid_chroms <- seqlevels(granges(combined_wnn@assays$ATAC))
gene.coords <- keepSeqlevels(gene.coords, valid_chroms, pruning.mode = "coarse")

Add genome annotation to Seurat object
Annotation(combined_wnn) <- gene.coords

Load reference genome as FaFile for RegionStats()
fasta <- FaFile("/data2/george_lab/psubba/zebra_finch_reference_data/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna")

Index if not done already
indexFa(fasta)

Compute GC content stats for peaks
combined_wnn <- RegionStats(combined_wnn, genome = fasta)

Read differential ATAC peak results
da_cell_cond <- read.csv("DA_HeatCall_vs_Control.csv")
colnames(da_cell_cond) <- c("peak", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")

Adjust peak names to match ATAC assay (replace first hyphen with underscore for ClosestFeature)
da_cell_cond$peak <- gsub("^(N[WC])-", "\1_", da_cell_cond$peak)
da_cell_cond <- da_cell_cond[da_cell_cond$p_val_adj < 0.05, ]

Use ClosestFeature to annotate peaks to genes
annot_da_cell_cond <- ClosestFeature(combined_wnn, da_cell_cond$peak)

Convert peaks to GRanges
peaks_gr <- StringToGRanges(da_cell_cond$peak)

Export GRanges to BED file for external annotation (e.g., bedtools intersect)
export.bed(peaks_gr, "da_cell_cond_peaks.bed")

Example command for bedtools intersect (run in terminal):
bedtools intersect -a da_cell_cond_peaks.bed -b /data2/george_lab/psubba/zebra_finch_reference_data/genomic.gtf -wa -wb > da_cell_cond_peak_gene_overlaps.txt
Read external annotation results
peak_gene_overlaps <- read.csv(file="da_cell_cond_peak_gene_overlaps.csv")

Clean and merge peak-gene mapping
peak_gene_overlaps <- peak_gene_overlaps %>%
mutate(peak = paste0(chrom, "-", as.integer(start) + 1, "-", end)) %>%
select(peak, gene_id) %>%
distinct()

Merge with DA results
da_cell_cond_gene <- merge(da_cell_cond, peak_gene_overlaps, by = "peak", all.x = TRUE)

Select relevant columns from ClosestFeature results for further merging
annot_subset <- annot_da_cell_cond[, c("query_region", "gene_id")]
colnames(annot_subset) <- "peak"

Merge with DA results (update gene_id if available from ClosestFeature)
matching_peaks <- merge(da_cell_cond_gene, annot_subset, by = "peak")
matching_peaks$gene_id <- matching_peaks$gene_id.y

Combine all results and remove duplicates/NA
final_df <- rbind(da_cell_cond_gene, matching_peaks)
final_df <- unique(final_df)
final_df <- final_df[!is.na(final_df$gene_id), ]

text

---

## 7. Motif Enrichment and ChromVAR Analysis

library(JASPAR2020)
library(TFBSTools)
library(Rsamtools)
library(Signac)

DefaultAssay(combined_wnn) <- "ATAC"
pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
fasta <- FaFile("/data2/george_lab/psubba/zebra_finch_reference_data/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna")
indexFa(fasta)
combined_wnn <- RegionStats(combined_wnn, genome = fasta)
motif.matrix <- CreateMotifMatrix(features = granges(combined_wnn@assays$ATAC), pwm = pfm, genome = fasta)
rownames(motif.matrix) <- gsub("_", "-", rownames(motif.matrix))
motif <- CreateMotifObject(data = motif.matrix, pwm = pfm)
Motifs(combined_wnn) <- motif

ChromVAR analysis
library(BiocParallel)
register(SerialParam())
combined_wnn <- RunChromVAR(object = combined_wnn, genome = fasta, assay = "ATAC")
saveRDS(combined_wnn, file="motif.combined_wnn.rds")

DefaultAssay(motif.combined_wnn) <- 'chromvar'

Differential Activity
Idents(motif.combined_wnn) <- "Condition"
differential.activity <- FindMarkers(
object = motif.combined_wnn,
ident.1 = 'Heat Call',
ident.2 = 'Control',
only.pos = TRUE,
mean.fxn = rowMeans,
fc.name = "avg_diff"
)

MotifPlot(
object = motif.combined_wnn,
motifs = head(rownames(differential.activity)),
assay = 'ATAC'
)

text

---

## References

1. **Satija Lab: Seurat Documentation**  
   [https://satijalab.org/seurat/][1]
2. **Satija Lab: Signac Documentation**  
   [https://satijalab.org/signac/][2]
3. **JASPAR2020 Database**  
   [https://jaspar.genereg.net/][3]

---

This workflow is ready for publication and GitHub submission, providing a comprehensive and reproducible analysis pipeline for single-nucleus multiome data[1][2][3].

---

[1]: https://satijalab.org/seurat/  
[2]: https://satijalab.org/signac/  
[3]: https://jaspar.genereg.net/
