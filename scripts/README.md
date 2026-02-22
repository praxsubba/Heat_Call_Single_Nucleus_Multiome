# README.md `scripts`

# scripts/

Analysis pipeline for **Subba et al. 2026** —
*Multiome Profiling Reveals Astrocyte and Neuroendocrine Targets of
Prenatal Acoustic Programming in Zebra Finch Embryos*

Single-nucleus RNA-seq + ATAC-seq (10x Multiome) analysis of zebra finch
embryos exposed to heat call acoustic stimulation. All scripts are written
in R and designed to be run sequentially on an HPC cluster (SLURM).

---

## Pipeline Overview
Raw 10x Multiome data
│
▼
01 RNA QC ──► 02 ATAC QC ──► 03 WNN Integration ──► 04 Cell Type Annotation
│
┌───────────────────────────────┤
▼ ▼
05 Differential Gene 06 Differential Chromatin
Expression Accessibility
│ │
└──────────────┬────────────────┘
▼
07 Peak Annotation
│
┌──────────────┼──────────────┐
▼ ▼ ▼
08 ChromVAR 09 hdWGCNA 10 Cicero
Motif Analysis Networks Chromatin Rewiring
│
┌────────────┴──────────────┐
▼ ▼
11 TF Regulatory 12 Pseudotime
Network GRN Trajectory Analysis


---

## Scripts

### 01 · `01_RNA_QC_and_Processing.R`
RNA quality control, ambient RNA removal (SoupX), doublet detection
(DoubletFinder), normalization (SCTransform), and PCA.

| | |
|---|---|
| **Input** | Raw 10x Cellranger output (`filtered_feature_bc_matrix/`) |
| **Output** | `results/01_RNA/01_seurat_rna_qc.rds` |
| **Key packages** | Seurat, SoupX, DoubletFinder, scater |

---

### 02 · `02_ATAC_QC_and_Processing.R`
ATAC-seq quality control (TSS enrichment, nucleosome signal, FRiP),
peak calling, LSI dimensionality reduction.

| | |
|---|---|
| **Input** | Cellranger ATAC fragments + peaks; `results/01_RNA/` |
| **Output** | `results/02_ATAC/02_seurat_atac_qc.rds` |
| **Key packages** | Signac, Seurat, GenomicRanges |

---

### 03 · `03_WNN_Integration_and_Clustering.R`
Weighted Nearest Neighbor (WNN) joint RNA + ATAC integration,
UMAP embedding, and graph-based clustering across 42 WNN clusters.

| | |
|---|---|
| **Input** | `results/01_RNA/`, `results/02_ATAC/` |
| **Output** | `results/03_WNN/03_combined_wnn.rds` |
| **Key packages** | Seurat, Signac |

---

### 04 · `04_Cell_Type_Annotation.R`
Manual and marker-based cell type annotation of WNN clusters.
Assigns `cluster_annotation` metadata column used throughout the pipeline.

| | |
|---|---|
| **Input** | `results/03_WNN/03_combined_wnn.rds` |
| **Output** | `results/04_Annotation/04_combined_wnn_annotated.rds` |
| **Key packages** | Seurat, ggplot2 |

---

### 05 · `05_Differential_Gene_Expression.R`
Pseudobulk differential gene expression (DESeq2) between Heat Call and
Control conditions per cell type. Produces per-cluster DEG tables and
volcano plots.

| | |
|---|---|
| **Input** | `results/04_Annotation/` |
| **Output** | `results/05_DEG/` — per-cluster DEG CSVs, volcano PDFs |
| **Key packages** | Seurat, DESeq2, ggrepel |

---

### 06 · `06_Differential_Chromatin_Accessibility.R`
Differential chromatin accessibility analysis (logistic regression) between
Heat Call and Control conditions per cell type. Annotates peaks with
nearest gene and regulatory context.

| | |
|---|---|
| **Input** | `results/04_Annotation/` |
| **Output** | `results/06_DA/` — per-cluster DA peak CSVs, `06_combined_wnn_DA.rds` |
| **Key packages** | Signac, Seurat, GenomicRanges |

---

### 07 · `07_Peak_Annotation.R`
Peak genomic context annotation (promoter, intron, distal intergenic),
distance-to-TSS distributions, and peak–gene overlap analysis.

| | |
|---|---|
| **Input** | `results/06_DA/06_combined_wnn_DA.rds`; zebra finch GTF |
| **Output** | `results/07_PeakAnnotation/` — annotated peak tables, distribution PDFs |
| **Key packages** | ChIPseeker, TxDb (custom zebra finch), clusterProfiler |

---

### 08 · `08_ChromVAR_Motif_Analysis.R`
Transcription factor motif enrichment in ATAC-seq peaks using ChromVAR
and JASPAR 2022 motif database. Differential motif activity between
conditions per cell type.

| | |
|---|---|
| **Input** | `results/06_DA/06_combined_wnn_DA.rds`; `data/ZebraFinchMotifsWithTFGenes.rds` |
| **Output** | `results/08_ChromVAR/08_astroobj_chromvar.rds`; motif deviation PDFs |
| **Key packages** | chromVAR, JASPAR2022, motifmatchr, Signac |

---

### 09 · `09_hdWGCNA_Coexpression_Networks.R`
Weighted gene co-expression network analysis (hdWGCNA) on astrocyte
subset. Constructs metacell-based co-expression modules, computes module
eigengenes, performs module–trait correlation, and assigns TF regulons
(strategy A). Generates Figure 4 TF-module delta heatmap and network.

| | |
|---|---|
| **Input** | `results/08_ChromVAR/08_astroobj_chromvar.rds`; `AstrohdWGCNARNAWorkspaceNov2025.RData`; `AstroTFNetworkCheckpoint.rds`; `subclusters.rds` |
| **Output** | `results/09_hdWGCNA/` — module trait CSV, regulon score RDS files, differential regulon CSVs, Fig4A/B PDFs, Astro-M7 violin |
| **Key packages** | hdWGCNA, WGCNA, Seurat |
| ⚠️ **Note** | Network construction is computationally intensive — run on HPC. Pre-computed workspace `AstrohdWGCNARNAWorkspaceNov2025.RData` skips this step. |

---

### 10 · `10_Cicero_Chromatin_Rewiring.R`
Cicero co-accessibility network analysis to identify condition-specific
chromatin rewiring in astrocytes. Integrates rewired peaks with RNA
differential expression, performs promoter anchoring on HC-specific
connections, and ranks genes by dual chromatin + RNA evidence.

| | |
|---|---|
| **Input** | `29May2025DAconditioncellcluster.RData`; Cicero output CSVs (`rewiredhc.csv`, `rewiredctrl.csv`, `connsclusterannotated.csv`) |
| **Output** | `results/10_Cicero/` — `HCpromotergenerankwithRNA.csv`, `rewiredhcpromoteranchored.csv`, SupplFig2, Fig2, ManuscriptTables/, GOenrichmentLists/ |
| **Key packages** | Cicero, Signac, GenomicRanges, ggpubr |

---

### 11 · `11_TF_Regulatory_Network_GRN.R`
TF regulon GRN visualization and per-module deep dive, extending
Script 09. Loads pre-computed regulon scores and differential regulon
tables, performs Astro-M7 module-specific TF filtering (JASPAR TRUE-TF
validated), generates regulon bar plots, TF network plots, UMAP
visualizations, and GO enrichment gene lists. Produces the ASCL1 focal
network (Figure 3D) and TOM-scaffold GRN (Figure 3E).

| | |
|---|---|
| **Input** | `AstroTFNetworkCheckpoint.rds`; `subclusters.rds`; `results/09_hdWGCNA/` regulon score RDS + differential regulon CSV; `data/jaspar_mapping.rds` |
| **Output** | `results/11_TF_GRN/` — `0GlobalAnalysis/` plots, `Astro-M7Analysis/` per-TF figures and GO lists, `ASCL1_TFNetwork_strategyCdepth2Cor_labeled.pdf` (Fig 3D), `ASCL1_TOMScaffoldNetwork.pdf` (Fig 3E) |
| **Key packages** | hdWGCNA, Seurat, ggraph, igraph, ggrepel |

---

### 12 · `12_Pseudotime_Trajectory_Analysis.R`
Slingshot pseudotime trajectory inference on the WNN UMAP embedding,
followed by tradeSeq GAM fitting per condition to identify
pseudotime-dynamic genes. Per-cluster pseudotime shift tested by
Wilcoxon rank-sum. Smooth expression plots for key differentiation
regulators (ASCL1, DLL1, NOTCH1, MAML3, NFIA, FABP7) along normalized
pseudotime. Generates Figure 5A (trajectory UMAP) and Figure 5B
(smooth expression plots).

| | |
|---|---|
| **Input** | `29May2025DAconditioncellcluster.RData` |
| **Output** | `results/12_Pseudotime/` — `figures/Fig5A_TrajectoryUMAP.pdf`, `figures/Fig5B_SmoothExpressionPlots.pdf`, `tables/PerClusterPseudotimeShift.csv`, `tables/tradeSeq_AssociationTest_*.csv` |
| **Key packages** | slingshot, tradeSeq, SingleCellExperiment, viridis |

---

## Dependencies

```r
# CRAN / Bioconductor
install.packages(c("Seurat", "Signac", "tidyverse", "cowplot",
                   "patchwork", "ggrepel", "igraph", "ggraph",
                   "scales", "viridis", "Matrix"))

BiocManager::install(c("hdWGCNA", "WGCNA", "chromVAR",
                       "JASPAR2022", "motifmatchr",
                       "ChIPseeker", "DESeq2",
                       "slingshot", "tradeSeq",
                       "SingleCellExperiment",
                       "GenomicRanges", "IRanges",
                       "BiocParallel"))

# From GitHub
devtools::install_github("smorabit/hdWGCNA")
devtools::install_github("cole-trapnell-lab/cicero-release",
                         ref = "monocle3")

Citation
Subba P et al. (2026). Multiome Profiling Reveals Astrocyte and
Neuroendocrine Targets of Prenatal Acoustic Programming in Zebra Finch
Embryos. [Journal TBD].
