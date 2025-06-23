#### 1. Data Processing and Quality Control (QC) ####
# RNA (snRNA-seq) 

library(Seurat)
library(SeuratObject)
library(stringr)

# Load mitochondrial gene list
mt_genes <- read.delim("/data2/george_lab/psubba/zebra_finch_reference_data/mt_chromosome.gtf", header = FALSE, comment.char = "#")
mt_gene_ids <- unique(str_extract(mt_genes$V9, "(?<=gene_id\\s)[^\";]+")) %>% na.omit() %>% as.character()

# Create Seurat objects for each sample
sample_ids <- paste0("S", 1:8)
seurat_objects <- list()
for (id in sample_ids) {
  data_dir <- paste0("/data2/george_lab/psubba/", id, "_combined_analysis/outs")
  data_raw <- Read10X_h5(paste0(data_dir, "/filtered_feature_bc_matrix.h5"))
  if (is.list(data_raw) && "Gene Expression" %in% names(data_raw)) {
    data <- data_raw$`Gene Expression`
  } else {
    data <- data_raw
  }
  seurat_obj <- CreateSeuratObject(counts = data, project = id, min.cells = 3, min.features = 200)
  mt_gene_ids_seurat <- mt_gene_ids %>% str_replace_all("_", "-")
  seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), mt_gene_ids_seurat))
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_objects[[id]] <- seurat_obj
}

# Add metadata (Condition, Sex)
metadata <- data.frame(
  Sample = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"),
  Condition = c("Heat Call", "Heat Call", "Control", "Control", "Heat Call", "Heat Call", "Control", "Control"),
  Sex = c("Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female")
)
for (i in seq_along(seurat_objects)) {
  seurat_objects[[i]]$Condition <- metadata$Condition[i]
  seurat_objects[[i]]$Sex <- metadata$Sex[i]
}

# Merge all samples
merged_seurat <- merge(seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = sample_ids, project = "Merged")

# QC Filtering
merged_seurat_filt <- subset(merged_seurat, nFeature_RNA > 200 & nFeature_RNA < 3500 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mt < 5)
