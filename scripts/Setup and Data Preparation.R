# Set working directory (change as appropriate)
setwd("your/working/directory")

# Load required packages
library(WGCNA)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(pheatmap)
library(DESeq2)  # For VST normalization

# Avoid factor conversion issues
options(stringsAsFactors = FALSE)

# Increase memory limit for large datasets
options(future.globals.maxSize = 250000 * 1024^2)

# Step 1: Prepare gene length data from GFF annotation
my_tags <- c("Name", "Dbxref", "gene")
my_columns <- c("seqid", "start", "end", "strand", "type")
my_filter <- list(type = "gene")

dat <- readGFF("GCF_003957565.2_bTaeGut1.4.pri_genomic.gff", 
               tags = my_tags, 
               columns = my_columns, 
               filter = my_filter)

dat$gene_name <- dat$gene
dat$interval_start <- dat$start
dat$interval_stop <- dat$end

dat_1 <- as.data.frame(dat)
all_data_from_dat <- dplyr::select(dat_1, gene_name, start, end)

# Calculate gene length as end - start - 1
all_data_from_dat$geneLength <- (all_data_from_dat$end - all_data_from_dat$start) - 1

gene_name_length <- dplyr::select(all_data_from_dat, gene_name, geneLength)

# Step 2: Read and process raw gene counts data
data0_csv <- read.csv("WGCNA_gene_count_7_June_2023.csv", header = TRUE)

# Rename first column for merging
colnames(data0_csv)[1] <- "gene_name"

# Merge gene length data with counts (geneLength retained for reference)
data0 <- merge(data0_csv, gene_name_length, by = "gene_name", all.x = TRUE)

# Set rownames as gene names and remove gene_name column
colnames(data0)[1] <- ""
rownames(data0) <- data0[, 1]
data0 <- data0[, -1]

# Convert numeric columns to integer counts (raw counts)
data0[] <- lapply(data0, function(x) if(is.numeric(x)) as.integer(x) else x)

# Step 3: Load sample metadata
sample_metadata <- read.csv(file = "WGCNA_Heat_Call_Experimental_Design.csv")
colnames(sample_metadata) <- c('sample_ID', 'Sex', 'Condition')
rownames(sample_metadata) <- sample_metadata$sample_ID
