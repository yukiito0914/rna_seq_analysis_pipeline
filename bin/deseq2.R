#!/usr/bin/env Rscript

# Load Libraries
library(DESeq2)
library(readr)
library(dplyr)
library(tibble)  
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(msigdbr)
library(fgsea)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
count_file <- args[1]
metadata_file <- args[2]
output_results <- args[3]
output_top10 <- args[4]
summary_file <- args[5]
padj_threshold <- as.numeric(args[6])
annotation_file <- args[7]  

# Load count data
counts <- read_tsv(count_file, col_names = TRUE)
counts <- as.data.frame(counts)

# Ensure 'gene' column exists
if (!"gene" %in% colnames(counts)) {
  stop("Error: The 'gene' column is missing from the count file. Check the format of filtered_counts.csv.")
}

# Set rownames and remove 'gene' column
rownames(counts) <- counts$gene
counts <- counts[, -1]  # Remove 'gene' column

# Convert columns to numeric while preserving rownames
counts[] <- lapply(counts, as.numeric)
counts <- as.matrix(counts)

# Check if counts is a valid matrix
if (!is.matrix(counts)) {
  stop("Error: counts is not a matrix. Conversion failed.")
}

# Load metadata
metadata <- read_csv(metadata_file)
metadata <- as.data.frame(metadata)

# Ensure 'sample' and 'condition' columns exist
if (!"sample" %in% colnames(metadata)) {
  stop("Error: 'sample' column is missing from metadata.csv")
}
if (!"condition" %in% colnames(metadata)) {
  stop("Error: 'condition' column is missing from metadata.csv")
}

# Set rownames and ensure correct format
rownames(metadata) <- as.character(metadata$sample)
metadata <- metadata[, -1, drop = FALSE]  # Keep metadata as a data.frame

# Check if metadata matches counts
if (!all(rownames(metadata) %in% colnames(counts))) {
  stop("Error: Row names of metadata do not match column names of counts.")
}

# Create DESeq2 dataset and run DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# Convert results to dataframe and add gene_id column (from rownames)
res <- as.data.frame(res)
res <- rownames_to_column(res, var = "gene_id")

# Load gene annotation mapping file (gene_id, gene_name)
gene_mapping <- read_tsv(annotation_file, col_names = TRUE)

# Join annotation mapping to DESeq2 results by gene_id
res <- left_join(res, gene_mapping, by = "gene_id")

# Save full results (now including gene_id and gene_name)
write_csv(res, output_results)

# Extract top 10 significant genes (by padj) and output with gene_id and gene_name
top10_genes <- res %>% arrange(padj) %>% head(10)
write_csv(top10_genes, output_top10)

# Count significant genes
num_significant_genes <- sum(res$padj < padj_threshold, na.rm = TRUE)

# Write summary to file
summary_text <- paste("Number of significant genes (padj <", padj_threshold, "):", num_significant_genes, "\n")
writeLines(summary_text, summary_file)

# Perform normalization (VST)
vsd <- vst(dds, blind=FALSE)

# Perform PCA
pca_plot <- plotPCA(vsd, intgroup=c("condition"))

# Save PCA plot
png("PCA_plot.png", width = 800, height = 600)
print(pca_plot)
dev.off()

# Create and save a heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png("Heatmap.png", width = 800, height = 600)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# Remove duplicates in gene_name
res <- res %>% distinct(gene_name, .keep_all = TRUE)

# Remove NAs
res <- res %>% filter(!is.na(gene_name))

# Create ranks for FGSEA Analysis
ranks <- res$log2FoldChange
names(ranks) <- res$gene_name 
ranks <- ranks[!is.na(ranks) & names(ranks) != ""]
ranks <- ranks[!duplicated(names(ranks))]  # Remove duplicates
ranks <- sort(ranks, decreasing = TRUE)

# Get C2:CP (Canonical Pathways)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP") %>%
    split(x = .$gene_symbol, f = .$gs_name)

# Perform fgsea
fgsea_results <- fgsea(pathways = m_t2g, stats = ranks, minSize = 15, maxSize = 500)

# Generate a plot that displays top 10 most significant pathways
topPathways <- fgsea_results[order(padj)][1:10, ]
p <- ggplot(topPathways, aes(reorder(pathway, NES), NES)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    labs(title = "Top 10 Enriched Pathways", x = "Pathway", y = "Normalized Enrichment Score (NES)")

# Save fgsea plot
ggsave("Top10_Pathways.png", plot = p, dpi = 300)