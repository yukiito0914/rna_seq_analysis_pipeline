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

# Extracr signuficant genes
sig_genes <- res %>%
  filter(padj < padj_threshold) %>%
  select(gene_name) %>%
  distinct()

# Write significant genes to file for Enrichr analysis
write.table(sig_genes, "enrichr_input.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

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

# Vocano Plot
# Remove NA gene names
res <- res %>% filter(!is.na(gene_name))

# Find the maximum value of p-value corresponding to FDR (padj) < 0.01
pval_cutoff <- max(res$pvalue[res$padj < 0.01], na.rm = TRUE)
log_pval_cutoff <- -log10(pval_cutoff)

# Classify DEGs
res$diffexpressed <- "NO"
res$diffexpressed[res$padj < 0.01 & res$log2FoldChange > 0] <- "UP"
res$diffexpressed[res$padj < 0.01 & res$log2FoldChange < 0] <- "DOWN"

# Set color
custom_colors <- c("UP" = "#FF6F61",   # Soft Red
                   "DOWN" = "#5A5EBC", # Iris Blue
                   "NO" = "gray50")

# Set y-axis limit (clip values above 16)
res$log_pvalue <- pmin(-log10(res$pvalue), 16)

# Create volcano plot
volcano_plot <- ggplot(res, aes(x=log2FoldChange, y=log_pvalue, color=diffexpressed)) +
    geom_point(size = 1) + 
    theme_minimal() +
    scale_color_manual(values=custom_colors) + 
    geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype="dashed") +  
    labs(x="Expression difference (log2 FC)",
         y="Significance (-log10 p value)") +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black")  
    )

# Count upregulated and downregulated genes
num_upregulated <- sum(res$diffexpressed == "UP", na.rm = TRUE)
num_downregulated <- sum(res$diffexpressed == "DOWN", na.rm = TRUE)

# Create summary text
regulation_summary <- paste0(
  "Number of upregulated genes (log2FC > 0, padj < 0.01):", num_upregulated, "\n",
  "Number of downregulated genes (log2FC < 0, padj < 0.01):", num_downregulated, "\n"
)

# Save summary to file
writeLines(regulation_summary, "regulation_summary.txt")

# Save volcano plot
ggsave(filename = "fig3.c.png", plot = volcano_plot, width = 4, height = 4, dpi = 300, bg = "white")

# GSEA plot
# Calculate the percentage of DE genes in each pathway
fgsea_results <- fgsea_results %>%
  mutate(
    num_DE_genes = sapply(leadingEdge, length),  # Number of DE genes in each pathway
    num_total_genes = size,  # Total number of genes in the pathway
    percentage_DE_genes = (num_DE_genes / num_total_genes) * 100  # Calculate percentage
  )

# Select the top 5 upregulated pathways
top_up <- fgsea_results %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  head(5)

# Select the top 5 downregulated pathways
top_down <- fgsea_results %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  head(5)

# Add a new column to label Up/Down pathways
top_up$regulation <- "Upregulated"
top_down$regulation <- "Downregulated"

# Combine upregulated and downregulated pathways
top_pathways <- bind_rows(top_up, top_down)

# Ensure correct Y-axis order within each regulation group based on NES
top_pathways <- top_pathways %>%
  mutate(regulation = factor(regulation, levels = c("Upregulated", "Downregulated"))) %>% # Ensure order of facets
  group_by(regulation) %>%
  mutate(pathway = reorder(pathway, NES)) %>%  # Order pathways by NES
  ungroup()

# Create bar plot with vertically aligned Upregulated and Downregulated pathways
gsea_bar <- ggplot(top_pathways, aes(x = percentage_DE_genes, y = pathway, fill = padj)) +
  geom_bar(stat = "identity") +  # Create bar plot
  theme_minimal() +  
  scale_fill_gradient(low = "blue", high = "red", name = "Adjusted p-value", trans = "reverse") +  # Adjust color scale based on padj
  labs(
    title = "GSEA Pathways with Percentage of DE Genes",
    x = "Percentage of DE genes in category / all genes in category (%)",
    y = "Pathway"
  ) +
  facet_grid(regulation ~ ., scales = "free_y") +  # Align Upregulated and Downregulated vertically
  xlim(0, max(top_pathways$percentage_DE_genes, na.rm = TRUE)) +  # Ensure the same X-axis scale for both facets
  theme(legend.position = "right")

  # Save GSEA bar plot
ggsave(filename = "fig3.f.png", plot = gsea_bar, width = 8, height = 6, dpi = 300, bg = "white")