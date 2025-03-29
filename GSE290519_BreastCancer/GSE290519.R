# setwd("F:/Practice work/R/GSE290519_BreastCancer")

# package installation

# Install CRAN packages
install.packages(c("dplyr", "readr", "tidyverse", "tidyr", "ggplot2", "pheatmap"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "DESeq2"))

# load library
library(GEOquery)
library(dplyr)
library(readr)
library(tidyverse)
library(tidyr)
library(DESeq2)
library(ggplot2)
library(pheatmap)

# metadata createing
gse <- getGEO(GEO = "GSE290519", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
metadata.modified <- metadata %>%
  select(23,2,47,49) %>%
  setNames(c("celltype", colnames(metadata)[2], "duration(day)", "treatment"))

#load dataset

dat <- read_table("F:/Practice work/R/GSE290519_BreastCancer/GSE290519_Raw_counts_all.txt")
dat.long <- dat %>%
  gather(key = "sample", value = "FPKM", -Gene)

# join dataframes = dat.long + metadata.modified

dat.long <- dat.long %>%
 left_join(., metadata.modified, by = c("sample" = "celltype"))

# dataset set as matrix for DESeq perform
dat.wide <- dat %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# leveling the metadata
metadata.modified$treatment = factor(metadata.modified$treatment, levels = c("Vehicle", "DEX", "Untreated"))


# Ensure sample names match between count data and metadata
metadata.modified <- metadata.modified %>%
  filter(celltype %in% colnames(dat.wide)) %>%
  arrange(match(celltype, colnames(dat.wide)))

rownames(metadata.modified) <- metadata.modified$celltype
metadata.modified$treatment[is.na(metadata.modified$treatment)] <- "DEX" # missing DEX force entry

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = dat.wide,
  colData = metadata.modified,
  design = ~ treatment
)

# Run DESeq2 pipeline
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("treatment", "DEX", "Vehicle"), alpha = 1e-5) # Compare DEX vs. Vehicle
summary(res)

# ========================================
# Visualization
# ========================================

# MA plot
plotMA(res)

# Volcano plot
res_df <- as.data.frame(res)
res_df$significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")
res_sig <- res[which(res$padj < 0.05), ]
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  scale_color_manual(values = c("red", "darkgreen")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(Adjusted p-value)")

# heatmap

# Top 50 genes based on both log2FC & significance
top_genes <- rownames(res_sig[order(res_sig$padj, res_sig$log2FoldChange, decreasing = c(FALSE, TRUE))[1:50], ])


# Apply variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) # Set blind = FALSE if you want to account for the experimental design
vsd_top_genes <- assay(vsd)[top_genes, ]

rownames(metadata.modified) <- metadata.modified$celltype
all(colnames(vsd_top_genes) %in% rownames(metadata.modified))
annotation_col <- metadata.modified[colnames(vsd_top_genes), , drop = FALSE]

# Create heatmap
pheatmap(vsd_top_genes, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row", 
         annotation_col = annotation_col[, "treatment", drop = FALSE], 
         main = "Heatmap of Top 50 Differentially Expressed Genes based on both log2FC & significance")


# PCA Analysis for Sample Clustering

vsd <- vst(dds, blind = TRUE) # Variance stabilizing transformation
metadata.modified$treatment <- as.factor(metadata.modified$treatment)
plotPCA(vsd, intgroup = "treatment")

write.csv(as.data.frame(res), "DESeq2_results.csv")
write.csv(as.data.frame(res_sig), "DESeq2_significant_results.csv")






