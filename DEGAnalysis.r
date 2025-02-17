#This code is written by Le Phuc Hoang Anh aka Hannah Le
#An intellectual property of ChewLab_NTU
#Last_modified on Feb_02_2025
#Input data required:*.rds
#System requirement: RAM > 64.0 Gb
#LoadLibraries
#If Libraries not installed, used BiocManager::install("Library_name") #nolint
library(hdf5r)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggplot2)
library(pheatmap)
##Part 1:DEG Analysis using FindMarkers() in Seurat v.5.0
# Load Seurat objects
seurat_object1 <- readRDS("combined_mir.rds")
seurat_object2 <- readRDS("combined_neg.rds")
seurat_object3 <- readRDS("combined_sham.rds")

# Add a "condition" metadata to each Seurat object
seurat_object1$condition <- "Mir"
seurat_object2$condition <- "NegMir"
seurat_object3$condition <- "Sham"

# Merge the Seurat objects
seurat_merged <- merge(seurat_object1, y = c(seurat_object2, seurat_object3), add.cell.ids = c("Mir", "NegMir", "Sham"))
seurat_merged_joined <- JoinLayers(seurat_merged)

# Assuming your Seurat object is named `seurat_obj`
# Extract the prefix from the cell ID (before the first underscore "_")
seurat_merged_joined$sample_ID <- sub("(^[^_]+_[^_]+).*", "\\1", colnames(seurat_merged_joined))

# View the first few entries to check the annotation
head(seurat_merged_joined$sample_ID)

# Set the identity to "conditionn" to compare between conditions
Idents(seurat_merged_joined) <- "condition"

# Perform differential expression analysis (DEG) between conditions
# Example: comparing Mir vs NegMir
DEG_results_mir_vs_negmir <- FindMarkers(seurat_merged_joined, ident.1 = "Mir", ident.2 = "NegMir")

# Example: comparing Mir vs Sham
DEG_results_mir_vs_sham <- FindMarkers(seurat_merged_joined, ident.1 = "Mir", ident.2 = "Sham")

# Example: comparing NegMir vs Sham
DEG_results_negmir_vs_sham <- FindMarkers(seurat_merged_joined, ident.1 = "NegMir", ident.2 = "Sham")

# View the results for Mir vs NegMir
head(DEG_results_mir_vs_negmir)

write.csv(DEG_results_mir_vs_negmir,"DEG_results_mir_vs_negmir.csv", quote = FALSE)
write.csv(DEG_results_mir_vs_sham,"DEG_results_mir_vs_sham.csv", quote = FALSE)
write.csv(DEG_results_negmir_vs_sham,"DEG_results_negmir_vs_sham.csv", quote = FALSE)

##Part 2:DATA Visualization
data1 <- read.csv("DEG_results_mir_vs_negmir.csv")
data2 <- read.csv("DEG_results_mir_vs_sham.csv")
data3 <- read.csv("DEG_results_negmir_vs_sham.csv")
data <- data1
# Calculate -log10(p_val_adj) for the plot
data$logP <- -log10(data$p_val_adj)

# Define significance threshold
significance_threshold <- 0.05
log2fc_threshold <- 1

data_sorted <- data[order(abs(data$avg_log2FC), data$p_val_adj), ]

# Get the top 5 genes based on the highest absolute log2FC and lowest p-values
top_genes <- head(data_sorted, 5)

# Create a new column for classification of significance
data$Significance <- "Not Very Significant"
data$Significance[data$p_val_adj < significance_threshold & data$avg_log2FC > log2fc_threshold] <- "Upregulated"
data$Significance[data$p_val_adj < significance_threshold & data$avg_log2FC < -log2fc_threshold] <- "Downregulated"

# Create Volcano Plot using ggplot2
Plot <- ggplot(data, aes(x=avg_log2FC, y=logP, color=Significance)) +
  geom_point(size=3) +
  scale_color_manual(values=c("blue", "gray", "red")) +
  theme_minimal() +
  labs(title="Volcano Plot", x="Average log2 Fold Change (avg_log2FC)", y="-log10(Adjusted p-value)") +
  theme(plot.title = element_text(hjust = 5)) +
  geom_text(data = top_genes, aes(label = rownames(top_genes)), vjust = -1, hjust = 0, size = 5, color = "black") +
  theme(legend.title = element_blank()) + ggtitle("MiR vs. NegMiR")

ggsave("Volcano_DEG_results_mir_vs_negmir.png", Plot, height = 8, width = 8)


# Extract the top DEGs from the results of DEG analysis
top_genes_mir_vs_negmir <- tail(data1[order(data1$p_val),], 30)
top_genes_mir_vs_sham <- tail(data2[order(data2$p_val),], 30)
top_genes_negmir_vs_sham <- tail(data3[order(data3$p_val),], 30)

top_genes_mir_vs_negmir2 <- tail(data1[order(data1$p_val),], 10)
top_genes_mir_vs_sham2 <- tail(data2[order(data2$p_val),], 10)
top_genes_negmir_vs_sham2 <- tail(data3[order(data3$p_val),], 10)

# Combine top DEGs across comparisons (optional)
top_genes <- unique(c(top_genes_mir_vs_negmir$X, top_genes_mir_vs_sham$X, top_genes_negmir_vs_sham$X, "Mbp", "Pdgfra"))
top_genes2 <- unique(c(top_genes_mir_vs_negmir2$X, top_genes_mir_vs_sham2$X, top_genes_negmir_vs_sham2$X, "Mbp", "Pdgfra"))

## Assuming avg_expression contains the aggregated data
avg_expression <- AggregateExpression(seurat_merged_joined, 
                                             group.by = "sample_ID",  # Group by sample (not condition)
                                             return.seurat = FALSE)  # Return as a matrix (not Seurat object)


# Extract the expression matrix of top DEGs from the Seurat object
expression_data <- avg_expression$RNA  # "RNA" is the default assay

# Now, extract the data for the top DEGs from the expression matrix
expression_data_top <- expression_data[top_genes, ]

# Remove rows or columns with all NA/NaN values (if applicable)
expression_data_clean <- na.omit(expression_data_top)

# Generate the heatmap of the top DEGs expression data
Plot2 <- pheatmap(expression_data_clean, 
         scale = "row",  # Scale each row (gene) to mean=0, sd=1
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         clustering_distance_rows = "none",  # Skip clustering for rows
         clustering_distance_cols = "none",  # Skip clustering for columns
         clustering_method = "none",
         show_rownames = TRUE,  # Show gene names
         show_colnames = TRUE,  # Show condition names
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Blue to red color scale
         main = "Gene Expression Heatmap Across Conditions",
         fontsize_row = 8,
         fontsize_col = 10)

ggsave("HeatmapExpression.png", Plot2, height = 12, width = 6)


# Generate the violin plot for the top 10 genes of each comparison
violin_plot <- VlnPlot(seurat_merged_joined, 
                       features = top_genes2, cols = c("#ff5874", "#34acdc", "black"),  # List of top 10 genes
                       group.by = "condition", # Group by condition (or sample)
                       pt.size = 0.1)        # This will allow custom ggplot adjustments

# Customize the plot using ggplot2 functions
Plot3 <- violin_plot+ 
  ggtitle("Top 10 DEGs Expression Across Conditions") +  # Add title
  ylab("Expression") +  # Set Y-axis label

ggsave("ViolinPlot_topgenes.png", Plot3, height = 10, width = 8)