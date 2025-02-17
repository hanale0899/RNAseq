##This code is to prepare merged RDS files representing for each biological group for shinyApp## #nolint
##Version 09.12.2024 by Hannah Le##
##Have fun making science works!##

#Load libraries##nolint
library(hdf5r)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggplot2)



#miRs samples
data_dir1 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Tier2_scRNA/B1_counts_matrix" # nolint
data_dir2 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Tier2_scRNA/B3_counts_matrix"# nolint
data_dir3 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Tier2_scRNA/S50_counts_matrix" #nolint

#Neg_miRs samples
data_dir4 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Tier2_scRNA/B5_counts_matrix" # nolint
data_dir5 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Tier2_scRNA/B7_counts_matrix" # nolint
data_dir6 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Tier2_scRNA/S45_counts_matrix" # nolint


#Sham samples
data_dir7 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Tier2_scRNA/B10_counts_matrix" # nolint
data_dir8 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Tier2_scRNA/B11_counts_matrix"# nolint


list.files(data_dir7) # Should show barcodes.tsv, genes.tsv, and matrix.mtx

#Read files using Read10X from Seurat# #nolint
data1 <- Read10X(data.dir = data_dir1)
data2 <- Read10X(data.dir = data_dir2)

data3 <- Read10X(data.dir = data_dir3)
data4 <- Read10X(data.dir = data_dir4)
data5 <- Read10X(data.dir = data_dir5)

data6 <- Read10X(data.dir = data_dir6)
data7 <- Read10X(data.dir = data_dir7)
data8 <- Read10X(data.dir = data_dir8)

#Creat Seurat Object for the analysis# #nolint
seurat_object1 = CreateSeuratObject(counts = data1)
seurat_object2 = CreateSeuratObject(counts = data2)

seurat_object3 = CreateSeuratObject(counts = data3)
seurat_object4 = CreateSeuratObject(counts = data4)
seurat_object5 = CreateSeuratObject(counts = data5)

seurat_object6 = CreateSeuratObject(counts = data6)
seurat_object7 = CreateSeuratObject(counts = data7)
seurat_object8 = CreateSeuratObject(counts = data8)

##QC for percentage of mitochondrial genes
seurat_object1[["percent.mt"]] <- PercentageFeatureSet(seurat_object1, pattern = "^MT-") # nolint

# Filter out low-quality cells
seurat_obj_filtered1 <- subset(seurat_object1, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint
                       percent.mt < 5)                # Less than 10% mitochondrial genes for human genome # nolint


seurat_object2[["percent.mt"]] <- PercentageFeatureSet(seurat_object2, pattern = "^MT-") # nolint

# Filter out low-quality cells
seurat_obj_filtered2 <- subset(seurat_object2, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint
                       percent.mt < 5)                # Less than 10% mitochondrial genes for human genome # nolint

seurat_object3[["percent.mt"]] <- PercentageFeatureSet(seurat_object3, pattern = "^MT-") # nolint: line_length_linter.

# Filter out low-quality cells
seurat_obj_filtered3 <- subset(seurat_object3, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint: line_length_linter.
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint: line_length_linter.
                       percent.mt < 5)                # Less than 10% mitochondrial genes for human genome # nolint: line_length_linter.

seurat_object4[["percent.mt"]] <- PercentageFeatureSet(seurat_object4, pattern = "^MT-") # nolint: line_length_linter.

# Filter out low-quality cells
seurat_obj_filtered4 <- subset(seurat_object4, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint: line_length_linter.
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint: line_length_linter.
                       percent.mt < 5)                # Less than 10% mitochondrial genes for human genome # nolint: line_length_linter.

seurat_object5[["percent.mt"]] <- PercentageFeatureSet(seurat_object5, pattern = "^MT-") # nolint: line_length_linter.

# Filter out low-quality cells
seurat_obj_filtered5 <- subset(seurat_object5, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint: line_length_linter.
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint: line_length_linter.
                       percent.mt < 5)                # Less than 10% mitochondrial genes for human genome # nolint: line_length_linter.

seurat_object6[["percent.mt"]] <- PercentageFeatureSet(seurat_object6, pattern = "^MT-") # nolint: line_length_linter.
# Filter out low-quality cells
seurat_obj_filtered6 <- subset(seurat_object6, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint: line_length_linter.
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint: line_length_linter.
                       percent.mt < 5)                # Less than 10% mitochondrial genes for human genome # nolint: line_length_linter.

seurat_object7[["percent.mt"]] <- PercentageFeatureSet(seurat_object7, pattern = "^MT-") # nolint: line_length_linter.
# Filter out low-quality cells
seurat_obj_filtered7 <- subset(seurat_object7, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint: line_length_linter.
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint: line_length_linter.
                       percent.mt < 5)                # Less than 10% mitochondrial genes for human genome # nolint: line_length_linter.

seurat_object8[["percent.mt"]] <- PercentageFeatureSet(seurat_object8, pattern = "^MT-") # nolint: line_length_linter.
# Filter out low-quality cells
seurat_obj_filtered8 <- subset(seurat_object8, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint: line_length_linter.
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint: line_length_linter.
                       percent.mt < 5)                # Less than 10% mitochondrial genes for human genome # nolint: line_length_linter.
# Normalize the data
seurat_obj_norm1 <- NormalizeData(seurat_obj_filtered1, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.
seurat_obj_norm2 <- NormalizeData(seurat_obj_filtered2, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.
seurat_obj_norm3 <- NormalizeData(seurat_obj_filtered3, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.
seurat_obj_norm4 <- NormalizeData(seurat_obj_filtered4, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.
seurat_obj_norm5 <- NormalizeData(seurat_obj_filtered5, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.
seurat_obj_norm6 <- NormalizeData(seurat_obj_filtered6, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.
seurat_obj_norm7 <- NormalizeData(seurat_obj_filtered7, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.
seurat_obj_norm8 <- NormalizeData(seurat_obj_filtered8, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.

#export raw counts to csv
# Assuming seurat_obj is your Seurat object
# Extract the count matrix from the Seurat object
count_matrix1 <- GetAssayData(seurat_obj_norm1, layer = "counts")
count_matrix2 <- GetAssayData(seurat_obj_norm2, layer = "counts")
count_matrix3 <- GetAssayData(seurat_obj_norm3, layer = "counts")
count_matrix4 <- GetAssayData(seurat_obj_norm4, layer = "counts")
count_matrix5 <- GetAssayData(seurat_obj_norm5, layer = "counts")
count_matrix6 <- GetAssayData(seurat_obj_norm6, layer = "counts")
count_matrix7 <- GetAssayData(seurat_obj_norm7, layer = "counts")
count_matrix8 <- GetAssayData(seurat_obj_norm8, layer = "counts")
# Convert to a data frame
count_df1 <- as.data.frame(as.matrix(count_matrix1))
count_df2 <- as.data.frame(as.matrix(count_matrix2))
count_df3 <- as.data.frame(as.matrix(count_matrix3))
count_df4 <- as.data.frame(as.matrix(count_matrix4))
count_df5 <- as.data.frame(as.matrix(count_matrix5))
count_df6 <- as.data.frame(as.matrix(count_matrix6))
count_df7 <- as.data.frame(as.matrix(count_matrix7))
count_df8 <- as.data.frame(as.matrix(count_matrix8))
# Write to CSV (with gene names as row names)
write.csv(count_df1, file = "scRNA_count_matrix_samp5.csv", quote = FALSE)
write.csv(count_df2, file = "scRNA_count_matrix_samp8.csv", quote = FALSE)
write.csv(count_df3, file = "scRNA_count_matrix_samp2.csv", quote = FALSE)
write.csv(count_df4, file = "scRNA_count_matrix_samp7.csv", quote = FALSE)
write.csv(count_df5, file = "scRNA_count_matrix_samp9.csv", quote = FALSE)
write.csv(count_df6, file = "scRNA_count_matrix_samp1.csv", quote = FALSE)

# Write to rds (with gene names as row names)
saveRDS(count_matrix1, file = "scRNA_count_matrix_B1.rds")
saveRDS(count_matrix2, file = "scRNA_count_matrix_B3.rds")
saveRDS(count_matrix3, file = "scRNA_count_matrix_S50.rds")
saveRDS(count_matrix4, file = "scRNA_count_matrix_B5.rds")
saveRDS(count_matrix5, file = "scRNA_count_matrix_B7.rds")
saveRDS(count_matrix6, file = "scRNA_count_matrix_S45.rds")
saveRDS(count_matrix7, file = "scRNA_count_matrix_B10.rds")
saveRDS(count_matrix8, file = "scRNA_count_matrix_B11.rds")
# (Optional) Ensure common genes between all objects
common_genes_mir <- Reduce(intersect, list(rownames(seurat_obj_norm1), rownames(seurat_obj_norm2), rownames(seurat_obj_norm3)))

# Subset to common genes
seurat_obj_norm1 <- subset(seurat_obj_norm1, features = common_genes_mir)
seurat_obj_norm2<- subset(seurat_obj_norm2, features = common_genes_mir)
seurat_obj_norm3<- subset(seurat_obj_norm3, features = common_genes_mir)
# Merge the multiple Seurat objects into one
#4_weeks_samples #nolint
combined_seurat_mir <- merge(seurat_obj_norm1, y = c(seurat_obj_norm2,seurat_obj_norm3), 
                         add.cell.ids = c("Sample1", "Sample2","Sample3"))

# Normalize the combined object (optional)
combined_seurat_mir_norm <- NormalizeData(combined_seurat_mir)

# Scale and run PCA (optional)
combined_seurat_mir <- FindVariableFeatures(combined_seurat_mir)
combined_seurat_mir <- ScaleData(combined_seurat_mir, features = VariableFeatures(combined_seurat_mir)) #nolint
combined_seurat_mir <- RunPCA(combined_seurat_mir, features = VariableFeatures(combined_seurat_mir)) #nolint

MIR <- GetAssayData(combined_seurat_mir_norm, layer = "counts")
write.csv(combined_seurat_mir, file = "RNAexpressionMirs.csv")
# Perform clustering (optional)
combined_seurat_mir <- FindNeighbors(combined_seurat_mir, dims = 1:20)
combined_seurat_mir <- FindClusters(combined_seurat_mir, resolution = 1.1)

combined_seurat_mir_umap <- RunUMAP(combined_seurat_mir, dims = 1:20)
a <- DimPlot(combined_seurat_mir_umap , reduction = "umap")
ggsave(a, file="rawclustermir.jpg", width=6, height=6)
# Save the merged Seurat object
saveRDS(combined_seurat_mir_norm, file = "combined_mir.rds")

#8_weeks_samples #nolint
combined_seurat_neg <- merge(seurat_obj_norm4, y = c(seurat_obj_norm5, seurat_obj_norm6), 
                         add.cell.ids = c("Sample1", "Sample2", "Sample3"))

# Normalize the combined object (optional)
combined_seurat_neg <- NormalizeData(combined_seurat_neg)
combined_seurat_neg_norm <- NormalizeData(combined_seurat_neg)

# Scale and run PCA (optional)
combined_seurat_neg <- FindVariableFeatures(combined_seurat_neg)
combined_seurat_neg <- ScaleData(combined_seurat_neg, features = VariableFeatures(combined_seurat_neg)) #nolint
combined_seurat_neg <- RunPCA(combined_seurat_neg, features = VariableFeatures(combined_seurat_neg)) #nolint

# Perform clustering (optional)
combined_seurat_neg <- FindNeighbors(combined_seurat_neg, dims = 1:20)
combined_seurat_neg <- FindClusters(combined_seurat_neg, resolution = 1.1)

combined_seurat_neg_umap <- RunUMAP(combined_seurat_neg, dims = 1:20)
b <- DimPlot(combined_seurat_neg_umap , reduction = "umap")
ggsave(b, file="rawclusterneg.jpg", width=6, height=6)
# Save the merged Seurat object
saveRDS(combined_seurat_neg_norm, file = "combined_neg.rds")

#12_weeks_samples #nolint
combined_seurat_sham <- merge(seurat_obj_norm7, y = seurat_obj_norm8, 
                         add.cell.ids = c("Sample1", "Sample2"))

# Normalize the combined object (optional)
combined_seurat_sham <- NormalizeData(combined_seurat_sham)
combined_seurat_sham_norm <- NormalizeData(combined_seurat_sham)

# Scale and run PCA (optional)
combined_seurat_sham <- FindVariableFeatures(combined_seurat_sham)
combined_seurat_sham <- ScaleData(combined_seurat_sham, features = VariableFeatures(combined_seurat_sham)) #nolint
combined_seurat_sham <- RunPCA(combined_seurat_sham, features = VariableFeatures(combined_seurat_sham)) #nolint

# Perform clustering (optional)
combined_seurat_sham <- FindNeighbors(combined_seurat_sham, dims = 1:20)
combined_seurat_sham <- FindClusters(combined_seurat_sham, resolution = 1.1)

combined_seurat_sham_umap <- RunUMAP(combined_seurat_sham, dims = 1:20)
c <- DimPlot(combined_seurat_sham_umap , reduction = "umap")
ggsave(c, file="rawclustersham.jpg", width=6, height=6)
# Save the merged Seurat object
saveRDS(combined_seurat_sham_norm, file = "combined_sham.rds")



