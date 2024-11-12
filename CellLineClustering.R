# Visualize the expression of neural cell markers in UMAP
FeaturePlot(seurat_obj, features = c("SOX2", "S100B", "OLIG2", "AIF1","TUBB3", "MAP2", "SMIM32", "FOXP1", "CALB1","NTRK2","RUNX1","UCHL1"))

# Identify variable features (genes) to focus on
seurat_obj_norm_1 <- FindVariableFeatures(seurat_obj_norm, selection.method = "vst", nfeatures = 10000)
# Identify differentially expressed genes for each cluster
marker_genes <- c("SOX2", "S100B", "OLIG2", "AIF1","TUBB3", "MAP2", "SMIM32", "FOXP1", "CALB1","NTRK2","RUNX1","UCHL1")

seurat_obj_markers <- subset(seurat_obj, features = marker_genes)

# Run PCA on the subset of marker genes
seurat_obj_markers <- RunPCA(seurat_obj_markers, features = marker_genes)

# Visualize the PCA results
DimPlot(seurat_obj_markers, reduction = "pca")

# Find neighbors (based on PCA dimensions, usually the first 10–20 PCs)
seurat_obj_markers <- FindNeighbors(seurat_obj_markers, dims = 1:10)

# Run clustering (using a resolution that suits your data, usually between 0.5 and 1.5)
seurat_obj_markers <- FindClusters(seurat_obj_markers, resolution = 1.1)

# Visualize the clusters in PCA space
DimPlot(seurat_obj_markers, reduction = "pca", group.by = "seurat_clusters")

# Visualize marker gene expression across clusters
FeaturePlot(seurat_obj_markers, features = marker_genes)

# Or, use violin plots to see marker gene expression in each cluster
VlnPlot(seurat_obj_markers, features = marker_genes, group.by = "seurat_clusters")

# View the current cluster labels
table(seurat_obj_markers$seurat_clusters)

# Manually change the cluster labels (for example, rename cluster 0 to 'Neurons', cluster 1 to 'Glial Cells', etc.)
seurat_obj_markers$seurat_clusters <- as.character(seurat_obj_markers$seurat_clusters)

##Based on the highlighted marker of the cluster, assign the cell line annotation
seurat_obj_markers$seurat_clusters[seurat_obj_markers$seurat_clusters == "0"] <- "Unknown"
seurat_obj_markers$seurat_clusters[seurat_obj_markers$seurat_clusters == "1"] <- "Maturing astrocytes/oligodendrocytes"
seurat_obj_markers$seurat_clusters[seurat_obj_markers$seurat_clusters == "2"] <- "Microglia"
seurat_obj_markers$seurat_clusters[seurat_obj_markers$seurat_clusters == "3"] <- "Microglia"
seurat_obj_markers$seurat_clusters[seurat_obj_markers$seurat_clusters == "4"] <- "Microglia/Sensory neurons"
seurat_obj_markers$seurat_clusters[seurat_obj_markers$seurat_clusters == "5"] <- "Motor neurons"
seurat_obj_markers$seurat_clusters[seurat_obj_markers$seurat_clusters == "6"] <- "Microglia"

library(ggplot2)

# Generate the DimPlot object
p <- DimPlot(seurat_obj_markers, reduction = "umap", group.by = "seurat_clusters")

# Add a custom title using ggtitle()
p <- p + ggtitle("Cell type clustering")

ggsave(p, file="cluster9.jpg", width=8, height=6)
