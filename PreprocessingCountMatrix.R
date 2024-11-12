#BiocManager::install('Seurat')
#install.packages("hdf5r")
library(hdf5r)
library(Seurat)
library(SeuratObject)
data_dir <- "/Users/littlevanilla/Downloads/Sample_9/Sample_9_counts_matrix"
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
#if (FALSE) {
  # For output from CellRanger < 3.0
  #data_dir <- 'path/to/data/directory'
  #list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
  #expression_matrix <- Read10X(data.dir = data_dir)
  #seurat_object = CreateSeuratObject(counts = expression_matrix)
  
  # For output from CellRanger >= 3.0 with multiple data types
  #data_dir <- 'path/to/data/directory'
  #list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
  data <- Read10X(data.dir = data_dir)
  seurat_object = CreateSeuratObject(counts = data)
  #seurat_object[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)
#}
  
BiocManager::install('Matrix')
##QC for percentage of mitochondrial genes
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

# Filter out low-quality cells
seurat_obj_filtered <- subset(seurat_object, 
                     subset = nFeature_RNA > 200 &  # At least 200 genes detected
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets)
                       percent.mt < 10)                # Less than 10% mitochondrial genes for human genome

# Normalize the data
seurat_obj_norm <- NormalizeData(seurat_obj_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify the top 10000 variable genes
seurat_obj <- FindVariableFeatures(seurat_obj_norm, selection.method = "vst", nfeatures = 10000)

# Scale the data (standardize the expression of each gene across all cells)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

# Run PCA on variable features
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj_norm))

# Visualize the PCA results
DimPlot(seurat_obj, reduction = "pca")

# Find neighbors and clusters (using the first 20 principal components)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Visualize clusters using UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
DimPlot(seurat_obj, reduction = "umap")

# Identify differentially expressed genes between clusters
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View the top markers
head(markers)

# Save Seurat object to disk
saveRDS(seurat_obj_norm, file = "processed_seurat_object.rds")

