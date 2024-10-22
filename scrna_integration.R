library(Seurat)


# Set the base directory
base_dir <- "/home/ranya/data/out"

# Identify all sample directories
sample_dirs <- list.dirs(path = base_dir, full.names = TRUE, recursive = FALSE)

# Initialize a list to store Seurat objects
seurat_objects <- list()

# Process each directory
for (dir in sample_dirs) {
  # Define the path to the filtered feature-barcode matrix .h5 file
  h5_file <- file.path(dir, "outs/filtered_feature_bc_matrix", "filtered_feature_bc_matrix.h5")
  
  # Check if the .h5 file exists
  if (file.exists(h5_file)) {
    # load the data matrix
    message(paste("Loading data from:", h5_file))
    seurat_obj <- Read10X(data.dir = dirname(h5_file))
    
    # Store the Seurat object in the list
    seurat_objects[[basename(dir)]] <- seurat_obj
  } else {
    message(paste("No .h5 data file found in:", dir))
  }
}

# Create seurat objects normalize and identify variable features 
seurat_RSD <- lapply(seurat_objects, CreateSeuratObject)
seurat_RSD <- lapply( seurat_RSD, NormalizeData)
seurat_RSD <- lapply( seurat_RSD,FindVariableFeatures)
# Select features that are frequently variable across all objects
features <- SelectIntegrationFeatures(seurat_RSD)
# Perform integration
chd.anchors<-FindIntegrationAnchors(object.list = seurat_RSD, anchor.features = features)
chd.combined <- IntegrateData(anchorset = chd.anchors)
chd.combined <- ScaleData(chd.combined,verbose = FALSE)
chd.combined <-RunPCA(chd.combined, npcs = 30, verbose = FALSE)
chd.combined <- RunUMAP(chd.combined, reduction = "pca", dims = 1:30)
# Computing nearest neighbor graph
# Computing SNN
chd.combined <- FindNeighbors(chd.combined, reduction = "pca", dims = 1:30)
# Return a seurat object where idents have been updated with new cluster info, 
# Change resolution value to obtain larger or smaller number of clusters
chd.combined  <- FindClusters(chd.combined , resolution = 0.3)

#To display all clusters
levels(chd.combined@meta.data$seurat_clusters)

