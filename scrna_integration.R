library(Seurat)

# Set the base directory
SRR19266863_dir <- "/home/ranya/data/out/SRR19266863/outs/filtered_feature_bc_matrix"

# Read file
SRR19266864_obj <- Read10X(data.dir = SRR19266864_dir)

#Create RSD
obj1 <-CreateSeuratObject(counts=SRR19266863_obj)

#Normalize data
