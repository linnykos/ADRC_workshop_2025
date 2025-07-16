rm(list=ls())
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(scCustomize)

SeuratDisk::Convert("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/adata_simplified.h5ad", 
                    dest = "h5seurat", 
                    overwrite = TRUE)

# Load the converted file as a Seurat object
seurat_obj <- SeuratDisk::LoadH5Seurat("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/adata_simplified.h5seurat")
seurat_obj

# Now put the raw counts back in
raw_mat <- Seurat::ReadMtx(
  mtx = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/adata_raw_10x/matrix.mtx",
  features = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/adata_raw_10x/features.tsv",
  cells = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/adata_raw_10x/barcodes.tsv"
)

# Make sure the raw counts has the same features and barcodes
raw_mat <- raw_mat[SeuratObject::Features(seurat_obj),
                   Seurat::Cells(seurat_obj)]

SeuratObject::LayerData(seurat_obj,
                        layer = "counts",
                        assay = "RNA") <- raw_mat

# Put in the metadata
metadata <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/adata_obs.csv")
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, metadata)

# Since there's already an scVI and UMAP, we can directly visualize the data
scCustomize::DimPlot_scCustom(seurat_obj,
                              reduction = "umap",
                              group.by = "Supertype")

scCustomize::DimPlot_scCustom(seurat_obj,
                              reduction = "umap",
                              group.by = "ADNC")

