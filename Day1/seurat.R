rm(list=ls())
library(Seurat)
pbmc.data <- Seurat::Read10X(data.dir = "/Users/kevinlin/Downloads/filtered_gene_bc_matrices/hg19")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc