# https://satijalab.org/seurat/articles/pbmc3k_tutorial 

rm(list=ls())
library(Seurat)
pbmc.data <- Seurat::Read10X(data.dir = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/data/seurat-pbmc/filtered_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k")
pbmc

# see https://satijalab.org/seurat/articles/pbmc3k_tutorial#identification-of-highly-variable-features-feature-selection
pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# plot variable features with and without labels
plot1 <- Seurat::VariableFeaturePlot(pbmc)
top10 <- head(Seurat::VariableFeatures(pbmc), 10)
plot2 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
