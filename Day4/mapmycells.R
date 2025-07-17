rm(list=ls())

library(Seurat)
library(DESeq2)

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/seaad_microglia.RData")

mapmycell_results <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/mapmycells_results/ec43c19b-1693-42c8-9200-423d649aa8cf_10xWholeHumanBrain(CCN202210140)_HierarchicalMapping_UTC_1752762564994.csv",
                              skip = 4)
rownames(mapmycell_results) <- mapmycell_results[,"cell_id"]
mapmycell_results <- mapmycell_results[,-1]
colnames(mapmycell_results) <- paste0("mapmycells_", colnames(mapmycell_results))
mapmycell_results <- mapmycell_results[Seurat::Cells(seurat_obj),]

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, mapmycell_results)

# Visualizing the predicted labels
scCustomize::DimPlot_scCustom(seurat_obj,
                              reduction = "umap",
                              group.by = "mapmycells_supercluster_name")

p1 <- scCustomize::FeaturePlot_scCustom(seurat_obj,
                                        reduction = "umap",
                                        features = "mapmycells_supercluster_bootstrapping_probability")
p2 <- Seurat::VlnPlot(seurat_obj,
                      features = "mapmycells_supercluster_bootstrapping_probability")
p1 + p2
