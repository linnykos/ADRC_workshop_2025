# from https://satijalab.org/seurat/articles/integration_introduction
# https://satijalab.org/seurat/articles/seurat5_integration
rm(list=ls())
library(Seurat)
library(SeuratData)
library(patchwork)

options(future.globals.maxSize = 8000 * 1024^2) # https://github.com/satijalab/seurat/issues/1845
set.seed(10)

# install dataset
SeuratData::InstallData("ifnb")

# load dataset
ifnb <- SeuratData::LoadData("ifnb")
ifnb_original <- ifnb # save this for later, for pedagogical comparison
e
###############################

# Step 1: Split the RNA measurements into two layers one for control cells, one for stimulated cells
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb

# Step 2: "Process" each split
ifnb <- Seurat::SCTransform(ifnb)
ifnb <- Seurat::RunPCA(ifnb, npcs = 30, verbose = F)
ifnb <- Seurat::IntegrateLayers(
  object = ifnb,
  method = Seurat::RPCAIntegration,
  normalization.method = "SCT",
  verbose = F
)

ifnb <- Seurat::RunUMAP(ifnb, 
                        reduction = "integrated.dr", 
                        dims = 1:30, 
                        reduction.name = "umap.dr")
p1 <- Seurat::DimPlot(
  ifnb,
  reduction = "umap.dr",
  group.by = c("stim"),
  combine = FALSE, label.size = 2
)



##############################

#