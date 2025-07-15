# https://satijalab.org/seurat/articles/sctransform_vignette.html

rm(list=ls())
library(Seurat)
library(scCustomize)

pbmc.data <- Seurat::Read10X(data.dir = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/data/seurat-pbmc/filtered_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data, 
                                   project = "pbmc3k")
pbmc

# For this tutorial, let's show how you can include percent.mt and cell cycling scores 
## in the SCTransform
# Store mitochondrial percentage in object meta data
pbmc <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
# Cell-cycling. see more: https://satijalab.org/seurat/articles/cell_cycle_vignette.html
pbmc <- Seurat::NormalizeData(pbmc) #Seurat::CellCycleScoring requires the "data" slot filled
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pbmc <- Seurat::CellCycleScoring(pbmc, 
                                 s.features = s.genes, 
                                 g2m.features = g2m.genes)

head(pbmc@meta.data)

# run sctransform
pbmc <- Seurat::SCTransform(pbmc, 
                            vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), 
                            verbose = TRUE)

pbmc

#######################

# https://satijalab.org/seurat/articles/pbmc3k_tutorial 

pbmc <- Seurat::RunPCA(pbmc, 
                       features = Seurat::VariableFeatures(pbmc),
                       verbose = FALSE)
pbmc

pbmc[["pca"]]

cell_embedding <- pbmc[["pca"]]@cell.embeddings
cell_embedding[1:5,1:5]
dim(cell_embedding)
Seurat::ElbowPlot(pbmc, ndims = 50)

########################

pbmc <- Seurat::RunUMAP(pbmc, dims = 1:8)

Seurat::FeaturePlot(pbmc,
                    features = "percent.mt",
                    reduction = "umap")
scCustomize::FeaturePlot_scCustom(pbmc,
                                  features = "percent.mt",
                                  reduction = "umap")

genes <- Seurat::VariableFeatures(pbmc)[1:4]
genes
scCustomize::FeaturePlot_scCustom(pbmc,
                                  features = genes,
                                  reduction = "umap")

