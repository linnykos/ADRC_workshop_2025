# https://satijalab.org/seurat/articles/sctransform_vignette.html
rm(list=ls())
library(Seurat)
library(scCustomize)

pbmc.data <- Seurat::Read10X(data.dir = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/data/seurat-pbmc/filtered_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data, 
                                   project = "pbmc3k")
pbmc

# Run the pipeline to get the UMAP
pbmc <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc <- Seurat::NormalizeData(pbmc) #Seurat::CellCycleScoring requires the "data" slot filled
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pbmc <- Seurat::CellCycleScoring(pbmc, 
                                 s.features = s.genes, 
                                 g2m.features = g2m.genes)

# run sctransform
pbmc <- Seurat::SCTransform(pbmc, 
                            vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), 
                            verbose = TRUE)
pbmc <- Seurat::RunPCA(pbmc, 
                       features = Seurat::VariableFeatures(pbmc),
                       verbose = FALSE)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:8)

# Now run scDblFinder
Seurat::DefaultAssay(pbmc) <- "RNA"
sce <- Seurat::as.SingleCellExperiment(pbmc)
sce.dbl <- scDblFinder::scDblFinder(sce, clusters=NULL) 

pbmc$scDblFinder.score <- sce.dbl@colData$scDblFinder.score
pbmc$scDblFinder.class <- sce.dbl$scDblFinder.class

scCustomize::FeaturePlot_scCustom(pbmc, 
                                  features = "percent.mt",
                                  reduction = "umap")
scCustomize::DimPlot_scCustom(pbmc, 
                              group.by = "scDblFinder.class",
                              reduction = "umap")
scCustomize::DimPlot_scCustom(pbmc, 
                              split.by = "scDblFinder.class",
                              reduction = "umap",
                              split_seurat = TRUE)

scCustomize::FeaturePlot_scCustom(pbmc,
                                  features = "scDblFinder.score",
                                  reduction = "umap")
Seurat::VlnPlot(pbmc, 
                features = "nCount_RNA",
                split.by = "scDblFinder.class")
