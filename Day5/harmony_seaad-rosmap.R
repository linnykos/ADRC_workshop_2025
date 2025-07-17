rm(list=ls())

library(Seurat)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)

Sys.setenv(R_MAX_VSIZE = 16e9)

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/seaad_microglia.RData")
seaad <- seurat_obj
rm(list = "seurat_obj"); gc(TRUE)
seaad <- Seurat::FindVariableFeatures(seaad, nfeatures = 2000)
seaad <- Seurat::ScaleData(seaad)
seaad <- Seurat::RunPCA(seaad, 
                        features = Seurat::VariableFeatures(seaad),
                        verbose = FALSE)
seaad <- Seurat::RunUMAP(seaad, 
                         dims = 1:30, 
                         verbose = FALSE)
scCustomize::DimPlot_scCustom(seaad,
                              reduction = "umap",
                              group.by = "Supertype")

# from https://compbio.mit.edu/microglia_states/
# Files:
# 1) https://personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.counts.rds?dl=0
# 2) https://personal.broadinstitute.org/cboix/sun_victor_et_al_data/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds?dl=0
# 3) https://cells.ucsc.edu/rosmap-ad-aging-brain/microglia-states/meta.tsv
# construct the ROSMAP object
rosmap_count <- readRDS("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/ROSMAP.ImmuneCells.6regions.snRNAseq.counts.rds")
rosmap_meta <- readRDS("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds")
rosmap_meta2 <- read.csv("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/ROSMAP.meta.tsv",
                         sep = "\t")
rownames(rosmap_meta2) <- rosmap_meta2[,"cellName"]

rosmap <- Seurat::CreateSeuratObject(counts = rosmap_count, 
                                     meta.data = rosmap_meta,
                                     min.cells = 3, 
                                     min.features = 200)
microglia_state <- rosmap_meta2[Seurat::Cells(rosmap), "State"]
rosmap$State <- microglia_state
rosmap <- subset(rosmap, !is.na(rosmap$State))
rm(list = "rosmap_count"); gc(TRUE)

# Do a simple processing of rosmap
rosmap <- Seurat::NormalizeData(rosmap)
rosmap <- Seurat::FindVariableFeatures(rosmap, 
                                       selection.method = "vst", 
                                       nfeatures = 2000)
rosmap <- Seurat::ScaleData(rosmap)
rosmap <- Seurat::RunPCA(rosmap, 
                         features = Seurat::VariableFeatures(rosmap),
                         verbose = FALSE)
rosmap <- Seurat::RunUMAP(rosmap, 
                          dims = 1:30, 
                          verbose = FALSE)
scCustomize::DimPlot_scCustom(rosmap,
                              reduction = "umap",
                              group.by = "State")

########################################

# Let's now try to integrate
# Using a simpler workflow from https://satijalab.org/seurat/archive/v4.3/integration_introduction

# Add dataset-specific prefixes so barcodes don't collide
seaad  <- Seurat::RenameCells(seaad,  add.cell.id = "seaad")
rosmap <- Seurat::RenameCells(rosmap, add.cell.id = "rosmap")

common_genes <- intersect(rownames(seaad), rownames(rosmap))
seaad  <- subset(seaad,  features = common_genes)[common_genes, ]
rosmap <- subset(rosmap,  features = common_genes)[common_genes, ]

anchors <- Seurat::FindTransferAnchors(
  reference           = seaad,
  query               = rosmap,
  normalization.method = "LogNormalize",
  reference.reduction  = "pca",   
  reduction            = "pcaproject",
  dims                 = 1:30      
)

rosmap <- Seurat::IntegrateEmbeddings(
  anchorset           = anchors,
  reference           = seaad,
  query               = rosmap,
  reductions          = "pcaproject",      
  new.reduction.name  = "integrated.pca",
  dims.to.integrate   = 1:30
)

# Now visualize
# copy the PCA reduction so both objects agree on the name
seaad[["integrated.pca"]] <- seaad[["pca"]]
seaad <- Seurat::DietSeurat(
  object      = seaad,
  features    = common_genes,
  assays      = "RNA",
  counts      = TRUE,   data = FALSE,  scale.data = FALSE,
  dimreducs   = "integrated.pca"
)
rosmap <- Seurat::DietSeurat(
  object      = rosmap,
  features    = common_genes,
  assays      = "RNA",
  counts      = TRUE,   data = FALSE,  scale.data = FALSE,
  dimreducs   = "integrated.pca"
)

combo <- merge(seaad, 
               rosmap)
dataset_vec <- rep(NA, length(Seurat::Cells(combo)))
names(dataset_vec) <- Seurat::Cells(combo)
dataset_vec[Seurat::Cells(seaad)] <- "SEAAD"
dataset_vec[Seurat::Cells(rosmap)] <- "ROSMAP"
combo$dataset <- dataset_vec

# Manually construct the integrated.pca assay
dimred_mat <- matrix(NA,
                     nrow = length(Seurat::Cells(combo)),
                     ncol = ncol(rosmap[["integrated.pca"]]@cell.embeddings),
                     dimnames = list(
                       Seurat::Cells(combo),
                       colnames(rosmap[["integrated.pca"]]@cell.embeddings)
                     ))
dimred_mat[Seurat::Cells(seaad),] <- seaad[["integrated.pca"]]@cell.embeddings[,colnames(dimred_mat)]
dimred_mat[Seurat::Cells(rosmap),] <- rosmap[["integrated.pca"]]@cell.embeddings[,colnames(dimred_mat)]
combo[["integrated.pca"]] <- Seurat::CreateDimReducObject(dimred_mat)

combo <- Seurat::RunUMAP(
  object     = combo,
  reduction  = "integrated.pca",
  dims       = 1:30,
  reduction.name = "umap.integrated"
)

Seurat::Idents(combo) <- "dataset"
scCustomize::DimPlot_scCustom(combo,
                              reduction = "umap.integrated",
                              group.by = "dataset")
scCustomize::DimPlot_scCustom(combo,
                              reduction = "umap.integrated",
                              split.by = "dataset",
                              split_seurat = TRUE)

scCustomize::DimPlot_scCustom(combo,
                              reduction = "umap.integrated",
                              split.by = "Supertype",
                              split_seurat = TRUE)

