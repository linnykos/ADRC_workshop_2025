rm(list=ls())
library(Seurat)
library(scDblFinder)
pbmc.data <- Seurat::Read10X(data.dir = "/Users/kevinlin/Downloads/filtered_gene_bc_matrices/hg19")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k")
pbmc

# https://satijalab.org/seurat/reference/as.singlecellexperiment
Seurat::DefaultAssay(pbmc) <- "RNA"
sce <- Seurat::as.SingleCellExperiment(pbmc)

# This procedure is a bit "better" if you have cell cluster information
## This code below though is simply for pedagogical demonstration
# https://bioconductor.org/books/3.15/OSCA.advanced/doublet-detection.html
set.seed(10)
sce.dbl <- scDblFinder::scDblFinder(sce, clusters=NULL) # takes about 1 minute

# Take a look at the doublet scores
head(sce.dbl@colData$scDblFinder.score)
quantile(sce.dbl@colData$scDblFinder.score)
hist(sce.dbl@colData$scDblFinder.score)

# scDblFinder also enumerates the cells with "too high" of a doublet score
table(sce.dbl$scDblFinder.class)

pbmc$scDblFinder.score <- sce.dbl@colData$scDblFinder.score
pbmc$scDblFinder.class <- sce.dbl$scDblFinder.class

pbmc <- subset(pbmc, scDblFinder.class == "singlet")
pbmc
