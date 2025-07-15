# https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html
rm(list=ls())
library(SoupX)

# there's a few things you need to do:
# 1) make sure there's a folder called "filtered_feature_bc_matrix" and "raw_feature_bc_matrix"
# 2) both folders need the files "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtx.gz"

# Load data and estimate soup profile
sc <- SoupX::load10X("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/data/seurat-pbmc/")

# the following lines is NOT how I could advise to run SoupX.
## This is purely pedagogical to run SoupX without a clustering
sc$tod <- sc$toc
sc <- SoupX::estimateSoup(sc, soupRange=c(0,1000))
sc <- SoupX::setContaminationFraction(sc, 0.1)

# Clean the data: takes about 2 minutes
out <- SoupX::adjustCounts(sc, clusters = FALSE)

# investigating the change
cntSoggy <- rowSums(sc$toc > 0)+1
cntStrained <- rowSums(out > 0)+1
mostZeroed <- sort((cntSoggy - cntStrained)/cntSoggy)
tail(mostZeroed)
plot(mostZeroed)

out[1:5,1:5]

############################

# if you had a Seurat object, then this is how you can integrate the 
## SoupX output into the Seurat object

pbmc.data <- Seurat::Read10X(data.dir = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/data/seurat-pbmc/filtered_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k")
pbmc

out <- out[,Seurat::Cells(pbmc)]
pbmc[["RNA_cleaned"]] <- Seurat::CreateAssayObject(counts = out)
Seurat::DefaultAssay(pbmc) <- "RNA_cleaned"
pbmc

head(pbmc@meta.data)

plot(x = pbmc$nCount_RNA, 
     y = pbmc$nCount_RNA_cleaned, 
     asp = TRUE,
     pch = 16,
     col = rgb(0.5, 0.5, 0.5, 0.5))
lines(c(0, 1e6), c(0, 1e6), col = "red", lwd = 2, lty = 2)

count_raw <- SeuratObject::LayerData(pbmc,
                                     layer = "counts",
                                     assay = "RNA")
gene_counts_raw <- Matrix::rowSums(count_raw)
count_soupx <- SeuratObject::LayerData(pbmc,
                                       layer = "counts",
                                       assay = "RNA_cleaned")
gene_counts_soupx <- Matrix::rowSums(count_soupx)

plot(x = gene_counts_raw, 
     y = gene_counts_soupx, 
     asp = TRUE,
     pch = 16,
     col = rgb(0.5, 0.5, 0.5, 0.5))
lines(c(0, 1e6), c(0, 1e6), col = "red", lwd = 2, lty = 2)



