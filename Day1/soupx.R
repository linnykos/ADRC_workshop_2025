# https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html
rm(list=ls())

# there's a few things you need to do:
# 1) make sure there's a folder called "filtered_feature_bc_matrix" and "raw_feature_bc_matrix"
# 2) both folders need the files "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtx.gz"

library(SoupX)
# Load data and estimate soup profile
sc = SoupX::load10X("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/data/seurat-pbmc/")

sc$tod <- sc$toc
sc <- SoupX::estimateSoup(sc, soupRange=c(0,1000))
sc = SoupX::setContaminationFraction(sc, 0.1)

# Clean the data
# takes about 2 minutes
out = SoupX::adjustCounts(sc, clusters = FALSE)

# investigating the change
cntSoggy = rowSums(sc$toc > 0)+1
cntStrained = rowSums(out > 0)+1
mostZeroed = sort((cntSoggy - cntStrained)/cntSoggy)
tail(mostZeroed)
plot(mostZeroed)

out[1:5,1:5]
