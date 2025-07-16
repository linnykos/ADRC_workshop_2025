rm(list=ls())

library(Seurat)
library(DESeq2)

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/out/ADRC_workshop_2025/seaad_microglia.RData")

seurat_obj <- subset(seurat_obj, ADNC %in% c("High", "Intermediate", "Low", "Not AD"))
seurat_obj <- subset(seurat_obj, assay == "10x 3' v3")
seurat_obj <- subset(seurat_obj, Supertype == "Micro-PVM_2")

table(seurat_obj$donor_id)

# Find variable genes for simplicity
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, 
                                           selection.method = "vst", 
                                           nfeatures = 2000)

pseudo_seurat <- Seurat::AggregateExpression(seurat_obj, 
                                             assays = "RNA", 
                                             return.seurat = TRUE, 
                                             group.by = c("donor_id", "APOE4.status", "ADNC", "sex"))
Seurat::VariableFeatures(pseudo_seurat) <- Seurat::VariableFeatures(seurat_obj)
pseudo_seurat
head(pseudo_seurat@meta.data)

###########

mat_pseudobulk <- SeuratObject::LayerData(pseudo_seurat,
                                          layer = "counts",
                                          assay = "RNA",
                                          features = Seurat::VariableFeatures(pseudo_seurat))
metadata_pseudobulk <- pseudo_seurat@meta.data

# do some simple processing of the metadata
metadata_pseudobulk[,"donor_id"] <- factor(metadata_pseudobulk[,"donor_id"])
metadata_pseudobulk[,"APOE4.status"] <- factor(metadata_pseudobulk[,"APOE4.status"])
metadata_pseudobulk[,"sex"] <- factor(metadata_pseudobulk[,"sex"])
adnc_vec <- rep("NonAD", nrow(metadata_pseudobulk))
adnc_vec[metadata_pseudobulk[,"ADNC"] %in% c("High", "Intermediate")] <- "AD"
adnc_vec <- factor(adnc_vec)
metadata_pseudobulk[,"ADNC"] <- relevel(adnc_vec, ref = "NonAD")

mat_pseudobulk[1:5,1:5]
dim(mat_pseudobulk)
summary(metadata_pseudobulk)
dim(metadata_pseudobulk)

# do DESeq2
dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_pseudobulk,
                                      colData = metadata_pseudobulk,
                                      design = ~ ADNC + sex)
dds <- DESeq2::DESeq(dds)
deseq2_res <- DESeq2::results(dds, name="ADNC_AD_vs_NonAD")

# pick a multiple-testing adjusted threshold
tmp <- subset(deseq2_res, padj >= 0.05)
y_threshold <- min(tmp$pvalue)

EnhancedVolcano::EnhancedVolcano(deseq2_res,
                                 lab = rownames(deseq2_res),
                                 x = "log2FoldChange",
                                 y = "pvalue",
                                 pCutoff = y_threshold)





