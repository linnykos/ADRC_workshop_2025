# https://github.com/constantAmateur/SoupX/blob/master/R/load10X.R
dataDir = "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/sumie-katie/data/seurat-pbmc"
cellIDs=NULL
channelName=NULL
readArgs=list()
includeFeatures=c('Gene Expression')
verbose=TRUE

#Work out which version we're dealing with
isV3 = dir.exists(file.path(dataDir,'raw_feature_bc_matrix'))
isV7 = dir.exists(file.path(dataDir,'analysis','clustering','gene_expression_graphclust'))
isMulti = dir.exists(file.path(dataDir,'analysis', 'clustering', 'gex'))
tgt = file.path(dataDir,
                ifelse(isV3,'raw_feature_bc_matrix','raw_gene_bc_matrices'))
#Add the reference genome for the non-V3 ones
if(!isV3)
  tgt = file.path(tgt,list.files(tgt))
if(verbose)
  message(sprintf("Loading raw count data"))
dat = do.call(Read10X,c(list(data.dir=tgt),readArgs))
if(verbose)
  message(sprintf("Loading cell-only count data"))

if(!is.null(cellIDs)){
  #This is now redundant as we require the version of Seurat that does not strip the suffix
  ####Do the same sripping that Seurat does on IDs
  ###if(all(grepl('\\-1$',cellIDs)))
  ###  cellIDs = gsub('\\-1$','',cellIDs)
  #Check we have the IDs
  if(!all(cellIDs %in% colnames(dat)))
    stop("Not all supplied cellIDs found in raw data.")
  datCells = dat[,match(cellIDs,colnames(dat))]
}else{
  #Work out which ones contain cells
  tgt = file.path(dataDir,
                  ifelse(isV3,'filtered_feature_bc_matrix','filtered_gene_bc_matrices'))
  if(!isV3) tgt = file.path(tgt,list.files(tgt))
  datCells = do.call(Read10X,c(list(data.dir=tgt),readArgs))
  #If it's a list of multiple types, have to decide what to include and collapse to one matrix.
  if(is.list(dat)){
    dat = do.call(rbind,dat[includeFeatures])
    datCells = do.call(rbind,datCells[includeFeatures])
  }
}
