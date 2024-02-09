library(Seurat)
library(ggplot2)
library(ggpubr)

options(future.globals.maxSize = 4000 * 1024^20)
Sys.setenv('R_MAX_VSIZE'=2500000000000)

dir<-read.csv(file="sample.lst")
dataList=dir$orig.ident

dataList <- lapply(X = dataList, FUN = function(x) {
    print (x)
    x <- readRDS(file=paste0( x,".rds"))
})

#print ("data")

obj.list = dataList
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
prep.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)

#print ("integrate")

anchors <- FindIntegrationAnchors(object.list = prep.list, 
    reference = c(3, 12),
#   20), 
    reduction = "rpca",
    dims = 1:30, normalization.method = "SCT",
    anchor.features = features)

#print ("anchor")

rm(prep.list, obj.list, dataList)

combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30)

#print ("reduced")

combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30) 
combined.sct <- FindClusters(combined.sct, resolution = 0.5)

#print ("clustered")

reference_seurat <-readRDS("/data/breadhea/qwang178/NOMIS/scRNA/Drop8/Seurat/norm/reference.sct.rds")
sim.anchors <- FindTransferAnchors(reference = reference_seurat, query = combined.sct, normalization.method = "SCT",
                                   dims = 1:30)
##replace Group with the actual column name from meta
predictions <- TransferData(anchorset = sim.anchors, refdata = reference_seurat$broad.cell.type,
                            dims = 1:30)
combined.sct <- AddMetaData(object = combined.sct, metadata = predictions)

write.csv(file="combined.csv",combined.sct@meta.data)                                  

saveRDS(file="combined.sct.rds", combined.sct)

