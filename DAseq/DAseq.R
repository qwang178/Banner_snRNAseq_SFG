library(Seurat)
library(DAseq)
library(ggpubr)

options(future.globals.maxSize = 4000 * 1024^20)
Sys.setenv('R_MAX_VSIZE'=2500000000000)

cluster<-readRDS("combined.sct.dx.rds")
meta<-cluster@meta.data
meta.dx<-meta[,c(1,36)]
labels_res <- meta.dx[meta.dx$updatedDX == "AD", "orig.ident"]
labels_nores <- meta.dx[meta.dx$updatedDX == "Control", "orig.ident"]
labels_res<-as.character(unique(labels_res))
labels_nores<-as.character(unique(labels_nores))

da_cells <- getDAcells(
  X = cluster@reductions$pca@cell.embeddings[,1:40],
  k.vector = seq(100, 4000, 500),
  cell.labels = as.character(cluster@meta.data$orig.ident),
  labels.1 = labels_res, labels.2 = labels_nores,
  plot.embedding = cluster@reductions$umap@cell.embeddings)

plot1<-da_cells$pred.plot

da_cells <- updateDAcells(
  X = da_cells, pred.thres = c(-0.8,0.8),
  plot.embedding = cluster@reductions$umap@cell.embeddings)

plot2<-da_cells$da.cells.plot

da_regions <- getDAregion(
  X = cluster@reductions$pca@cell.embeddings[,1:40],
  da.cells = da_cells,
  cell.labels =as.character(cluster@meta.data$orig.ident),
  labels.1 = labels_res,
  labels.2 = labels_nores,
  resolution = 0.01,
  plot.embedding = cluster@reductions$umap@cell.embeddings)

plot3<-da_regions$da.region.plot

cluster$DA=da_regions$da.region.label

python2use = "/packages/7x/anaconda3/4.4.0/bin/python"
GPU = 4

STG_markers <- STGmarkerFinder(
  X = as.matrix(cluster@assays$RNA@data[,da_cells$cell.idx]),
  da.regions = da_regions,
  lambda = 1.5, n.runs = 5, return.model = T,
  python.use = python2use, GPU = GPU
)

i = 1
while (i <=  NROW(da_regions$DA.stat)) {
    markers<-get(as.character(i),STG_markers$da.markers)
    write.csv(file=paste("integrated.100.DA",i,".markers.csv",sep=""),markers,quote=F)
    i = i+1
}

plot <- ggarrange(plotlist =list(plot1, plot2, plot3),ncol = 3,nrow = 1, 
          labels = c("A", "B","C"))

saveRDS(file="cluster.DA.rds",cluster)
save(list=c("STG_markers", "da_cells", "da_regions"), file="integrated.100.DAseq.RData")
ggsave(file="integrated.100.pdf",plot, h=10, w=40)
