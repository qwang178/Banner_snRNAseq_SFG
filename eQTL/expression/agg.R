library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(Seurat)
options(future.globals.maxSize = 4000 * 1024^20)
Sys.setenv('R_MAX_VSIZE'=2500000000000)
library(scran)

celltype<-readRDS("Oli.rds")
DefaultAssay(celltype)<-"RNA"
celltype.sce <- as.SingleCellExperiment(celltype)

clusters <- quickCluster(celltype.sce)
celltype.sce <- computeSumFactors(celltype.sce, clusters=clusters)
summary(sizeFactors(celltype.sce))
celltype.sce <- logNormCounts(celltype.sce)
agg.sce <- aggregateAcrossCells(celltype.sce, ids= celltype.sce$orig.ident, statistics="mean", use.assay.type = "logcounts")
agg<-logcounts(agg.sce)

NOMIS.annotation<-read.csv(file="../NOMIS.annotation.tsv",sep="")
agg.complete<-agg[NOMIS.annotation$feature_id,]
PC.complete<-prcomp(agg.complete)
write.csv(file="common.PC.complete.csv",PC.complete$rotation)
write.table(agg.complete, file = "e.tsv", sep = "\t", row.names = TRUE, col.names = TRUE,quote=F)

PC<-PC.complete$rotation[,1:20]
cov<-read.csv(file="../covariates.tsv",sep="")
cov0<-merge(PC, cov, by=0, all=TRUE)
write.table(cov0, file = "covariates0.tsv", sep = "\t", row.names = F, col.names = TRUE,quote=F)
