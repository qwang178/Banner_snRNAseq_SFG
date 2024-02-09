library(Seurat)
library(DAseq)
library(ggpubr)
library(SingleCellExperiment)

options(future.globals.maxSize = 4000 * 1024^20)
Sys.setenv('R_MAX_VSIZE'=2500000000000)

Oli<-readRDS("Oli.rds")
Idents(Oli) <- "DA"
Oli_DEG.DA<-FindMarkers(Oli, slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "DA", ident.2 = "nonDA",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="Oli_DEG.DA.csv", Oli_DEG.DA, quote=F)

Ex<-readRDS("Ex.rds")
Idents(Ex) <- "DA"
Ex_DEG.DA<-FindMarkers(Ex, slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "DA", ident.2 = "nonDA",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="Ex_DEG.DA.csv", Ex_DEG.DA, quote=F)

In<-readRDS("In.rds")
Idents(In) <- "DA"
In_DEG.DA<-FindMarkers(In, slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "DA", ident.2 = "nonDA",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="In_DEG.DA.csv", In_DEG.DA, quote=F)

Ast<-readRDS("Ast.rds")
Idents(Ast) <- "DA"
Ast_DEG.DA<-FindMarkers(Ast, slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "DA", ident.2 = "nonDA",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="Ast_DEG.DA.csv", Ast_DEG.DA, quote=F)

Mic<-readRDS("Mic.rds")
Idents(Mic) <- "DA"
Mic_DEG.DA<-FindMarkers(Mic, slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "DA", ident.2 = "nonDA",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="Mic_DEG.DA.csv", Mic_DEG.DA, quote=F)

Opc<-readRDS("Opc.rds")
Idents(Opc) <- "DA"
Opc_DEG.DA<-FindMarkers(Opc,slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "DA", ident.2 = "nonDA",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="Opc_DEG.DA.csv", Opc_DEG.DA, quote=F)

EndPer<-readRDS("EndPer.rds")
Idents(EndPer) <- "DA"
EndPer_DEG.DA<-FindMarkers(EndPer, slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "DA", ident.2 = "nonDA",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="EndPer_DEG.DA.csv", EndPer_DEG.DA, quote=F)

