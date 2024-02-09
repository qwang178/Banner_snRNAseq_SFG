library(Seurat)
library(DAseq)
library(ggpubr)
library(SingleCellExperiment)

options(future.globals.maxSize = 4000 * 1024^20)
Sys.setenv('R_MAX_VSIZE'=2500000000000)

Oli<-readRDS("Oli.rds")
Idents(Oli) <- "updatedDX"
Oli_DEG.DX<-FindMarkers(Oli, slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "AD", ident.2 = "Control",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="Oli_DEG.DX.csv", Oli_DEG.DX, quote=F)

Ex<-readRDS("Ex.rds")
Idents(Ex) <- "updatedDX"
Ex_DEG.DX<-FindMarkers(Ex, slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "AD", ident.2 = "Control",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="Ex_DEG.DX.csv", Ex_DEG.DX, quote=F)

In<-readRDS("In.rds")
Idents(In) <- "updatedDX"
In_DEG.DX<-FindMarkers(In, slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "AD", ident.2 = "Control",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="In_DEG.DX.csv", In_DEG.DX, quote=F)

Ast<-readRDS("Ast.rds")
Idents(Ast) <- "updatedDX"
Ast_DEG.DX<-FindMarkers(Ast, slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "AD", ident.2 = "Control",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="Ast_DEG.DX.csv", Ast_DEG.DX, quote=F)

Mic<-readRDS("Mic.rds")
Idents(Mic) <- "updatedDX"
Mic_DEG.DX<-FindMarkers(Mic, slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "AD", ident.2 = "Control",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="Mic_DEG.DX.csv", Mic_DEG.DX, quote=F)

Opc<-readRDS("Opc.rds")
Idents(Opc) <- "updatedDX"
Opc_DEG.DX<-FindMarkers(Opc,slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "AD", ident.2 = "Control",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="Opc_DEG.DX.csv", Opc_DEG.DX, quote=F)

EndPer<-readRDS("EndPer.rds")
Idents(EndPer) <- "updatedDX"
EndPer_DEG.DX<-FindMarkers(EndPer, slot = "data", assay = "RNA",  test.use = "MAST", ident.1 = "AD", ident.2 = "Control",logfc.threshold = 0.01,min.pct = 0.01)
write.csv(file="EndPer_DEG.DX.csv", EndPer_DEG.DX, quote=F)

