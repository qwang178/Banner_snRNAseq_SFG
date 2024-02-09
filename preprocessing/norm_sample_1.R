library(Seurat)
library(ggplot2)
library(sctransform)
library(ggpubr)
library(DoubletFinder)
library(loomR)
library(patchwork)
library(SeuratData)
library(SeuratDisk)

data <- Read10X(data.dir = "/data/breadhea/qwang178/NOMIS/scRNA/Drop8/cellranger_outputs/sample_1/outs/filtered_feature_bc_matrix/")
data <- CreateSeuratObject(counts = data, project = "sample_1", min.cells = 3, min.features = 200)
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")

data <- PercentageFeatureSet(data, "^RP[SL]", col.name = "percent.ribo")
data <- PercentageFeatureSet(data, "^HB[^(P)]", col.name = "percent.hb")
data <- PercentageFeatureSet(data, "PECAM1|PF4", col.name = "percent.plat")

x <- log10(sort(data@meta.data$nCount_RNA))
y <- seq_along(x)
fit <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = data.frame(x, y))

b<-coef(fit)["xmid"]                                                            
c<-coef(fit)["scal"]                                                            
                                                                                
Nmin <- b - c * log(4)                                                          
Nmax <- b + c * log(4)                                                          
                                                                                
Nmin <- 10^Nmin                                                                 
Nmax <- 10^Nmax                                                                 

data <- subset(x = data, 
       subset = (nCount_RNA >= Nmin) & 
                (nCount_RNA <= Nmax) & 
                (log10GenesPerUMI > 0.80))

selected_mito <- WhichCells(data, expression = percent.mt < 5)
selected_ribo <- WhichCells(data, expression = percent.ribo < 5)
selected_hb <- WhichCells(data, expression = percent.hb < 0.1) 
data <- subset(data, cells = selected_mito)
data <- subset(data, cells = selected_ribo)
data <- subset(data, cells = selected_hb)  

data.filt <- data

# Filter MALAT1
data.filt <- data.filt[!grepl("MALAT1", rownames(data.filt)), ]

# Filter Mitocondrial
data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]

# Filter Ribossomal gene (optional if that is a problem on your data) data.filt
data.filt <- data.filt[ ! grepl('^RP[SL]', rownames(data.filt)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
data.filt <- data.filt[!grepl("^HB[^(P)]", rownames(data.filt)), ]

sample_1 <- data.filt

# cell cycle identification and regression
cell_cycle_markers = read.csv(file="/data/breadhea/qwang178/NOMIS/scRNA/Drop8/Seurat/norm/cell_cycle_markers.csv")
s_genes <- cell_cycle_markers %>%
        dplyr::filter(phase == "S") %>%
        dplyr::pull("gene_name")
      
g2m_genes <- cell_cycle_markers %>%
        dplyr::filter(phase == "G2/M") %>%
        dplyr::pull("gene_name")

sample_1 <- NormalizeData(sample_1)
sample_1 <- CellCycleScoring(sample_1, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)

# run sctransform
sample_1 <- SCTransform(sample_1,  vst.flavor = "v2",
                        vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA", "S.Score", "G2M.Score"),
                        verbose = FALSE)

#Perform dimensionality reduction by PCA and UMAP embedding                     
# These are now standard steps in the Seurat workflow for visualization and clustering
sample_1 <- RunPCA(sample_1, verbose = FALSE)                                   
sample_1 <- RunUMAP(sample_1, dims = 1:30, verbose = FALSE)                     

# run parameter optimization with paramSweep

sweep.res <- paramSweep_v3(sample_1, PCs = 1:10, sct = T) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats) 

BCmetric=bcmvn$BCmetric
pK=as.numeric(as.character(bcmvn$pK))
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

# define the expected number of doublet cells.
nExp <- round(ncol(sample_1) * 0.075/10000 * ncol(sample_1))  # expect doublets
sample_1 <- doubletFinder_v3(sample_1, pN = 0.25, pK = pK_choose, nExp = nExp, PCs = 1:10, sct = T)

# remove doublets
DF.name = colnames(sample_1@meta.data)[grepl("DF.classification", colnames(sample_1@meta.data))]
sample_1 = sample_1[, sample_1@meta.data[, DF.name] == "Singlet"]

# clustering
sample_1 <- FindNeighbors(sample_1, dims = 1:30, verbose = FALSE)
sample_1 <- FindClusters(sample_1, verbose = FALSE)

# annotation by Mathys
reference_seurat <-readRDS("/data/breadhea/qwang178/NOMIS/scRNA/Drop8/Seurat/norm/reference.sct.rds")
sim.anchors <- FindTransferAnchors(reference = reference_seurat, query = sample_1, normalization.method = "SCT",
                                   dims = 1:30)

##replace Group with the actual column name from meta
predictions <- TransferData(anchorset = sim.anchors, refdata = reference_seurat$broad.cell.type,
                            dims = 1:30)
sample_1 <- AddMetaData(object = sample_1, metadata = predictions)

# save files
write.csv(file="sample_1.csv",sample_1@meta.data)                                  
saveRDS(file="sample_1.rds",sample_1)
SaveH5Seurat(sample_1, filename = "sample_1.h5Seurat")
Convert("sample_1.h5Seurat", dest = "h5ad")

