
options(stringsAsFactors = FALSE)

library(reshape2)
library(plyr)
library(data.table)
library(openxlsx)
library(doParallel)
library(limma)
library(parallel)

all_SEG <- read.xlsx(xlsxFile = "Mic.SEG.xlsx", sheet = "all.SEG")
DA9_SEG <- read.xlsx(xlsxFile = "Mic.SEG.xlsx", sheet = "da.SEG")

MEGENA_df <- read.xlsx(xlsxFile = "MEGENA.hub.connectivity.xlsx", sheet = "Mic", colNames = FALSE)

querySets <- list(
                  DA9_MEGENA_hubs = unique(MEGENA_df$X1),
                  DA9_SEG = unique(DA9_SEG$X1)
                   )

source("gene_set_enrichment_fxn.R")

# Enrichment analysis -----------------------------------------------------

classesToGenerate <- grep("Pathways|Ontologies|Transcription|Cell_Types|Diseases_Drugs", dir("enrichr/data/RData/Jan2023", full.names = TRUE), value = TRUE) #gene set libraries downloaded from EnrichR Jan 2023

nCores <- 10

gene_bg <- union(MEGENA_df$X1, union(all_SEG$X1, DA9_SEG$X1)) #run enrichments against union of genes detected in All microglia + DA9 microglia, so that enrichments denote a comparative difference with microglia

enrDir <- "results/enrichments/DA9_microglia/"
dir.create(enrDir, recursive = TRUE)
minSize <- 3
minFDR <- 0.05
convertID <- FALSE

for(category_i in 1:length(classesToGenerate)) {
  
  file_i <- paste0(enrDir, basename(classesToGenerate[category_i]))
  
  if(exists(file_i) == FALSE){
    
    load(classesToGenerate[category_i])
    enrList <- mclapply(X = gmtList, mc.silent = TRUE, mc.cores = nCores, mc.preschedule = TRUE, function(library_i) {
      setpairs2enrichments(set2 = querySets, set1 = library_i, backgroundGenes = gene_bg, minSize = minSize, minFDR = minFDR, convertID = convertID, mapSource = NULL, mapTarget = NULL, printOverlap = TRUE)
    })
    names(enrList) <- names(gmtList)
    save(file = file_i, x = enrList)
    
    #})
  }
}


# Output results ----------------------------------------------------------

resFiles <- dir(enrDir, full.names = TRUE)

resList <-
  lapply(resFiles, function(res_i) {
    load(res_i)
    enr_df <- ldply(enrList, rbind, .id = "Library")
    enr_df
  })
names(resList) <- gsub(basename(resFiles), pattern = ".RData", replacement = "")

openxlsx::write.xlsx(file = "results/Banner_SFG_DA9_microglia_vs_microglia_enrichments.xlsx",x = resList, colNames = TRUE, rowNames = FALSE)
