library(MEGENA)
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)

n.cores <- 14;
doPar <- TRUE;
method = "spearman"
FDR.cutoff = 0.05
module.pval = 0.05
hub.pval = 0.05

library(DAseq)
library(ggpubr)

cluster_da<-readRDS("cluster_da9.rds")
cnt=cluster_da@assays$RNA@counts

min.cells = ncol(cluster_da) * 0.05
nc = rowSums(cnt >= 1E-320,na.rm = TRUE)
ii = which(nc >= min.cells) # index of expressed genes
print(length(ii))

min.genes = 100;
ng = colSums(cnt[ii,] > 1E-320,na.rm = TRUE)
jj = which(ng >= min.genes)
print(length(jj))

datExpr = cluster_da@assays$RNA@data[ii,jj]
ijw <- calculate.correlation(datExpr = as.matrix(datExpr), doPerm = 10, method = method, FDR.cutoff = FDR.cutoff)

if (doPar)
{
    cl = parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    # check how many workers are there
    cat(paste("number of cores to use:", foreach::getDoParWorkers(), "\n"))
}

el <- calculate.PFN(ijw[, 1:3], doPar = doPar, num.cores = n.cores)
write.table(el, file = "MEGENA_Network.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

g <- graph.data.frame(el,directed = F)

MEGENA.output <- do.MEGENA(g,
                             mod.pval = module.pval, hub.pval = hub.pval, remove.unsig = TRUE,
                             min.size = 20,#max.size = vcount(g)/2,
                             doPar = TRUE, num.cores = n.cores, n.perm = 100,
                             save.output = TRUE)

save(MEGENA.output, file = "MEGENA_output.RData")
q()
