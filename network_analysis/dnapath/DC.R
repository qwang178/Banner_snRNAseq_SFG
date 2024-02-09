library(MEGENA)
library(Matrix)
library(rstatix)
library(dnapath)
library(Seurat)
source("init_required_packages_and_files.R")                 

method = "spearman"
FDR.cutoff = 0.05
nodes_cor<-read.csv(file="all.nodes.cor.txt",sep="",header=F)
colnames(nodes_cor)=c("var1","var2","rho")                                  

network_inference = function(x) {
    ijw <- calculate.correlation(datExpr = as.matrix(x), doPerm = 10, method = method, FDR.cutoff = FDR.cutoff)
    colnames(ijw)=c("var1","var2","rho")
    ijw_cor=rbind(nodes_cor,ijw)
    ijw_cor_mtx<-as.matrix(cor_spread(ijw_cor, value = "rho"))
    ijw_cor_mtx<-trimws(ijw_cor_mtx)
    rownames(ijw_cor_mtx)=ijw_cor_mtx[,1]
    ijw_cor_mtx=ijw_cor_mtx[,-1]
    class(ijw_cor_mtx) <- "numeric"
    ijw_cor_mtx[is.na(ijw_cor_mtx)] <- 0
    scores <- standardize_scores(ijw_cor_mtx, ignore_zeroes = F, robust = F)    
    rownames(scores)<-rownames(ijw_cor_mtx)
    colnames(scores)<-colnames(ijw_cor_mtx)
    return(scores)
}
sub_nodes<-read.csv(file="../../MEGENA/Mic/sub/c1_2.nodes.txt", header=T)
da_nodes<-read.csv(file="../../MEGENA/Mic/da/c1_5.nodes.txt", header=T)
colnames(da_nodes)=c("node")
colnames(sub_nodes)=c("node")
nodes<-unique(rbind(da_nodes,sub_nodes))
nodes_cor<-subset(nodes_cor, var1 %in% nodes$node)
nodes_cor<-nodes_cor[order(nodes_cor$var1),]

sub<-readRDS("../../comp_objects/Mic.sub.rds")
sub<-subset(sub,comp_sub=="sub")
subE<-sub@assays$RNA@data
sube<-subE[rownames(subE) %in% nodes$node,]
da<-readRDS("../../da_objects/cluster_da9.rds")
daE = da@assays$RNA@data
dae<-daE[rownames(daE) %in% nodes$node,]
score_da<-network_inference(dae)
score_sub<-network_inference(sube)

counts_joined<-cbind(dae,sube)                                                  
d_gene <- d_genesC(score_sub, score_da, lp = 2)
d_edge <- d_edgesC(score_sub, score_da, lp = 2)
d_gene_ranks <- order(d_gene, decreasing = TRUE)
d_edge_ranks <- order(d_edge, decreasing = TRUE)
p <- nrow(counts_joined)
N <- ncol(counts_joined)
pval_gene <- rep(1, p)
pval_gene_mono <- rep(1, p)
pval_edge_mono <- rep(1, p * (p - 1) / 2)
pval_edge <- rep(1, p * (p - 1) / 2)
n <- c(ncol(daE) - 1, ncol(subE) - 1)

n_perm=1000
if(n[1] == n[2] & choose(N, n[1]) / 2 <= n_perm) {
    permutations <- combn(1:N, n[1])
    permutations <- permutations[, 1:(ncol(permutations) / 2)] 
    permutations <- permutations[, -1]
    # cat("Iterating over all", ncol(permutations), "permutations.\n")
  } else if(n[1] != n[2] & choose(N, n[1]) <= n_perm) {
    permutations <- combn(1:N, n[1])
    permutations <- permutations[, -1]
    # cat("Iterating over all", ncol(permutations), "permutations.\n")
  } else {
    permutations <- cbind(sapply(1:n_perm, function(x) sample(1:N, n[1])))
    # cat("Iterating over", ncol(permutations), "random permutations.\n")
  }

for(i in 1:n_perm) {
      network_1 <- permutations[, i]
      scores_1 <- network_inference(counts_joined[,network_1])
      scores_1[which(is.na(scores_1))] <- 0
      scores_2 <- network_inference(counts_joined[,-network_1])
      scores_2[which(is.na(scores_2))] <- 0
pval_gene <- pval_gene + (d_genesC(scores_1, scores_2, lp = 2) >= d_gene) 
pval_edge <- pval_edge + (d_edgesC(scores_1, scores_2, lp = 2) >= d_edge) 
index <- d_gene_ranks
n_index <- length(index)
        d_gene_perm <- d_genesC(scores_1, scores_2, lp = 2)
        d_gene_perm <- cummax(d_gene_perm[rev(index)])[n_index:1][order(index)]
        pval_gene_mono <- pval_gene_mono +    (d_gene_perm >= d_gene)
index <- d_edge_ranks
n_index <- length(index)
        d_edge_perm <- d_edgesC(scores_1, scores_2, lp = 2)
        d_edge_perm <- cummax(d_edge_perm[rev(index)])[n_index:1][order(index)]
        pval_edge_mono <- pval_edge_mono +    (d_edge_perm >= d_edge)
print(i)
}

pval_gene <- pval_gene/ (n_perm + 1)
pval_edge <- pval_edge/ (n_perm + 1)
rownames(pval_gene)=rownames(score_da)
rownames(d_gene)=rownames(score_da)
DC<-cbind(d_gene,pval_gene)
colnames(DC)=c("DC","pval")
write.csv(file="DC.csv",DC,quote=F)

