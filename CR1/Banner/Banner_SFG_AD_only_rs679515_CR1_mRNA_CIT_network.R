
options(stringsAsFactors = FALSE)

library(limma)
library(edgeR)
library(plyr)
library(data.table)
library(cit)
library(doParallel)
library(pbapply)
library(WGCNA)
library(parallel)

# Read in Oligo snRNAseq expression, aggreagted to sample level -------------------------------------------------
exprCounts <- fread(input = "Oli.counts.bySample.csv", data.table = FALSE)
rownames(exprCounts) <- exprCounts$V1
exprCounts <- exprCounts[,-1]
y <- DGEList(exprCounts)
y <- calcNormFactors(y)
keep <- filterByExpr(y)
v <- voom(y)[keep,]
normedExpr <- v$E

# Read in table of sample ID, diagnosis & rs679515 dosage ----------------------------------------------------------

summary_df <- fread(file = "Banner_rs679515_dosage.txt", data.table = FALSE, showProgress = TRUE)

# Check other genes that correlate with CR1 to scope CIT analysis -------------------------------
sample_scope_i <- intersect(summary_df$Sample_ID[which(summary_df$updatedDX == "AD")], colnames(v$E))


corList <- corAndPvalue(x = v$E["CR1",sample_scope_i], y = t(v$E[,sample_scope_i]), method = "pearson")
cor_df <- data.frame(melt(corList$cor)[,2:3], melt(corList$p)[,3])
colnames(cor_df) <- c("Symbol", "Cor", "Pvalue")
cor_df <- cor_df[order(cor_df$Pvalue, decreasing = FALSE),]
cor_df$FDR_BH <- p.adjust(cor_df$Pvalue, method = "BH")

candidate_T <- as.character(cor_df$Symbol[which(cor_df$FDR_BH <= 0.05)])[-1]


nCores <- 20
nPerm <- 1000

mediator2trait_CIT_list <- mclapply(X = candidate_T, mc.silent = TRUE, mc.cores = nCores, mc.preschedule = TRUE, mc.set.seed = 12345, function(T_i) {
                                     l_i <- summary_df$dosage[match(sample_scope_i, summary_df$Sample_ID)]
                                     g_i <- v$E["CR1",sample_scope_i]
                                     t_i <- v$E[T_i,sample_scope_i]

                                     fdr.cit(list(res1 = cit.cp(L = l_i, G = g_i, T = t_i, n.perm = nPerm)))
})
names(mediator2trait_CIT_list) <- candidate_T

mediator2trait_CIT_df <- ldply(mediator2trait_CIT_list, rbind, .id = "Transgene")
mediator2trait_CIT_df <- mediator2trait_CIT_df[order(mediator2trait_CIT_df$q.cit, decreasing = FALSE),]
mediator2trait_CIT_df$gene_name <- sym2name2entrez$gene_name[match(mediator2trait_CIT_df$Transgene, sym2name2entrez$symbol)]

mediator2trait_CIT_df <- data.frame(mediator2trait_CIT_df, cor_df[match(mediator2trait_CIT_df$Transgene, cor_df$Symbol),c("Cor", "Pvalue", "FDR_BH")])
mediator2trait_CIT_df <- data.frame(Variant = "1:207577223", Cisgene = "CR1", mediator2trait_CIT_df)

colnames(mediator2trait_CIT_df) <- c("Variant", "Cisgene", "Transgene","p.cit","q.cit","q.cit.ll","q.cit.ul","q.TaL","q.ll.TaL","q.ul.TaL","q.TaGgvL","q.ll.TaGgvL","q.ul.TaGgvL","q.GaLgvT","q.ll.GaLgvT","q.ul.GaLgvT","q.LiTgvG","q.ll.LiTgvG","q.ul.LiTgvG","p_TassocL","p_TassocGgvnL","p_GassocLgvnT","p_LindTgvnG","gene_name","Cis2transgene_Cor","Cis2transgene_Pvalue","FDR_BH")
mediator2trait_CIT_df <- mediator2trait_CIT_df[,c("Variant", "Cisgene", "Transgene","gene_name", "p.cit","q.cit","q.cit.ll","q.cit.ul","q.TaL","q.ll.TaL","q.ul.TaL","q.TaGgvL","q.ll.TaGgvL","q.ul.TaGgvL","q.GaLgvT","q.ll.GaLgvT","q.ul.GaLgvT","q.LiTgvG","q.ll.LiTgvG","q.ul.LiTgvG","p_TassocL","p_TassocGgvnL","p_GassocLgvnT","p_LindTgvnG","Cis2transgene_Cor","Cis2transgene_Pvalue","FDR_BH")]

o_i <- list(CR1_AD_GWAS_CIT = mediator2trait_CIT_df)

openxlsx::write.xlsx(file = "CR1_AD_GWAS_rs679515_AD_only_CIT_summary_v1.xlsx", x = o_i, colNames = TRUE, rowNames = FALSE)

