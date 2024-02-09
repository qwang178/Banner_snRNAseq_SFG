library(rhdf5)
library(qvalue)
library(dplyr)

observedFeatures <- 0
results <- NULL
snpAnnotation <- NULL
featureAnnotation <- NULL

tmp=h5dump(file = "perm_results_22.h5")
if(length(tmp)>0){
    for (j in names(tmp)) tmp[[j]][["feature"]] <- j
    observedFeatures = observedFeatures+length(tmp)
    df <- bind_rows(tmp)
    if(nrow(df)>0){
      results = rbind(results,df)
    }
}
results["QTL"] <- paste(results$snp_id, results$feature,sep="-")
if(length(which(duplicated(results$QTL)))>0){
  results <- results[-which(duplicated(results$QTL)),]
}
write.table("permutationInformation.txt",x = results,sep="\t",row.names=F,quote=F)

library(matrixStats)
qtl<-read.csv(file="qtl_results_22.txt",sep="")
perm<-left_join(qtl,results, by=c("feature_id"="feature","snp_id"))
write.table("permutationInformation.sorted.txt",x = perm,sep="\t",row.names=F,quote=F)

p<-perm[,22:121]
p_min<-rowMins(as.matrix(p))
write.table(file="Permutation.pValues.top.txt",p_min, row.names=F,quote=F,col.names=F)

