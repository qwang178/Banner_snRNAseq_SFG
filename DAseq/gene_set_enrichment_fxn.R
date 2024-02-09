
setpairs2enrichments <- function(set1, set2, backgroundGenes, minSize, printOverlap = FALSE, alternative = "over", minFDR = 1, convertID = FALSE, mapSource = NULL, mapTarget = NULL){
  
  #require("HTSanalyzeR")
  
  set1 <- lapply(set1, function(x) intersect(x, backgroundGenes))
  set2 <- lapply(set2, function(x) intersect(x, backgroundGenes))
  
  set1 <- set1[which(sapply(set1,length) >= minSize)]
  set2 <- set2[which(sapply(set2,length) >= minSize)]
  
  if(printOverlap == FALSE){
    
    modEnrPerGeneSetClass <- vector("list", length(set1))
    names(modEnrPerGeneSetClass) <- names(set1)
    
    ###############
    for(setType_i in 1:length(set1)){
      
      datAll <- data.frame(matrix(data = NA, nrow = 0,ncol = 9))
      dat <- multiHyperGeoTest(set2, backgroundGenes, set1[[setType_i]], pAdjustMethod = "BH", verbose = TRUE, minGeneSetSize = minSize)
      dat <- data.frame(feature = rownames(dat), dat,row.names = NULL)
      
      if(nrow(dat) > 0){
        
        dat <- data.frame(module = names(set1)[setType_i], dat,row.names=NULL)
        
        datAll <- rbind(datAll, dat)
        
      }
      
      modEnrPerGeneSetClass[[setType_i]] <- datAll
      
    }
    
    names(modEnrPerGeneSetClass) <- names(set1)
    
  }
  
  if(printOverlap == TRUE){
    
    modEnrPerGeneSetClass <- vector("list", length(set1))
    names(modEnrPerGeneSetClass) <- names(set1)
    
    ###############
    for(setType_i in 1:length(set1)){
      
      datAll <- data.frame(matrix(data = NA, nrow = 0,ncol = 11))
      dat <- multiHyperGeoTest(set2, backgroundGenes, set1[[setType_i]], pAdjustMethod = "BH", verbose = TRUE, minGeneSetSize = minSize)
      dat <- data.frame(feature = rownames(dat), dat,row.names = NULL)
      
      set2_reordered <- set2[match(dat$feature, names(set2))]
      
      if(nrow(dat) > 0){
        
        overlap <- sapply(set2_reordered, function(x) intersect(x, set1[[setType_i]]))
        overlap[sapply(overlap, function(x) length(x) == 0)] <- ""
        
        dat <- data.frame(module = names(set1)[setType_i], dat,enrichment_overlap = unlist(sapply(overlap, function(x) paste(sort(x), collapse = "|"))), row.names=NULL)
        
        if(convertID == TRUE){
          
          dat$enrichment_overlap_convertedID <- sapply(overlap, function(x) paste(sort(mapTarget[match(x, mapSource)]), collapse = "|"))
          
        }
        
        datAll <- rbind(datAll, dat)
        
      }
      
      modEnrPerGeneSetClass[[setType_i]] <- datAll
      
    }
    
    names(modEnrPerGeneSetClass) <- names(set1)
    
  }
  
  
  allEnr <- plyr::ldply(modEnrPerGeneSetClass[as.numeric(which(sapply(modEnrPerGeneSetClass,nrow) > 0))], rbind)
  colnames(allEnr)[1] <- "set class"
  allEnr$foldchange <- allEnr$Observed.Hits/allEnr$Expected.Hits
  allEnr <- allEnr[order(allEnr[,paste("Pvalue", alternative, sep = "_")], decreasing = FALSE),]
  allEnr <- allEnr[which(allEnr[,paste("Adjusted_Pvalue", alternative, sep = "_")] <= minFDR),]
  
  return(allEnr)
  
}

hyperGeoTest <- function (geneSet, universe, hits) {
  N <- length(universe)
  geneSet <- intersect(geneSet[[1]], universe)
  m <- length(geneSet)
  Nm <- N - m
  overlap <- intersect(geneSet, hits)
  k <- length(overlap)
  n <- length(hits)
  
  #see: http://mengnote.blogspot.qa/2012/12/calculate-correct-hypergeometric-p.html
  HGTresults_overrepresented <- phyper(k - 1, m, Nm, n, lower.tail = F) #over
  HGTresults_underrepresented <- phyper(k, m, Nm, n, lower.tail= TRUE) #under
  HGTresults_twotailed <- 2 * min(HGTresults_overrepresented, HGTresults_underrepresented)
  HGTresults_twotailed[HGTresults_twotailed > 1] <- 1
  ex <- (n/N) * m
  if (m == 0) 
    HGTresults <- NA
  hyp.vec <- c(N, m, n, ex, k, HGTresults_overrepresented, HGTresults_underrepresented, HGTresults_twotailed)
  names(hyp.vec) <- c("Universe Size", "Gene Set Size", "Total Hits", 
                      "Expected Hits", "Observed Hits", "Pvalue_over", "Pvalue_under", "Pvalue_twosided")
  return(hyp.vec)
}

multiHyperGeoTest <- function (collectionOfGeneSets, universe, hits, minGeneSetSize = 15, pAdjustMethod = "BH", verbose = TRUE) {
  l.GeneSet <- length(collectionOfGeneSets)
  geneset.size <- unlist(lapply(lapply(collectionOfGeneSets, 
                                       intersect, y = universe), length))
  if (all(geneset.size < minGeneSetSize)) 
    stop(paste("The largest number of overlapped genes of gene ", 
               "sets with universe is: ", max(geneset.size), ", which is < ", 
               minGeneSetSize, "!\n", sep = ""))
  geneset.filtered <- which(geneset.size >= minGeneSetSize)
  if (verbose) 
    pb <- txtProgressBar(style = 3)
  results <- t(sapply(geneset.filtered, function(i) {
    if (verbose) 
      setTxtProgressBar(pb, i/l.GeneSet)
    hyperGeoTest(collectionOfGeneSets[i], universe, hits)
  }))
  if (verbose) 
    close(pb)
  if (length(results) > 0) {
    adjPvals_over <- p.adjust(results[, "Pvalue_over"], method = pAdjustMethod)
    results <- cbind(results, adjPvals_over)
    colnames(results)[ncol(results)] <- "Adjusted_Pvalue_over"
    
    adjPvals_under <- p.adjust(results[, "Pvalue_under"], method = pAdjustMethod)
    results <- cbind(results, adjPvals_under)
    colnames(results)[ncol(results)] <- "Adjusted_Pvalue_under"
    
    adjPvals_twosided <- p.adjust(results[, "Pvalue_twosided"], method = pAdjustMethod)
    results <- cbind(results, adjPvals_twosided)
    colnames(results)[ncol(results)] <- "Adjusted_Pvalue_twosided"
    
    results <- results[order(results[, "Adjusted_Pvalue_twosided"]), 
                       , drop = FALSE]
  }
  else {
    results <- matrix(, nrow = 0, ncol = 7)
    colnames(results) <- c("Universe Size", "Gene Set Size", 
                           "Total Hits", "Expected Hits", "Observed Hits", "Pvalue_over","Pvalue_under","Pvalue_twosided","Adjusted_Pvalue_over","Adjusted_Pvalue_under","Adjusted_Pvalue_twosided")
  }
  return(results)
}
