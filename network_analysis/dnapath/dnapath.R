library(dnapath)
library(clustifyr)
#gl <- gmt_to_list(path = "data/reactome/ReactomePathways1.gmt",cutoff=20)


source("normalized_to_dna.R")
source("dna_to_summary.R")
source("init_required_packages_and_files.R")

dir_counts <- "./data/counts/"                                                  
dir_output <- "./output/"                                                       
dir_log <- "./logs/"                                                            

gl <- gmt_to_list(path = "/data/breadhea/qwang178/ref/human/gmt/msigdb2.v7.5.1.symbols.gmt")      

fn_network_inference <- function(fn = cor) {
  function(x) {
    # We will only work with genes that have nonzero variation. 
    index <- which(apply(x, 2, sd) > 0)
    if(length(index) < 2) return(NA)
    
    # Center and scale the gene expression profile by column.
    x[, index] <- scale(x[, index])
    
    # Estimate the association network, then standardize the matrix of scores.
    scores <- fn(x[, index])
    diag(scores) <- 0
    
    # Complete the association matrix A.
    p <- ncol(x)
    if(length(index) == p) {
      A <- scores
    } else {
      A <- matrix(0, p, p)
      A[index, index] <- scores
    }
    colnames(A) <- colnames(x)
    return(A)
  }
}

counts_files_list <- paste0(dir_counts, list.files(dir_counts))
#Load and filter counts.
counts_pair <- lapply(counts_files_list, function(csv_file) {
  create_counts(csv_file) %>%
    apply_filter_percent_zero_at_most(0.80) %>%
    (function(counts) {
      cat("\t\t", nrow(counts), "genes remaining\n")
      return(counts)
    })
})
names(counts_pair) <- gsub(".csv", "", list.files(dir_counts))
#filtered_counts_name <- save_counts(counts_pair,
#                                    file_name = NULL,
#                                    dir_output = dir_output,
#                                    make_subdir = TRUE)
#filtered_counts_name <- "zeroes0.8"
#dir_filtered <- paste0(dir_output, filtered_counts_name, "/")

lp_set <- c(1,2)
network_inference = fn_network_inference(corC)
pathway_list<-gl

results_by_lp <-  dna(counts_pair, pathway_list, n_perm = 100,network_inference = network_inference, lp_set = lp_set)

for(k in 1:length(results_by_lp)) {
    results_list <- results_by_lp[[k]]
    
    results_list <- results_list[!sapply(results_list, is.null)]
    results_list <- results_list[!sapply(results_list, function(x) is.na(x$p_value_path))]
    
    save_file <- paste0("output/zeroes0.8/", "lp", lp_set[k], ".rds")
    cat("Saving dna results to", save_file, "\n")
    saveRDS(results_list, save_file)
  } # End for k in lp_set

