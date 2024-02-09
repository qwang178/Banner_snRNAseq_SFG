
library(coloc)
library(dplyr)
library(pbapply)
library(data.table)
library(liftOver)
library(GenomicRanges)

#system("wget -P . https://ctg.cncr.nl/documents/p1651/PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis.txt.gz")

AD_dat <-fread(input = "PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis.txt.gz", data.table = FALSE)

eQTL_dat <-fread(input = "qtl_results_all.txt", data.table = FALSE)
CR1_eQTL_dat <- eQTL_dat[which(eQTL_dat$feature_id == "CR1"),]

CR1_GRange <- GRanges(seqnames = paste("chr", CR1_eQTL_dat$snp_chromosome, sep = ""), ranges = IRanges(start = CR1_eQTL_dat$snp_position, end = CR1_eQTL_dat$snp_position))

#Need to convert our GRCh38 eQTL coordinates to GRCh37/Hg19
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
CR1_GRange_Hg19 = liftOver(CR1_GRange, ch)
CR1_eQTL_dat$Hg19_coords <- sapply(CR1_GRange_Hg19, function(x) start(x))

CR1_coords_Hg19 <- liftOver(GRanges(seqnames = "chr1", ranges = IRanges(start = 207496147, end = 207641765)), ch)
buffer_bp <- 250000

AD_CR1_bounds <- AD_dat[which(AD_dat$base_pair_location >= (start(CR1_coords_Hg19[[1]])-buffer_bp) & AD_dat$base_pair_location <= (end(CR1_coords_Hg19[[1]])+buffer_bp)),]
AD_CR1_bounds$SNP_ID_coordinates <- paste(AD_CR1_bounds$chromosome, AD_CR1_bounds$base_pair_location, sep = ":")

CR1_eQTL_dat$SNP_ID_coordinates <- paste(CR1_eQTL_dat$snp_chromosome, CR1_eQTL_dat$Hg19_coords, sep = ":")

SNP_scope <- intersect(AD_CR1_bounds$SNP_ID_coordinates, CR1_eQTL_dat$SNP_ID_coordinates)

AD_scoped <- AD_CR1_bounds[match(SNP_scope, AD_CR1_bounds$SNP_ID_coordinates),]
CR1_eQTL_scoped <- CR1_eQTL_dat[match(SNP_scope, CR1_eQTL_dat$SNP_ID_coordinates),]

#A genome-wide association study with 1,126,563 individuals identifies new risk loci for Alzheimer’s disease
#See Supp Note: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10243600/bin/NIHMS1871222-supplement-PGC-ALZ2SupplementaryNote.pdf
#The UK Biobank (UKB; www.ukbiobank.ac.uk)5 summary statistics for 46,613 cases and 318,246 controls were obtained from Jansen et al. (2019)1. 
#The 23andMe data consists of 3807 cases and 359,839 controls. Among the controls, there were 19,638 individuals between the age of 45-60 and 340,201 individuals over 60. There were 130 cases between 45-60 and 3677 cases over the age of 60.

N_samples <- 1126563-(46613 + 318246) - (3807 + 359839)


CaseNumbers <- list(deCODE = 7002,
                    HUNT = 1156,
                    BioVU = 600,
                    DemGene = 1638,
                    STSA = 320,
                    TwinGene = 224,
                    IGAP = 21982,
                    Finngen = 1798)

# > sum(unlist(CaseNumbers))
# [1] 34720

# > 34720/N_samples
# [1] 0.08722347

AD_dataset_list <- list(pvalues = AD_scoped$p_value,
                        s=0.087,
                        N=N_samples,
                        beta = AD_scoped$beta,
                        type = "cc"
                        )

eQTL_dataset_list <- list(pvalues = CR1_eQTL_scoped$p_value,
                        N=99,
                        beta = CR1_eQTL_scoped$beta,
                        type = "quant")



AD_coloc_list <- 
coloc.abf(dataset1 = AD_dataset_list, dataset2 = eQTL_dataset_list, MAF = CR1_eQTL_scoped$maf, p1 = 1e-04, p2 = 1e-04,
          p12 = 1e-05)

# > AD_coloc_list$summary
# nsnps    PP.H0.abf    PP.H1.abf    PP.H2.abf    PP.H3.abf    PP.H4.abf 
# 1.122000e+03 3.983172e-24 1.227152e-06 5.075669e-20 1.465198e-02 9.853468e-01 

o_i <- list(Colocalization_summary = AD_coloc_list$summary, Colocalization_results = AD_coloc_list$results)
#openxlsx::write.xlsx(file = "NOMIS_CR1_Oligo_eQTL_AD_GWAS_colocalization_v1.xlsx", x = o_i, colNames =TRUE, rowNames = TRUE)
