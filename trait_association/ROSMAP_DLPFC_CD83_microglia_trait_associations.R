
library(dplyr)
library(pbapply)
library(data.table)
library(plyr)
library(broom)

Mathys_DA <- read.table(file = "Mathys_DA_path.csv", sep = ",", header = TRUE)
Mathys_DA$CD83_status <- ifelse(test = Mathys_DA$DA4 > 0, yes = "CD83_pos", no = "CD83_neg")

Fujita_DA <- read.table(file = "Fujita_DA_path.csv", sep = ",", header = TRUE)
Fujita_DA$CD83_status <- ifelse(test = Fujita_DA$DA_A > 0, yes = "CD83_pos", no = "CD83_neg")

shortDat <- read.csv(file = "ROSMAP_projid_annotations.csv")

Mathys_DA <- data.frame(Mathys_DA, shortDat[match(Mathys_DA$projid, shortDat$projid),])
Fujita_DA <- data.frame(Fujita_DA, shortDat[match(Fujita_DA$projid, shortDat$projid),])

Mathys_DA_AD <- Mathys_DA[which(Mathys_DA$DX == "AD"),]
Fujita_DA_AD <- Fujita_DA[which(Fujita_DA$DX == "AD"),]

Mathys_DA_CN <- Mathys_DA[which(Mathys_DA$DX == "CN"),]
Fujita_DA_CN <- Fujita_DA[which(Fujita_DA$updatedDX == "CN"),]

cohort2clinicopath_ks <- function(cohort_df, trait_col) {
  try(
    ks.test(cohort_df[[trait_col]][which(cohort_df$CD83_status == "CD83_pos")], cohort_df[[trait_col]][which(cohort_df$CD83_status == "CD83_neg")]), 
    silent = TRUE)
}
  
cohort2clinicopath_tTest <- function(cohort_df, trait_col) {
  try(
    t.test(log2(cohort_df[[trait_col]][which(cohort_df$CD83_status == "CD83_pos")]+0.5), log2(cohort_df[[trait_col]][which(cohort_df$CD83_status == "CD83_neg")]+0.5)), 
    silent = TRUE)
}

cohort2clinicopath_GLM <- function(cohort_df, trait_col) {
  cohort_df$CD83_binary <- ifelse(cohort_df$CD83_status == "CD83_pos", yes = 1, no = 0)
  form_i <- as.formula(object = paste("CD83_binary ~ ", trait_col, " + age_death + msex + DX", sep = ""))
  tryCatch(expr = summary(glm(formula = form_i, data = cohort_df, family = binomial))$coefficients, error = function(e) {NA})
}

runMultiTest <- function(AD_cohort_metadata, Full_cohort_metadata, ks_phenotypes, tTest_phenotypes, GLM_phenotypes) {
  
  #KS --------------------------------------------------------------
  
  DA_AD_ks <-
    pblapply(ks_phenotypes, function(phen_i)  {
      cohort2clinicopath_ks(cohort_df = AD_cohort_metadata, trait_col = phen_i)
    })
  names(DA_AD_ks) <- ks_phenotypes
  
  DA_AD_ks_df <- ldply(
    lapply(DA_AD_ks, function(x) {
      as.data.frame(tidy(x))
    }), rbind, .id = "Trait")
  
  DA_AD_ks_df$FDR_BH <- p.adjust(DA_AD_ks_df$p.value, method = "BH")
  
  
  # T-test of log2 transformed NP values  -----------------------------------
  
  
  DA_AD_tTest <-
    pblapply(tTest_phenotypes, function(phen_i)  {
      cohort2clinicopath_tTest(cohort_df = AD_cohort_metadata, trait_col = phen_i)
    })
  names(DA_AD_tTest) <- tTest_phenotypes
  
  DA_AD_tTest_df <- ldply(
    lapply(DA_AD_tTest, function(x) {
      as.data.frame(tidy(x))
    }), rbind, .id = "Trait")
  
  DA_AD_tTest_df$FDR_BH <- p.adjust(DA_AD_tTest_df$p.value, method = "BH")
  
  # GLM ---------------------------------------------------------------------

  DA_GLM <-
    pblapply(GLM_phenotypes, function(phen_i)  {
      cohort2clinicopath_GLM(cohort_df = Full_cohort_metadata, trait_col = phen_i)
    })
  names(DA_GLM) <- GLM_phenotypes
  
  DA_GLM_df <- ldply(lapply(DA_GLM, function(y) data.frame(Term = rownames(y), y)), rbind, .id = "Trait")
  DA_GLM_df <- DA_GLM_df[which(sapply(1:nrow(DA_GLM_df), function(x) length(unique(unlist(DA_GLM_df[x,c("Trait", "Term")])))) == 1),-2]
  DA_GLM_df$FDR_BH <- p.adjust(DA_GLM_df$Pr...z.., method = "BH")
  
  list(ks_df = DA_AD_ks_df, tTes_df = DA_AD_tTest_df, GLM_df = DA_GLM_df)
  
}

# Fujita AD -----------------------------------------------------------------

ks_traits_of_interest_i <- 
  c("age_death", "educ", "gpath", "braaksc", "plaq_d", "plaq_n", "dlbdx", "nft", "caa_4gp", "cogn_global", "pmi",   "amyloid",   "tangles")

t_test_traits_of_interest_i <- 
  c( "gpath", "braaksc", "plaq_d", "plaq_n", "dlbdx", "nft", "caa_4gp", "cogn_global", "amyloid",   "tangles")

GLM_traits_of_interest_i <- 
  c("gpath", "braaksc", "plaq_d", "plaq_n", "dlbdx", "nft", "caa_4gp", "cogn_global", "pmi", "amyloid", "tangles")

Fujita_AD_summary_list <- runMultiTest(AD_cohort_metadata = Fujita_DA_AD, Full_cohort_metadata = Fujita_DA, ks_phenotypes = ks_traits_of_interest_i, tTest_phenotypes = t_test_traits_of_interest_i, GLM_phenotypes = GLM_traits_of_interest_i)

# Mathys AD (remove samples already included in Fujita) -----------------------------------------------------------------

ks_traits_of_interest_i <-
  c("age_death", "educ", "gpath", "braaksc", "plaq_d", "plaq_n", "dlbdx", "nft", "caa_4gp",  "pmi",   "amyloid", "tangles") #slightly different traits vaailble across cohorts

t_test_traits_of_interest_i <-
  c( "gpath", "braaksc", "plaq_d", "plaq_n", "dlbdx", "nft", "caa_4gp", "amyloid",   "tangles")

GLM_traits_of_interest_i <-
  c("gpath", "braaksc", "plaq_d", "plaq_n", "dlbdx", "nft", "caa_4gp", "pmi", "amyloid", "tangles")

MathysNotFujita_AD <- runMultiTest(AD_cohort_metadata = Mathys_DA_AD[match(setdiff(Mathys_DA_AD$projid, Fujita_DA_AD$projid),Mathys_DA_AD$projid),], Full_cohort_metadata = Mathys_DA[match(setdiff(Mathys_DA$projid, Fujita_DA$projid),Mathys_DA$projid),], ks_phenotypes = ks_traits_of_interest_i, tTest_phenotypes = t_test_traits_of_interest_i, GLM_phenotypes = GLM_traits_of_interest_i)


# Combined unique set of Fujita & Mathys ----------------------------------------------------------------


ks_traits_of_interest_i <-
  c("age_death", "educ", "gpath", "braaksc", "plaq_d", "plaq_n", "dlbdx", "nft", "caa_4gp",  "pmi",   "amyloid", "tangles")#"cogn_global",

t_test_traits_of_interest_i <-
  c( "gpath", "braaksc", "plaq_d", "plaq_n", "dlbdx", "nft", "caa_4gp", "amyloid",   "tangles")

GLM_traits_of_interest_i <-
  c("gpath", "braaksc", "plaq_d", "plaq_n", "dlbdx", "nft", "caa_4gp", "pmi", "amyloid", "tangles")

sharedCols <- c("projid", "CD83_status", "msex", "DX", "age_death", "educ", "gpath", "braaksc", "plaq_d", "plaq_n", "dlbdx", "nft", "caa_4gp", "pmi",   "amyloid", "tangles")

Combined_AD <- rbind(Fujita_DA_AD[,sharedCols],Mathys_DA_AD[match(setdiff(Mathys_DA_AD$projid, Fujita_DA_AD$projid),Mathys_DA_AD$projid),sharedCols])
Combined_AD_CN <- rbind(Fujita_DA[,sharedCols],Mathys_DA[match(setdiff(Mathys_DA$projid, Fujita_DA$projid),Mathys_DA$projid),sharedCols])

Mathys_u_Fujita_AD <- runMultiTest(AD_cohort_metadata = Combined_AD, 
                                   Full_cohort_metadata = Combined_AD_CN, 
                                   ks_phenotypes = ks_traits_of_interest_i, tTest_phenotypes = t_test_traits_of_interest_i, GLM_phenotypes = GLM_traits_of_interest_i)

# Output key results ------------------------------------------------------

Fujita_AD_df <- rbind(Fujita_AD_summary_list$ks_df, Fujita_AD_summary_list$tTes_df[,c("Trait", "statistic", "p.value","FDR_BH", "method", "alternative")])
Mathys_not_Fujita_AD <- rbind(MathysNotFujita_AD$ks_df, MathysNotFujita_AD$tTes_df[,c("Trait", "statistic", "p.value","FDR_BH", "method", "alternative")])
Mathys_and_Fujita_AD <- rbind(Mathys_u_Fujita_AD$ks_df, Mathys_u_Fujita_AD$tTes_df[,c("Trait", "statistic", "p.value","FDR_BH", "method", "alternative")])

CD83_AD_vs_non_CD83_AD_summary <- rbind(data.frame(Cohort = "Fujita_AD", Fujita_AD_df),
                                        data.frame(Cohort = "Mathys_not_Fujita_AD", Mathys_not_Fujita_AD),
                                        data.frame(Cohort = "Mathys_and_Fujita_AD", Mathys_and_Fujita_AD)
                                        )
CD83_AD_vs_non_CD83_AD_summary <- data.frame(Comparison = "CD83_AD_vs_nonCD83_AD", CD83_AD_vs_non_CD83_AD_summary)



CD83_AD_and_CN_GLM_summary <- rbind(data.frame(Cohort = "Fujita_AD", Fujita_AD_summary_list$GLM_df),
                                        data.frame(Cohort = "Mathys_not_Fujita_AD", MathysNotFujita_AD$GLM_df),
                                        data.frame(Cohort = "Mathys_and_Fujita_AD", Mathys_u_Fujita_AD$GLM_df)
)
CD83_AD_and_CN_GLM_summary <- data.frame(Comparison = "CD83_status_binary ~ Trait + Age_death + Sex + AD_diagnosis", CD83_AD_and_CN_GLM_summary)

o_i <- list(ROSMAP_CD83_AD_vs_non_CD83_AD = CD83_AD_vs_non_CD83_AD_summary,
            ROSMAP_CD83_AD_and_CN_GLM = CD83_AD_and_CN_GLM_summary
            )

openxlsx::write.xlsx(file = "results/ROSMAP_Fujita_Mathys_CD83_microglia_trait_associations.xlsx", x = o_i, colNames =TRUE, rowNames = TRUE)
