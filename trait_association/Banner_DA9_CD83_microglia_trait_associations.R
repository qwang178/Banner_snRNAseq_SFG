
library(dplyr)
library(pbapply)
library(data.table)
library(plyr)
library(openxlsx)
library(broom)


Mic9_AD_df <- read.table(file = "DA9_meta.csv", sep = ",", header = TRUE) #Data frame of Sample ID, DA9/CD83 status and comorbidity, demographic & clinicopathological covariates
Mic9_AD_df <- Mic9_df[which(Mic9_df$updatedDX == "AD"),]

Mic9_AD_CD83_df <- Mic9_AD_df[which(Mic9_AD_df$Mic9 > 0),]
Mic9_AD_nonCD83_df <- Mic9_AD_df[which(Mic9_AD_df$Mic9 == 0),]

ks_traits <- c(c("expired_age","gender.1", "PMI", "PlaqueF","PlaqueT","PlaqueP","PlaqueH","PlaqueE","PlaqueTotal"),
            c("TangleF","TangleT","TangleP","TangleE","TangleH","TangleTotal"),
            c("obt", "brain_stem_ix_x", "brain_stem_lc", "bf_amygdala", "bf_nbm", "brain_stem_sn", "bf_trans", "bf_cing", "nctx_temporal", "nctx_frontal", "sum_lb_density", "nctx_parietal", "dementia_years")
)

ksRes <-
pblapply(ks_traits, function(q_i) {
  try(
    ks.test(Mic9_AD_CD83_df[[q_i]], Mic9_AD_nonCD83_df[[q_i]]), silent = TRUE)
})
names(ksRes) <- ks_traits

ksRes <- ksRes[sapply(ksRes, class) != "try-error"]


ksRes_df <- pblapply(ksRes, function(x) {
  as.data.frame(tidy(x))
})
ksRes_df <- ldply(ksRes_df, rbind, .id = "Trait")

# T-test of log2 transformed NP values  -----------------------------------
t_test_traits <-  c("PlaqueF","PlaqueT","PlaqueP","PlaqueH","PlaqueE","PlaqueTotal","TangleF","TangleT","TangleP","TangleE","TangleH","TangleTotal")

tRes <-
  pblapply(t_test_traits, function(q_i) {
    try(
      t.test(log2(Mic9_AD_CD83_df[[q_i]]+0.5), log2(Mic9_AD_nonCD83_df[[q_i]]+0.5), silent = TRUE)
    )
  })
names(tRes) <- t_test_traits

tRes <- tRes[sapply(tRes, class) != "try-error"]

tRes_df <- pblapply(tRes, function(x) {
  as.data.frame(tidy(x))
})
tRes_df <- ldply(tRes_df, rbind, .id = "Trait")

Banner_CD83_AD_NP_assoc <- rbind(ksRes_df, tRes_df[,c("Trait", "statistic", "p.value", "method", "alternative")])

# Comorbidity Dx -----------------------------------------------------------------

queryFields <- c("Circ.Wilis.Score.3", "Atrial.Fib", "Pneumonia","Any.acute.brain.inf.or.ischemia", 
                 "COPD","AnyApoE4", "BMI...30", "Cardiomegaly", "High.Heart.Weight", 
                 "Hypertension", "Diabetes", "Smoking", "Summer", "LBS", "PARKINSONISM.NOS",
                 "diverticulosis_BodyPathDXSummary", "hepatic_congestion_BodyPathDXSummary", "hepatic_atrophy_BodyPathDXSummary", "coronary_stenosis_BodyPathDXSummary", "Cardiomegaly_BodyPathDXSummary", "atherosclerosis_BodyPathDXSummary")

fRes <- 
  pblapply(queryFields, function(trait_i) {
    fisher.test(table(Mic9_AD_df[,c(trait_i, "CD83_microglia")]))
  })
names(fRes) <- queryFields

fisherRes_df <- pblapply(fRes, function(x) {
  as.data.frame(tidy(x))
})
FisherRes_df <- ldply(fisherRes_df, rbind, .id = "Trait")
Banner_CD83_AD_Comorb_assoc <- FisherRes_df[order(FisherRes_df$p.value, decreasing = FALSE),]

# GLM based approach across all samples -----------------------------------
glm_traits <- c(c("PMI", "PlaqueF","PlaqueT","PlaqueP","PlaqueH","PlaqueE","PlaqueTotal"),
               c("TangleF","TangleT","TangleP","TangleE","TangleH","TangleTotal"),
               c("obt", "brain_stem_ix_x", "brain_stem_lc", "bf_amygdala", "bf_nbm", "brain_stem_sn", "bf_trans", "bf_cing", "nctx_temporal", "nctx_frontal", "sum_lb_density", "nctx_parietal", "dementia_years")
)

Mic9_df$CD83_status <- ifelse(test = Mic9_df$Mic9 > 0, yes = "CD83_pos", no = "CD83_neg")


glmList <-
  pblapply(1:length(glm_traits), function(prot_i) {
    Mic9_df$CD83_binary <- ifelse(Mic9_df$CD83_status == "CD83_pos", yes = 1, no = 0)
    form_i <- as.formula(object = paste("CD83_binary ~ ", glm_traits[prot_i], " + expired_age + Sex + updatedDX", sep = ""))
    tryCatch(expr = summary(glm(formula = form_i, data = Mic9_df, family = binomial))$coefficients, error = function(e) {NA})
  })
names(glmList) <- glm_traits

glmList <- glmList[!is.na(glmList)]

DA_GLM_df <- ldply(lapply(glmList, function(y) data.frame(Term = rownames(y), y)), rbind, .id = "Trait")
DA_GLM_df <- DA_GLM_df[which(sapply(1:nrow(DA_GLM_df), function(x) length(unique(unlist(DA_GLM_df[x,c("Trait", "Term")])))) == 1),-2]
DA_GLM_df$FDR_BH <- p.adjust(DA_GLM_df$Pr...z.., method = "BH")
Banner_CD83_All_samples_GLM <- DA_GLM_df

o_i$Banner_CD83_AD_Comorb_assoc <- Banner_CD83_AD_Comorb_assoc
o_i$Banner_CD83_AD_NP_assoc <- Banner_CD83_AD_NP_assoc
o_i$Banner_CD83_All_samples_GLM <- Banner_CD83_All_samples_GLM

openxlsx::write.xlsx(file = "results/Banner_DA9_CD83_microglia_trait_associations.xlsx", x = o_i, colNames =TRUE, rowNames = FALSE)
