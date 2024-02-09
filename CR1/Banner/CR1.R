library(ggpubr)
library(dplyr)

CR1_geno<-read.csv(file="CR1_geno.csv")
cov<-read.csv(file="covariates.csv")
CR1_geno<-left_join(CR1_geno,cov, by=c("sampleID"="Row.names"))

m1<-lm(CR1~rs679515_dosage+updatedDX + race + gender + expired_age + PMI,data=CR1_geno)
capture.output(summary(m1), file = "rs679515.output.txt")
m1<-lm(CR1~rs679515_dosage*updatedDX + race + gender + expired_age + PMI,data=CR1_geno)
capture.output(summary(m1), file = "rs679515.output.txt", append = T)

m2<-lm(CR1~rs9429780_dosage+updatedDX + race + gender + expired_age + PMI,data=CR1_geno)
capture.output(summary(m2), file = "rs9429780.output.txt", append = T)          
m2<-lm(CR1~rs9429780_dosage*updatedDX + race + gender + expired_age + PMI,data=CR1_geno)
capture.output(summary(m2), file = "rs9429780.output.txt", append = T)          

m3<-lm(CR1~rs11118328_dosage+updatedDX + race + gender + expired_age + PMI,data=CR1_geno)
capture.output(summary(m3), file = "rs11118328.output.txt", append =F)
m3<-lm(CR1~rs11118328_dosage*updatedDX + race + gender + expired_age + PMI,data=CR1_geno)
capture.output(summary(m3), file = "rs11118328.output.txt", append =T)

m_all<-lm(CR1~rs11118328_dosage+rs9429780_dosage+rs679515_dosage+updatedDX + race + gender + expired_age + PMI,data=CR1_geno)
capture.output(summary(m_all), file = "all.output.txt"

CR1_geno$rs679515_genotype=CR1_geno$rs679515_dosage
CR1_geno$rs11118328_genotype=CR1_geno$rs11118328_dosage
CR1_geno$rs9429780_genotype=CR1_geno$rs9429780_dosage

CR1_geno$rs679515_genotype=gsub("2","CC",CR1_geno$rs679515_genotype)
CR1_geno$rs679515_genotype=gsub("1","CT",CR1_geno$rs679515_genotype)
CR1_geno$rs679515_genotype=gsub("0","TT",CR1_geno$rs679515_genotype)

CR1_geno$rs11118328_genotype=gsub("0","CC",CR1_geno$rs11118328_genotype)
CR1_geno$rs11118328_genotype=gsub("1","CT",CR1_geno$rs11118328_genotype)
CR1_geno$rs11118328_genotype=gsub("2","TT",CR1_geno$rs11118328_genotype)

CR1_geno$rs9429780_genotype=gsub("2","CC",CR1_geno$rs9429780_genotype)
CR1_geno$rs9429780_genotype=gsub("1","CG",CR1_geno$rs9429780_genotype)
CR1_geno$rs9429780_genotype=gsub("0","GG",CR1_geno$rs9429780_genotype)

CR1_geno$rs9429780_genotype=factor(CR1_geno$rs9429780_genotype,levels=c("GG","CG","CC"))
CR1_geno$rs11118328_genotype=factor(CR1_geno$rs11118328_genotype,levels=c("CC","CT","TT"))
CR1_geno$rs679515__genotype=factor(CR1_geno$rs679515_genotype,levels=c("TT","CT","CC"))

p1<-
ggboxplot(CR1_geno,x="rs9429780_genotype",y="CR1",color="DX",palette = "jco",add = "jitter")
p2<-
ggboxplot(CR1_geno,x="rs679515_genotype",y="CR1",color="DX",palette = "jco",add = "jitter")
p3<-                                                                           
ggboxplot(CR1_geno,x="rs11118328_genotype",y="CR1",color="DX",palette = "jco",add = "jitter")

p_all_dx<-ggarrange(p1,p2,p3,nrow=1,common.legend=T)
ggsave(file="CR1_all3snps_dx.pdf",p_all_dx,w=12,h=4)

q()
