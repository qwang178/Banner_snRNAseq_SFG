library(tidyverse)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

all<-read.csv(file="../data_S.NOMIS1.csv")
all$cluster<-paste(all$cluster_annotated,all$da0,sep="_")

AD<-subset(all, da0=="0" & updatedDX=="AD")
CN<-subset(all, da0=="0" & updatedDX=="Control")
da<-subset(all, da0!="0")                                   

ave_da<-                                                                        
da %>%                                                                          
  group_by(cluster) %>%                                               
  summarise_at(vars(args[1]), list(name = mean), na.rm=T)                            
ave_da<-column_to_rownames(ave_da, var = "cluster") 

ave_AD<-                                                                        
AD %>%                                                                          
  group_by(cluster) %>%                                               
  summarise_at(vars(args[1]), list(name = mean), na.rm=T)                                 
ave_AD<-column_to_rownames(ave_AD, var = "cluster") 

ave_CN<-                                                                        
CN %>%                                                                          
  group_by(cluster) %>%                                               
  summarise_at(vars(args[1]), list(name = mean), na.rm=T)                                 
ave_CN<-column_to_rownames(ave_CN, var = "cluster") 

n_perm=10000
df<-data.frame(matrix(ncol = 2, nrow = 0))                  

for (p in rownames(ave_AD)) {
    
    print(p)
 
    AD_sub<-subset(AD,cluster==p)
    permutations <- cbind(sapply(1:n_perm, function(x) sample(1:dim(all)[1], dim(AD_sub)[1])))
    ave_AD_matrix<-cbind(sapply(1:n_perm, function(x) all[permutations[,x],args[1]]))
    ave_AD_mean<-rowMeans(t(ave_AD_matrix), na.rm=T)
    mean_AD<-mean(ave_AD_mean, na.rm=T)                                                        
    sd_AD<-sd(ave_AD_mean, na.rm=T)
    z_score_AD=(ave_AD[p,"name"]-mean_AD)/sd_AD                                                               
    df[nrow(df)+1,]=c(paste0(p,"_AD"),z_score_AD)

}

for (p in rownames(ave_CN)) {                                                   
                                                                                
    print(p)                                                                    

    CN_sub<-subset(CN,cluster==p)                                               
    permutations <- cbind(sapply(1:n_perm, function(x) sample(1:dim(all)[1], dim(CN_sub)[1])))
    ave_CN_matrix<-cbind(sapply(1:n_perm, function(x) all[permutations[,x],args[1]])) 
    ave_CN_mean<-rowMeans(t(ave_CN_matrix), na.rm=T)                                               
    mean_CN<-mean(ave_CN_mean, na.rm=T)                                                  
    sd_CN<-sd(ave_CN_mean, na.rm=T)                                                      
    z_score_CN=(ave_CN[p,"name"]-mean_CN)/sd_CN                                           
    df[nrow(df)+1,]=c(paste0(p,"_CN"),z_score_CN)

}                                                                               
                                                                                          
for (p in rownames(ave_da)) {                                                   

    print(p)

    da_sub<-subset(da,cluster==p)                                               
    permutations <- cbind(sapply(1:n_perm, function(x) sample(1:dim(all)[1], dim(da_sub)[1])))
    ave_da_matrix<-cbind(sapply(1:n_perm, function(x) all[permutations[,x],args[1]])) 
    ave_da_mean<-rowMeans(t(ave_da_matrix), na.rm=T)                                               
    mean_da<-mean(ave_da_mean, na.rm=T)                                                  
    sd_da<-sd(ave_da_mean, na.rm=T)                                                      
    z_score_da=(ave_da[p,"name"]-mean_da)/sd_da                                           
    df[nrow(df)+1,]=c(paste0(p,"_da"),z_score_da)
}

colnames(df)=c("cluster","Z_score")
rownames(df)=df$cluster
df<-df[-1]
write.csv(file="z_score.csv",df,quote=F)


