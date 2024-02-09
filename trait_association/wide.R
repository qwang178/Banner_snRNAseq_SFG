 library(tidyverse)
 all<-read.csv(file="combined.count.csv",header=F)
 colnames(all)=c("count","cluster","cellType","feature")
 all$combined<-paste(all$cellType,all$cluster,sep="_")
 wide = all %>% 
 spread(feature, count)
 wide[is.na(wide)]<-0
 rownames(wide)=wide$combined
 wide<-wide[-3]
 wide$sum=wide[,3]+wide[,4]
 write.csv(file="all.count.csv",wide,quote=F)

