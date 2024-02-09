kinship<-read.csv(file="plink2.rel",sep="",header=F)
id<-read.csv(file="plink2.rel.id",sep="")
colnames(kinship)=id$IID
rownames(kinship)=id$IID
write.table(kinship,file="kinship.txt",sep="\t",quote=F)
q()
