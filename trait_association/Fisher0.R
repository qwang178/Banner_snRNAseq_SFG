f <- function(x, output) {
testor=rbind(c(as.numeric(x[2]),as.numeric(x[3])),c(as.numeric(x[4]),as.numeric(x[5])))
capture.output(fisher.test(testor,alternative="greater"), file = paste(x[1],"tests1.txt",sep="."))
foo=fisher.test(testor,alternative="greater")
capture.output(foo$estimate[1], file=paste(x[1],"FisherOR.txt",sep="."))
capture.output(foo$p.value, file=paste(x[1],"p.txt",sep="."))
}
dir<-read.csv(file="dir.lst",header=F)
dataList=dir$V1
for(i in dataList){
    setwd(i)
    fisher<-read.csv(file="Fisher.csv",header=F,sep="")
    apply(fisher, 1, f)
}

