library(mgcv)
library(plotrix)
library(latex2exp)
library(cfdr)
library(dplyr)

cortex<-read.csv(file="/data/breadhea/qwang178/NOMIS/scRNA/Drop8/eQTL/data/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.GRCh38.csv")
qtl<-read.csv(file="qtl_results_all.txt", sep="")
inner_joined<-inner_join(qtl, cortex, by = c("ensembl_gene_id" = "gene", "snp_id" = "GRCh38"))

p_brca=inner_joined$empirical_feature_p_value; names(p_brca)=paste(inner_joined$snpLocId,inner_joined$ensembl_gene_id,sep="_")
p_oca=inner_joined$FDR; names(p_oca)=paste(inner_joined$snpLocId,inner_joined$ensembl_gene_id,sep="_")
z_brca=-qnorm(p_brca/2); z_oca=-qnorm(p_oca/2)
w=which(is.finite(z_brca+z_oca))
p_brca=p_brca[w]; p_oca=p_oca[w]
id=as.character(inner_joined$ensembl_gene_id)
fold=as.numeric(as.factor(id))
nfold=max(fold)

chr=as.character(inner_joined$chromosome)
pos=as.character(inner_joined$snp_postion)
gene=as.character(inner_joined$ensembl_gene_id)

sub_brca=which(z_oca^2 + z_brca^2 > 4^2 )
xbrca1=matrix(-1,length(sub_brca),2004); xbrca2=xbrca1
r_brca=rank(p_brca)
totalx=0
for (xfold in 1:max(fold)) {
inx=which(fold==xfold)
sub=intersect(sub_brca,inx)
if (length(sub)>0) {
  v1=vl(p_brca,p_oca,adj=T,indices=sub,fold=inx,mode=2,nv=2000)
  v2=vl(p_brca,p_oca,adj=F,indices=sub,fold=inx,mode=2,nv=2000)
  xbrca1[match(sub,sub_brca),]=v1$x; ybrca1=v1$y
  xbrca2[match(sub,sub_brca),]=v2$x; ybrca2=v2$y
} else v1=NULL
totalx=totalx+length(sub)
print(c(xfold,totalx))
}
sub_oca=which(z_oca^2 + z_brca^2 > 4^2 )
xoca1=matrix(-1,length(sub_oca),2004); xoca2=xoca1
r_oca=rank(p_oca)
totalx=0
for (xfold in 1:max(fold)) {
inx=which(fold==xfold)
sub=intersect(sub_oca,inx)
if (length(sub)>0) {
  v1=vl(p_oca,p_brca,adj=T,indices=sub,fold=inx,mode=2,nv=2000)
  v2=vl(p_oca,p_brca,adj=F,indices=sub,fold=inx,mode=2,nv=2000)
  xoca1[match(sub,sub_oca),]=v1$x; yoca1=v1$y
  xoca2[match(sub,sub_oca),]=v2$x; yoca2=v2$y
} else v1=NULL
totalx=totalx+length(sub)
print(c(xfold,totalx))
}

pars_brca=fit.2g(p_oca[which(p_brca> 0.5)])$pars
pi0_brca=pars_brca[1]
sigma_brca=pars_brca[2]
# Parametrisation of Q|H0: OCA|BRCA
pars_oca=fit.2g(p_brca[which(p_oca> 0.5)])$pars
pi0_oca=pars_oca[1]
sigma_oca=pars_oca[2]

# Integrate over L: BRCA|OCA
i_brca1=il(xbrca1,ybrca1,pi0_null=pi0_brca,sigma_null=sigma_brca,dist="norm") # adjusted cFDR
i_brca2=il(xbrca2,ybrca2,pi0_null=pi0_brca,sigma_null=sigma_brca,dist="norm") # unadjusted cFDR

# Integrate over L: OCA|BRCA
i_oca1=il(xoca1,yoca1,pi0_null=pi0_oca,sigma_null=sigma_oca,dist="norm") # adjusted cFDR
i_oca2=il(xoca2,yoca2,pi0_null=pi0_oca,sigma_null=sigma_oca,dist="norm") # unadjusted cFDR

xv1=xbrca1; xv2=xbrca2; yv1=ybrca1; yv2=ybrca2; sub=sub_brca
p=p_brca; q=p_oca; v1=i_brca1; v2=i_brca2; pi0=pi0_brca; sigma=sigma_brca; 
k=length(sub); # number of values for which L-regions were computed

alpha=0.1
r=rank(p)
hp=which(r <= max(r[which(p/r < alpha/length(p))]))

# Benjamini-Hochberg method on v-values
vbc1=rep(1,length(p)); vbc1[sub]=v1; rvbc1=rank(vbc1)
mvx1=max(rvbc1[which(vbc1/rvbc1 < alpha/length(vbc1))]); mvxl1=which(sub==match(mvx1,rvbc1))
h1=which(rvbc1 <= mvx1)

filtered<-inner_joined[h1,]                                                     
write.csv(file="filtered.csv",filtered)                                         
save(file = “cFDR.RData”)
q()
