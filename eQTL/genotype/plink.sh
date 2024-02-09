plink2 --vcf $1.vcf.gz --make-bed --out $1 
plink2 --bfile $1 --set-missing-var-ids @:# --make-bed --out $1.annotated
plink2 --bfile $1.annotated --maf 0.05 --make-bed --out $1.maf0.05
plink2 --bfile $1.maf0.05 --hwe 1e-6 --make-bed --out $1.hwe
plink2 --bfile $1.hwe --indep-pairwise 250 50 0.4 --out $1.pruned
plink2 --bfile ../beagled/common.annotated --extract $1.pruned.prune.in --make-bed --out $1.imputed
plink2 --make-rel square --bfile $1.imputed
Rscript ../kinship.R
mv kinship.txt $1.kinship.txt
