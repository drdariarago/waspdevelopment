## Script to detect Differential expression versus Splicing
## Workflow: load gene, exon and splicing DE, then apply the following rules:
## For DE: gene must be DE, all exons must be DE, all exon must have same sign of coefficient
## For Splicing: >1 exon must be DE, >1 exon must have splicing
## LIMMA alternative to percent spliced in (see Drosophila developmental transcriptome): Number of times the exon has positive value over background

## Load libraries and data
library(limma)
library(plyr)
library(reshape2)

load(file="./Output/BaseStats_fit2.RData")
load(file="./Output/Exon_BaseStats_FDRtool.RData")
load(file="./Output/Splicing_BaseStats_fit2_fdrtool.RData")

## Create table with 1 gene and factor per row and 1 logical condition per column
gene_stage<-expand.grid(colnames(fit2$coefficients),row.names(fit2$coefficients))
names(gene_stage)<-c("Factor","geneID")
data_splicingDE<-gene_stage[,c("geneID","Factor")]
row.names(data_splicingDE)<-paste(gene_stage$geneID, gene_stage$Factor, sep="_")

## Set fdr threshold to be used in all further steps
fdr_threshold<-0.00001

## Calculate gene DE based on fdr_threshold
data_splicingDE$GeneDE<-c(fit2$fdr<fdr_threshold)

## Calculate significant DE for all and any Exons in gene per condition
exonDE<-melt(fit2_exon$fdr)
exonDE$geneID<-as.factor(substr(exonDE$Var1,0,14))
names(exonDE)<-c("transcriptID","Factor","fdr","geneID")
exonDE_2<-ddply(exonDE, .variables=.(geneID, Factor), summarize, anyExonDE=any(fdr<fdr_threshold), allExonDE=all(fdr<fdr_threshold))
data_splicingDE<-merge(data_splicingDE, exonDE_2, by=c("geneID","Factor"))

## Calculate which genes have exons with same sign of coefficient in each condition
exonCoef<-fit2_exon$coefficients
row.names(exonCoef)<-fit2_exon$genes
exonCoef<-melt(exonCoef)
exonCoef$geneID<-as.factor(substr(exonCoef$Var1,0,14))
names(exonCoef)<-c("transcriptID","Factor","Coefficient","geneID")
exonCoef_2<-ddply(exonCoef, .variables=.(geneID, Factor), summarize, sameExonCoef=all(sign(Coefficient)==sign(Coefficient[1])))
data_splicingDE<-merge(data_splicingDE, exonCoef_2, by=c("geneID", "Factor"))

## Calculate which genes have alternatively spliced exons
splicingDE<-fit2_splicing$fdr
splicingDE<-melt(splicingDE)
splicingDE$geneID<-as.factor(substr(splicingDE$Var1, 0,14))
names(splicingDE)<-c("transcriptID","Factor","fdr","geneID")
splicingDE_2<-ddply(splicingDE, .(geneID, Factor), summarize, anyExonSpliced=any(fdr<fdr_threshold))
data_splicingDE<-merge(data_splicingDE, splicingDE_2, by=c("geneID","Factor"))

## Clean workspace and save data
rm(ls()[-grep("data_splicingDE", ls())])
save(data_splicingDE, file="./Output/SplicingDE.RData")


## Define which genes are DE
detect_DE<-function(data){data$GeneDE==T&data$allExonDE==T&data$sameExonCoef==T}
## And apply to our dataframe
data_splicingDE$Gene_DE_2<-detect_DE(data_splicingDE)

## Define which genes are alternatively spliced
detect_splicing<-function(data){data$anyExonDE==T&data$anyExonSpliced==T}
## And apply to our dataframe
data_splicingDE$Spliced<-detect_splicing(data_splicingDE)

## Save results to csv and R object
write.csv(data_splicingDE, file="./Output/SplicingDE.csv")
save(data_splicingDE, file="./Output/SplicingDE.RData")