## Second script for splicing detection
## Ideal workflow: cbind exon and gene data, annotate gene and exon "samples", specify gene-exon as main contrast with gene as intercept and stage/sex as further nested factors.
# run limmma and FDR
## Check for intercept coefficients: if ~1, then subtraction is an acceptable simplification

#### Load and merge exon and gene datasets

## Load exon data
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON

# filter only rows with samples (exclude means)
nasdevexon<-nasonia_devtesova.exon.r99[,grep("[.][123]",names(nasonia_devtesova.exon.r99))]
# Exclude gonads 
nasdevexon<-nasdevexon[,-grep("(ovaries|testes)", names(nasdevexon))]
# Append exon ID to colnames and convert to data.frame
colnames(nasdevexon)<-paste(colnames(nasdevexon), "exon", sep="_")
nasdevexon<-as.data.frame(nasdevexon)

## Load gene data
nasonia_devtesova.gene.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.gene.r99.score")
row.names(nasonia_devtesova.gene.r99)<-nasonia_devtesova.gene.r99$GENE

# filter only rows with samples (exclude means)
nasdevgene<-nasonia_devtesova.gene.r99[,grep("[.][123]",names(nasonia_devtesova.gene.r99))]
# Exclude gonads 
nasdevgene<-nasdevgene[,-grep("(ovaries|testes)", names(nasdevgene))]
# Append gene ID to colnames and convert to data.frame
colnames(nasdevgene)<-paste(colnames(nasdevgene), "gene", sep="_")
nasdevgene<-as.data.frame(nasdevgene)

## merge datasets by gene ID (defined as all characters before t in exon name)
nasdevexon$geneID<-unlist(regmatches(row.names(nasdevexon), gregexpr("^[^t]*",row.names(nasdevexon))))
nasdevgene$geneID<-row.names(nasdevgene)
nasdevdata<-merge(nasdevexon, nasdevgene, by="geneID", all.x=T)
row.names(nasdevdata)<-row.names(nasdevexon)
nasdevdata<-nasdevdata[,-grep("geneID", names(nasdevdata))]
  
# Sava dataset before analyses
write.csv(nasdevdata, file="./Output/Splicing2/nasdevsplicing_intercept_data.csv")


###### LIMMA analyses
library(limma)

## Reformat as list
nasdevlist_intercept<-list(E=as.matrix(nasdevdata), genes=as.character(row.names(nasdevdata)), targets=names(nasdevdata))
## BEAR does not add gene names and sample names to rows and column for some reason, adding them manually
row.names(nasdevlist_intercept$E)<-nasdevlist_intercept$genes
colnames(nasdevlist_intercept$E)<-nasdevlist_intercept$targets
# Converting to LIMMA EList format
nasdevEList_intercept<-new("EList", nasdevlist_intercept)
# remove list (disable for diagnostics)
rm(list=c("nasdevlist_intercept"))

## Create model matrix
stage<-factor(as.character(regmatches(nasdevEList_intercept$targets, gregexpr("^[[:alnum:]]+",nasdevEList_intercept$targets))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
stage_ordered<-ordered(stage, levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
# male<-grepl("[^[:alnum:]]male|testes", nasdevEList$targets)
# male<-factor(ifelse(male==T, "Male", "Female"),levels=c("Male", "Female"))
male<-factor(ifelse(grepl("[^[:alnum:]]male|testes", nasdevEList_intercept$targets)==T, "Male", "Female"),levels=c("Male", "Female"))
exon<-factor(ifelse(grepl("exon", nasdevEList_intercept$targets)==T, "exon","gene"), levels=c("gene","exon"))

# options(contrasts=c("contr.treatment", "contr.poly"))
library(car)
## calculate gene-based intercept for exon values, then apply nested stage/sex contrasts
design<-model.matrix(~0+stage/male+exon, contrasts.arg=list(exon="contr.treatment",stage="contr.Sum", male="contr.Sum"))
# check if data matrix is correct
as.data.frame(cbind(as.character(stage), as.character(male), design))

## fit models
fit_intercept<-lmFit(nasdevEList_intercept, design)
fit2_intercept<-eBayes(fit_intercept)




# Alternative version, fit new model to residuals extracted from gene-exon contrasts
design<-model.matrix(~0+exon, contrasts.arg=list(exon="contr.treatment"))


## Extract residuals for exons
nasdevsplicing<-residuals.MArrayLM(fit2_intercept, nasdevEList_intercept)
nasdevsplicing<-nasdevsplicing[,grep("exon", colnames(nasdevsplicing))]
  
#### Re-fit models to residuals after removal of intercept

## Reformat as list
nasdevlist_splicing<-list(E=as.matrix(nasdevsplicing), genes=as.character(row.names(nasdevsplicing)), targets=colnames(nasdevsplicing))
## BEAR does not add gene names and sample names to rows and column for some reason, adding them manually
row.names(nasdevlist_splicing$E)<-nasdevlist_splicing$genes
colnames(nasdevlist_splicing$E)<-nasdevlist_splicing$targets
# Converting to LIMMA EList format
nasdevEList_splicing<-new("EList", nasdevlist_splicing)
# remove list (disable for diagnostics)
rm(list=c("nasdevlist_splicing"))

## Create model matrix
stage<-factor(as.character(regmatches(nasdevEList_splicing$targets, gregexpr("^[[:alnum:]]+",nasdevEList_splicing$targets))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
stage_ordered<-ordered(stage, levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
male<-factor(ifelse(grepl("[^[:alnum:]]male|testes", nasdevEList_splicing$targets)==T, "Male", "Female"),levels=c("Male", "Female"))

# options(contrasts=c("contr.treatment", "contr.poly"))
library(car)
design<-model.matrix(~0+stage/male, contrasts.arg=list(stage="contr.Sum",male="contr.Sum"))
# check if data matrix is correct
as.data.frame(cbind(as.character(stage), as.character(male), design))

## fit models
fit_splicing<-lmFit(nasdevEList_splicing, design)
fit2_splicing<-eBayes(fit2_splicing)

## apply FDR 

