## load raw data, analyze with LIMMA and add Fdr and fdr for developmental time series and gonads vs adults
rm(list=ls())

# data import and formatting
nasonia_devtesova.gene.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.gene.r99.score")
row.names(nasonia_devtesova.gene.r99)<-nasonia_devtesova.gene.r99$GENE
# filter only rows with samples (exclude means)
nasoniadevgene<-nasonia_devtesova.gene.r99[,grep("[.][123]",names(nasonia_devtesova.gene.r99))]
# Exclude gonads 
nasoniadevgene<-nasoniadevgene[,-grep("(ovaries|testes)", names(nasoniadevgene))]

# Sava dataset before analyses
write.csv(nasoniadevgene, file="./Output/nasoniadevgene.csv")

## LIMMA analyses
## convert into EList object (list first, then create as EList)
library(limma)
nasoniadevlist<-list(E=as.matrix(nasoniadevgene), genes=as.character(row.names(nasoniadevgene)), targets=names(nasoniadevgene))
nasdevEList<-new("EList", nasoniadevlist)
rm(list=c("nasoniadevlist","nasoniadevgene","nasonia_devtesova.gene.r99"))

# Create complete design matrix
stage<-factor(as.character(regmatches(nasdevEList$targets, gregexpr("^[[:alnum:]]+",nasdevEList$targets))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
stage_ordered<-ordered(stage, levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
# male<-grepl("[^[:alnum:]]male|testes", nasdevEList$targets)
# male<-factor(ifelse(male==T, "Male", "Female"),levels=c("Male", "Female"))
male<-factor(ifelse(grepl("[^[:alnum:]]male|testes", nasdevEList$targets)==T, "Male", "Female"),levels=c("Female", "Male"))

# options(contrasts=c("contr.treatment", "contr.poly"))
library(car)
design<-model.matrix(~0+stage/male, contrasts.arg=list(stage="contr.Sum",male="contr.treatment"))
# check if data matrix is correct
as.data.frame(cbind(nasdevEList$targets, as.character(stage), as.character(male), design))

# NOTE: we can add time (in hours or relative) to pick up genes that have a linear variation across development, or any interaction with other factors
# NOTE: using orthogonal contrasts on stage_ordered would detect if each gene can be described by n-order parameters during development (constant, linear, quadratic, cubic)
# NOTE: we could complexify the nested design adding the comparisons between adult and embryo and adult and gonads, i.e.: 0+adult/stage/male or 0+gonad/adult/stage/male

# fit models
fit<-lmFit(nasdevEList, design)
fit2<-eBayes(fit)

# We can fetch either the p-value of differential expression for each stage from 
fit2$p.value
# or as logarithm of odds
fit2$lods
# effect sizes and signs can be retrieved using
fit2$coefficients


# we now apply FDR correction to each factor
library(fdrtool)
## This version applies FDR correction to every factor separately
#FDRfits<-apply(fit2$p.value, 2, function(x){fdrtool(x, statistic="pvalue")$qval})

# This version applies FDR correction to all p-values
FDRfitall<-fdrtool(as.vector(fit2$p.value), statistic="pvalue")
Fdrfitall<-matrix(FDRfitall$qval, ncol=ncol(fit2$p.value), nrow=nrow(fit2$p.value), dimnames=list(row.names(fit2$p.value), colnames(fit2$p.value))) # add global (tail area based)
fdrfitall<-matrix(FDRfitall$lfdr, ncol=ncol(fit2$p.value), nrow=nrow(fit2$p.value), dimnames=list(row.names(fit2$p.value), colnames(fit2$p.value))) # local fdr (density based)

# Save Fdr and fdr to base model
fit2$Fdr<-Fdrfitall
fit2$fdr<-fdrfitall

# write tables with results for Fdr 1E-4 (10 false positives expected) and FDR 1E-5 (1 false positive expected)
fitE4<-apply(Fdrfitall, c(1,2), function(x){x<0.0001})
fitE5<-apply(Fdrfitall, c(1,2), function(x){x<0.00001})
row.names(fitE4)<-fit2$genes
row.names(fitE5)<-fit2$genes

# and save as csv
write.csv(fitE4, file="./Output/NasoniaDev_limma_E4_unsigned.csv")
write.csv(fitE5, file="./Output/NasoniaDev_limma_E5_unsigned.csv")
write.csv(Fdrfitall, file="./Output/NasoniaDev_limma_qval_unsigned.csv")


# Add sign of sex bias
fitE4[,6:10]<-fitE4[,6:10]*sign(fit2$coefficients[,6:10])
fitE5[,6:10]<-fitE5[,6:10]*sign(fit2$coefficients[,6:10])
fitq<-Fdrfitall
fitq[,6:10]<-fitq[,6:10]*sign(fit2$coefficients[,6:10])
# and save as csv
write.csv(fitE4, file="./Output/NasoniaDev_limma_E4_signed.csv")
write.csv(fitE5, file="./Output/NasoniaDev_limma_E5_signed.csv")
write.csv(fitq, file="./Output/NasoniaDev_limma_qval_signed.csv")

# Save models as R object
save(file="./Output/BaseStats_fit2.RData", list=c("fit2"))


## Second limma model set for testing in gonad-enriched transcripts
# data import and formatting
nasonia_devtesova.gene.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.gene.r99.score")
row.names(nasonia_devtesova.gene.r99)<-nasonia_devtesova.gene.r99$GENE
# filter only rows with samples (exclude means)
nasoniadevgene<-nasonia_devtesova.gene.r99[,grep("[.][123]",names(nasonia_devtesova.gene.r99))]
# Exclude embryos 
nasoniagonadgene<-nasoniadevgene[,grep("(ovaries|testes|adult)", names(nasoniadevgene))]

# Sava dataset before analyses
write.csv(nasoniadevgene, file="./Output/nasoniagonadgene.csv")

# Convert to EList
## LIMMA analyses
## convert into EList object (list first, then create as EList)
library(limma)
nasoniagonadgene<-list(E=as.matrix(nasoniagonadgene), genes=as.character(row.names(nasoniagonadgene)), targets=names(nasoniagonadgene))
nasgonadEList<-new("EList", nasoniagonadgene)
rm(list=c("nasoniagonadgene","nasoniadevgene","nasonia_devtesova.gene.r99"))

# Create complete design matrix
# gonads<-grepl("(ovaries|testes)", nasgonadEList$targets)
# male<-grepl("(testes|_male)", nasgonadEList$targets)

adult<-as.numeric(grepl("male",nasgonadEList$targets))
adult<-as.factor(ifelse(adult==1, "adult", "gonads"))
male<-grepl("_male|testes", nasgonadEList$targets)
male<-as.factor(ifelse(male==T, "male", "female"))
adultmale<-as.factor(paste(adult, male, sep="_"))

design<-model.matrix(~0+adultmale)
colnames(design)<-levels(adultmale)
# Check model matrix
cbind(nasgonadEList$targets, design)
# fit limma model
gonadfit<-lmFit(nasgonadEList, design)
# apply relevant contrasts (adult males vs females, gonads vs whole animal)
cont.matrix<-makeContrasts(
  Male=adult_male-adult_female,
  Testes=gonads_male-adult_male,
  Ovaries=gonads_female-adult_female,
  Gonads_female=gonads_male-gonads_female,
  levels=design)

# # apply relevant contrasts (males vs females, whole adults (M and F) vs gonads (M and F), F gonads vs M gonads)
# cont.matrix<-makeContrasts(
#   Males=adult_female-adult_male,
#   Testes=adult_male-gonads_male,
#   Ovaries=adult_female-gonads_female,
#   Gonads=gonads_female-gonads_male,
#   levels=design)


# re-fit model to new contrasts
gonadfit2<-contrasts.fit(gonadfit, cont.matrix)
gonadfit2<-eBayes(gonadfit2)

# plot results
vennDiagram(decideTests(gonadfit2, method="global"))

# store q-values
library(fdrtool)
FDRfitallgonads<-fdrtool(as.vector(gonadfit2$p.value), statistic="pvalue")
Fdrfitallgonads<-matrix(FDRfitallgonads$qval, ncol=ncol(gonadfit2$p.value), nrow=nrow(gonadfit2$p.value), dimnames=list(row.names(gonadfit2$p.value), colnames(gonadfit2$p.value))) # add global (tail area based)
fdrfitallgonads<-matrix(FDRfitallgonads$lfdr, ncol=ncol(gonadfit2$p.value), nrow=nrow(gonadfit2$p.value), dimnames=list(row.names(gonadfit2$p.value), colnames(gonadfit2$p.value))) # local fdr (density based)

# Store Fdr and fdr in model parameters
gonadfit2$Fdr<-Fdrfitallgonads
gonadfit2$fdr<-fdrfitallgonads

# write tables with results for Fdr 1E-3 (10 false positives expected) and FDR 1E-4 (1 false positive expected)
fitE3gonads<-apply(Fdrfitallgonads, c(1,2), function(x){x<0.001})
fitE4gonads<-apply(Fdrfitallgonads, c(1,2), function(x){x<0.0001})
row.names(fitE3gonads)<-gonadfit2$genes
row.names(fitE4gonads)<-gonadfit2$genes
# and save as csv
write.csv(fitE3gonads, file="./Output/NasoniaGonads_limma_E3_unsigned.csv")
write.csv(fitE4gonads, file="./Output/NasoniaGonads_limma_E4_unsigned.csv")
write.csv(Fdrfitallgonads, file="./Output/NasoniaGonads_limma_qval_unsigned.csv")

# Save models as R object
save(file="./Output/BaseStats_gonadfit2.RData", list=c("gonadfit2"))

# Save workspace
save.image(file="./Output/BaseStats_out.RData")

# Clean dataset and reload models
rm(list=ls())
load(file="./Output/BaseStats_fit2.RData")
load(file="./Output/BaseStats_gonadfit2.RData")

#### outdated code and various snippets
# # remove gonad columns and add correct contrasts
# design2<-design[,-grep("stagetestes|stageovaries", colnames(design))]
# testes<-as.numeric(stage=="testes")
# testes<-testes-1*(stage=="adult"&male==T)
# ovaries<-as.numeric(stage=="ovaries")
# ovaries<-ovaries-1*(stage=="adult"&male==F)
# design2<-cbind(design2, testes, ovaries)
# # substitute 1/0 contrasts for stages in 1/-1
# design2[1:30,1:5]<-apply(design2[1:30,1:5],c(1,2),function(x) ifelse(x==0, -1, 4))

# # 95% confidence intervals for log2fold changes (using t distribution)
# confints<-qt(0.975, fit2$df.residual+fit2$df.prior)*fit2$stdev.unscaled*sqrt(fit2$s2.post)

# 
# # the code below shows that the CI is the same for all factors in a gene, which is ??
# apply(qt(0.975, fit2$df.residual+fit2$df.prior)*fit2$stdev.unscaled*sqrt(fit2$s2.post),1,function(x){
#   y<-x[1]
#   any(x!=y)
# })
# 
# # checking from the in-built function topTable
# testframe<-data.frame()
# tempframe<-as.data.frame(1:10)
# for(i in 1:ncol(fit2$coefficients)) {
#   tempframe[,1]<-topTable(fit2, confint=T, coef=i)[,1]
#   tempframe[,2]<-topTable(fit2, confint=T, coef=i)[,2]-topTable(fit2, confint=T, coef=i)[,3]
#   tempframe[,3]<-colnames(fit2$coefficients)[i]
#   testframe<-rbind(testframe, tempframe)
# }
# testframe[order(testframe[,1]),]
# # normal within same factor, not testable across factors (no gene in toptable in both), does not match with those calculated by hand