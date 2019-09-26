# Script for detecting which exons have statistically significant alternative splicing across sex/development
# Main script identical to BaseStats applied to exon database, but with 1 modification:
# Load Basestats data, subtract gene estimates from exon scores (flattened by geneID, stage and sex, merge and subtract, then unflatten)
# run limmma and FDR

#### Load LIMMA results
load(file="./Output/BaseStats_fit2.RData")
# Load exon data
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON
# filter only rows with samples (exclude means)
nasoniadevexon<-nasonia_devtesova.exon.r99[,grep("[.][123]",names(nasonia_devtesova.exon.r99))]
# Exclude gonads 
nasoniadevexon<-nasoniadevexon[,-grep("(ovaries|testes)", names(nasoniadevexon))]
# Sava dataset before analyses
write.csv(nasoniadevexon, file="./Output/nasoniadevexon.csv")

#### Subtract estimates from samples

# Annotate samples with male, stage and replicate
testData<-nasoniadevexon
testData$exonID<-row.names(testData)
# Gene names are defined as the portion of the exon name that goes from the beginning to the last non "t" character
testData$geneID<-unlist(regmatches(testData$exonID, gregexpr("^[^t]*",testData$exonID)))
library(reshape)
testData<-melt(testData)
testData$Male<-grepl("_male",testData$variable)
testData$Stage<-factor(as.character(regmatches(testData$variable, gregexpr("^[[:alnum:]]+",testData$variable))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
testData$replicate<-factor(as.character(regmatches(testData$variable, gregexpr("[[:digit:]]$",testData$variable))))
testData<-subset(testData, select=-c(variable))
# Load and annotate coefficients wit stage and male specific
testCoefs<-fit2$coefficients[row.names(fit2$coefficients)%in%testData$geneID,]
testCoefs<-melt(testCoefs)
testCoefs$Male<-grepl("male",testCoefs$X2)
testCoefs$Stage<-factor(as.character(substring(as.character(regmatches(testCoefs$X2, gregexpr("^[[:alnum:]]*",testCoefs$X2))),6)),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
testCoefs<-subset(testCoefs, select=-c(X2))
names(testCoefs)<-c("geneID","Coef","Male","Stage")

# Subtract male coefficients to males
testData2<-merge(testData[which(testData$Male==T),], testCoefs, by=c("geneID", "Stage","Male"), all.x=T)
testData2$value<-testData2$value-testData2$Coef
testData2<-testData2[,c("exonID","geneID","value","Male","Stage","replicate")]

# Subtract stage specific coefficients
testData3<-rbind(testData2, testData[which(testData$Male==F),])
testData3<-merge(testData3, testCoefs[which(testCoefs$Male==F),-which(colnames(testCoefs)=="Male")], by=c("geneID", "Stage"), all.x=T)
testData3$value<-testData3$value-testData3$Coef
testData3<-testData3[,c("exonID","geneID","Stage","Male","replicate","value")]

# Put into exon by sample format
testData3$Male<-as.factor(ifelse(testData3$Male==T, "male","female"))
nasoniadevsplicing2<-cast(testData3, formula=exonID~Stage+Male+replicate)
row.names(nasoniadevsplicing2)<-nasoniadevsplicing2$exonID
nasoniadevsplicing2<-subset(nasoniadevsplicing2, select=-exonID)

# Save as CSV and clean workspace
write.csv(nasoniadevsplicing2, file="./Output/Nasonia_Dev_Splicing.csv")
rm(list=(grep("test", ls(), value=T)))

#### Create LIMMA object
library(limma)
# Reformat as list
nasoniadevlist_splicing<-list(E=as.matrix(nasoniadevsplicing2), genes=as.character(row.names(nasoniadevsplicing2)), targets=names(nasoniadevsplicing2))
## BEAR does not add gene names and sample names to rows and column for some reason, adding them manually
row.names(nasoniadevlist_splicing$E)<-nasoniadevlist_splicing$genes
colnames(nasoniadevlist_splicing$E)<-nasoniadevlist_splicing$targets
# Converting to LIMMA EList format
nasdevEList_splicing<-new("EList", nasoniadevlist_splicing)
# remove list (disabled for diagnostics)
rm(list=c("nasoniadevlist_splicing"))

### MALE FEMALES CONTRASTS WRONG!, must re-code
# # Create complete design matrix
# # Create factors for matrix
# stage<-factor(as.character(regmatches(nasdevEList_splicing$targets, gregexpr("^[[:alnum:]]+",nasdevEList_splicing$targets))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
# stage_ordered<-ordered(stage, levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
# male<-factor(grepl("[^[:alnum:]]male|testes", nasdevEList_splicing$targets))
# # Draw matrix
# design<-model.matrix(~0+stage/male, contrasts.arg=list(stage="contr.sum",male="contr.sum"))
# 
# ## fit models
# fit_splicing<-lmFit(nasdevEList_splicing, design)
# fit2_splicing<-eBayes(fit_splicing)

# Save splicing models and EList
save(file="./Output/Splicing_Basestats_fit2_LIMMAonly.RData", list=c("fit2_splicing"))

###### we now apply FDR correction to each factor
library(fdrtool)
## This version applies FDR correction to every factor separately
#FDRfits<-apply(fit2$p.value, 2, function(x){fdrtool(x, statistic="pvalue")$qval})

# This version applies FDR correction to all p-values
FDRfitall_splicing<-fdrtool(as.vector(fit2_splicing$p.value), statistic="pvalue", verbose=T)
Fdrfitall_splicing<-matrix(FDRfitall_splicing$qval, ncol=ncol(fit2_splicing$p.value), nrow=nrow(fit2_splicing$p.value), dimnames=list(fit2_splicing$genes, colnames(fit2_splicing$p.value))) # add global (tail area based)
fdrfitall_splicing<-matrix(FDRfitall_splicing$lfdr, ncol=ncol(fit2_splicing$p.value), nrow=nrow(fit2_splicing$p.value), dimnames=list(fit2_splicing$genes, colnames(fit2_splicing$p.value))) # local fdr (density based)

# Save Fdr and fdr to base model
fit2_splicing$Fdr<-Fdrfitall_splicing
fit2_splicing$fdr<-fdrfitall_splicing

# and save as csv
write.csv(Fdrfitall_splicing, file="./Output/NasoniaDev_Splicing_limma_qval_unsigned.csv")

# Add sign of sex bias
fitq_splicing<-Fdrfitall_splicing
fitq_splicing[,6:10]<-fitq_splicing[,6:10]*sign(fit2_splicing$coefficients[,6:10])
# and save as csv
write.csv(fitq_splicing, file="./Output/NasoniaDev_Splicing_limma_qval_signed.csv")

# and save updated model
save(file="./Output/Splicing_BaseStats_fit2_fdrtool.RData", list=c("fit2_splicing"))

### Save list of GENES with alternative splicing in any treatment

# # Save lower Fdr per exon per condition for each gene to original models
library(plyr)
Fdr_splicing<-as.data.frame(fit2_splicing$Fdr)
Fdr_splicing$geneID<-unlist(regmatches(row.names(Fdr_splicing), gregexpr("^[^t]*",row.names(Fdr_splicing))))
Fdr_splicing<-ddply(Fdr_splicing, .variable=.(geneID), numcolwise(min))
row.names(Fdr_splicing)<-Fdr_splicing$geneID
Fdr_splicing<-as.matrix(subset(Fdr_splicing, select=-geneID))
# add to limma results
fit2_splicing$Fdr_splicing<-Fdr_splicing


# Select lowest EXON fdr for each GENE (significant alternative splicing in that treatment)
fdr_splicing<-as.data.frame(fit2_splicing$fdr)
fdr_splicing$geneID<-as.factor(unlist(regmatches(row.names(fdr_splicing), gregexpr("^[^t]*",row.names(fdr_splicing)))))
fdr_splicing<-ddply(fdr_splicing, .variable=.(geneID), numcolwise(min))
row.names(fdr_splicing)<-fdr_splicing$geneID
fdr_splicing<-as.matrix(subset(fdr_splicing, select=-geneID))
# add to limma results
fit2_splicing$fdr_splicing<-fdr_splicing
# and save as csv
write.csv(fdr_splicing, file="./Output/NasoniaDev_Gene_splicing_fdr.csv")


# Store position of minimum q-value exon per treatment
splicing_T<-as.data.frame(fit2_splicing$fdr)
splicing_T$geneID<-as.factor(unlist(regmatches(row.names(splicing_T), gregexpr("^[^t]*",row.names(splicing_T)))))
splicing_T<-ddply(splicing_T, .variable=.(geneID), numcolwise(which.min))
# Create dataframe of coefficients
splicing_coefs<-as.data.frame(fit2_splicing$coefficients)
splicing_coefs$geneID<-unlist(regmatches(row.names(splicing_coefs), gregexpr("^[^t]*",row.names(splicing_coefs))))
# Create data.frame with all geneIDs on array, set gene counter to zero
# Main loop: For all gene IDs in qvalue dataframe, add gene position in coef to min qval exon position to get row position in coef dataset
# Internal loop: For each minimum qvalue estimate value add one column to dataframe
splicing_coefs_2<-data.frame(row.names=unique(splicing_coefs$geneID))
splicing_coefs_2<-data.frame(matrix(nrow=length(unique(splicing_coefs$geneID)), ncol=10, dimnames=list(row=unique(splicing_coefs$geneID), column=colnames(splicing_coefs)[1:10])))
i=0
for (gID in row.names(splicing_coefs_2)){
  i=i+1 # row to be used for storing coefficients
  j=0 # column counter
  rows<-unlist(which(splicing_coefs$geneID%in%gID)[1]+splicing_T[which(splicing_T$geneID==gID),2:11]-1)
  for (r in rows){
    j<-j+1 # Column counter
    splicing_coefs_2[i, j]<-splicing_coefs[r, j] # add to the gene row the coefficient in the appropriate exon row, column by column
  }
}
colnames(splicing_coefs_2)<-gsub(".",":",colnames(splicing_coefs_2), fixed=T)
# add to limma results
fit2_splicing$coefficients_splicing<-as.matrix(splicing_coefs_2)
# Save workspace
save.image(file="./Output/Splicing_Basestats.RData")
# Save updated gene models
save(file="./Output/Splicing_Basestats_fit2_gene.RData", list=c("fit2_splicing"))

#### Repeat splicing modelling for Gonad specific isoforms

#### Load LIMMA results
load(file="./Output/BaseStats_gonadfit2.RData")
# Load exon data
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON
# filter only rows with samples (exclude means)
nasoniagonadexon<-nasonia_devtesova.exon.r99[,grep("[.][123]",names(nasonia_devtesova.exon.r99))]
# Exclude non-gonads and non-adults
nasoniagonadexon<-nasoniagonadexon[,grep("(adult|ovaries|testes)", names(nasoniagonadexon))]
# Sava dataset before analyses
write.csv(nasoniagonadexon, file="./Output/nasoniagonadexon.csv")

#### Subtract estimates from samples

# Annotate samples with male, stage and replicate
testData<-nasoniagonadexon
testData$exonID<-row.names(testData)
# Gene names are defined as the portion of the exon name that goes from the beginning to the last non "t" character
testData$geneID<-unlist(regmatches(testData$exonID, gregexpr("^[^t]*",testData$exonID)))

library(reshape)
testData<-melt(testData)
testData$adult<-as.numeric(grepl("male",testData$variable))
testData$adult<-as.factor(ifelse(testData$adult==1, "adult", "gonads"))
testData$male<-grepl("_male|testes", testData$variable)
testData$male<-as.factor(ifelse(testData$male==T, "male", "female"))
testData$replicate<-factor(as.character(regmatches(testData$variable, gregexpr("[[:digit:]]$",testData$variable))))
testData<-subset(testData, select=-c(variable))

# Load and annotate coefficients wit whole animal vs testes/ovaries and testes vs ovaries
testCoefs<-gonadfit2$coefficients[row.names(gonadfit2$coefficients)%in%testData$geneID,]
testCoefs<-melt(testCoefs)
names(testCoefs)<-c("geneID","variable","Coef")
testCoefs$Testes<-grepl("Testes",testCoefs$variable)
testCoefs$Ovaries<-grepl("Ovaries",testCoefs$variable)
testCoefs$Gonads<-grepl("Gonads",testCoefs$variable)
testCoefs<-subset(testCoefs, select=-c(variable))
# # Retain only coefficients for testes, ovaries and testes vs ovaries
testCoefs<-testCoefs[which(apply(testCoefs[,3:5], 1, function(x){any(x==T)})),]
# create factor testes, ovaries or gonads
testCoefs$Factor<-as.factor(with(testCoefs, ifelse(Testes==T, "Testes", ifelse(Ovaries==T, "Ovaries", ifelse(Gonads==T, "Gonads", NA)))))

# Subtract testes coefficients to testes, ovaries to ovaries
testData2<-testData[which(testData$adult=="gonads"),]
testData2$Factor<-as.factor(with(testData2, ifelse(adult=="gonads", ifelse(male=="male", "Testes", "Ovaries"), NA)))
testData2<-merge(testData2, testCoefs, by=c("geneID","Factor"),all.x=T)
testData2$value<-testData2$value-testData2$Coef
testData2<-testData2[,c("exonID","geneID","value","Factor","replicate")]

# Subtract gonad (testes-ovaries) specific coefficients to testes
testData3<-testData2[which(testData2$Factor=="Testes"),]
testData3<-merge(testData3, testCoefs[which(testCoefs$Gonads==T),], by=c("geneID"), all.x=T)
testData3$value<-testData3$value-testData3$Coef
testData3<-testData3[,c("exonID","geneID","value","Factor.x","replicate")]
names(testData3)<-c("exonID","geneID","value","Factor","replicate")

# Add ovaries and adult data
testData3<-rbind(testData3, testData2[which(testData2$Factor=="Ovaries"),])
testData3$adult<-"gonads"
testData3$male<-ifelse(grepl("Testes", testData3$Factor), "male","female")
testData3<-testData3[,c("exonID","geneID","value","adult","male","replicate")]
testData3<-rbind(testData3, testData[which(testData$adult=="adult"),])

# Put into exon by sample format
testData3$male<-as.factor(testData3$male)
testData3$adult<-as.factor(testData3$adult)
nasoniagonadsplicing2<-cast(testData3, formula=exonID~adult+male+replicate)
row.names(nasoniagonadsplicing2)<-nasoniagonadsplicing2$exonID
nasoniagonadsplicing2<-subset(nasoniagonadsplicing2, select=-exonID)

# Save as CSV and clean workspace
write.csv(nasoniagonadsplicing2, file="./Output/Nasonia_Gonad_Splicing.csv")
rm(list=(grep("test", ls(), value=T)))

#### Create LIMMA object
library(limma)
# Reformat as list
nasoniagonadlist_splicing<-list(E=as.matrix(nasoniagonadsplicing2), genes=as.character(row.names(nasoniagonadsplicing2)), targets=names(nasoniagonadsplicing2))
## BEAR does not add gene names and sample names to rows and column for some reason, adding them manually
row.names(nasoniagonadlist_splicing$E)<-nasoniagonadlist_splicing$genes
colnames(nasoniagonadlist_splicing$E)<-nasoniagonadlist_splicing$targets
# Converting to LIMMA EList format
nasgonadEList_splicing<-new("EList", nasoniagonadlist_splicing)
# remove list (disabled for diagnostics)
rm(list=c("nasoniagonadlist_splicing"))

# Create complete design matrix
# Create factors for matrix
male<-ifelse(grepl("_male", nasgonadEList_splicing$targets), "male","female")
gonads<-ifelse(grepl("gonads", nasgonadEList_splicing$targets), "gonads","adult")
adultmale<-as.factor(paste(gonads, male, sep="_"))
  
# Draw matrix
design<-model.matrix(~0+adultmale)
colnames(design)<-levels(adultmale)

# fit limma model
gonadfit_splicing<-lmFit(nasgonadEList_splicing, design)

# apply relevant contrasts (males vs females, whole adults (M and F) vs gonads (M and F), F gonads vs M gonads)
cont.matrix<-makeContrasts(
  Males=adult_female-adult_male,
  Testes=adult_male-gonads_male,
  Ovaries=adult_female-gonads_female,
  Gonads=gonads_female-gonads_male,
  levels=design)


# re-fit model to new contrasts
gonadfit2_splicing<-contrasts.fit(gonadfit_splicing, cont.matrix)
gonadfit2_splicing<-eBayes(gonadfit2_splicing)

# Save splicing models and EList
save(file="./Output/Splicing_Basestats_gonadfit2_LIMMAonly.RData", list=c("gonadfit2_splicing"))

###### we now apply FDR correction to each factor
library(fdrtool)
## This version applies FDR correction to every factor separately
#FDRfits<-apply(fit2$p.value, 2, function(x){fdrtool(x, statistic="pvalue")$qval})

# This version applies FDR correction to all p-values
FDRfitallgonads_splicing<-fdrtool(as.vector(gonadfit2_splicing$p.value), statistic="pvalue", verbose=T)
Fdrfitallgonads_splicing<-matrix(FDRfitallgonads_splicing$qval, ncol=ncol(gonadfit2_splicing$p.value), nrow=nrow(gonadfit2_splicing$p.value), dimnames=list(gonadfit2_splicing$genes, colnames(gonadfit2_splicing$p.value))) # add global (tail area based)
fdrfitallgonads_splicing<-matrix(FDRfitallgonads_splicing$lfdr, ncol=ncol(gonadfit2_splicing$p.value), nrow=nrow(gonadfit2_splicing$p.value), dimnames=list(gonadfit2_splicing$genes, colnames(gonadfit2_splicing$p.value))) # add global (tail area based)

# Save Fdr and fdr to base model
gonadfit2_splicing$Fdr<-Fdrfitallgonads_splicing
gonadfit2_splicing$fdr<-fdrfitallgonads_splicing

# and save as csv
write.csv(Fdrfitallgonads_splicing, file="./Output/NasoniaGonad_Splicing_limma_qval_unsigned.csv")
