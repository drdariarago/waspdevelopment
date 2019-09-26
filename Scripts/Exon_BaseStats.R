# Limma based analyses for Exon differential expression
## load raw data, analyze with LIMMA and add Fdr and fdr for developmental time series and gonads vs adults

# Load exon data
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON
# filter only rows with samples (exclude means)
nasoniadevexon<-nasonia_devtesova.exon.r99[,grep("[.][123]",names(nasonia_devtesova.exon.r99))]
# Exclude gonads 
nasoniadevexon<-nasoniadevexon[,-grep("(ovaries|testes)", names(nasoniadevexon))]

# Save as CSV and clean workspace
write.csv(nasoniadevexon, file="./Output/Nasonia_Dev_Exon.csv")
rm(nasonia_devtesova.exon.r99)

## Create LIMMA object
library(limma)
# Reformat as list
nasoniadevlist_exon<-list(E=as.matrix(nasoniadevexon), genes=as.character(row.names(nasoniadevexon)), targets=names(nasoniadevexon))
## BEAR does not add gene names and sample names to rows and column for some reason, adding them manually
row.names(nasoniadevlist_exon$E)<-nasoniadevlist_exon$genes
colnames(nasoniadevlist_exon$E)<-nasoniadevlist_exon$targets
# Converting to LIMMA EList format
nasdevEList_exon<-new("EList", nasoniadevlist_exon)
# remove list
rm(list=c("nasoniadevlist_exon"))

# Create complete design matrix
stage<-factor(as.character(regmatches(nasdevEList_exon$targets, gregexpr("^[[:alnum:]]+",nasdevEList_exon$targets))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
stage_ordered<-ordered(stage, levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
# male<-grepl("[^[:alnum:]]male|testes", nasdevEList$targets)
# male<-factor(ifelse(male==T, "Male", "Female"),levels=c("Male", "Female"))
male<-factor(ifelse(grepl("[^[:alnum:]]male|testes", nasdevEList_exon$targets)==T, "Male", "Female"),levels=c("Female", "Male"))

# options(contrasts=c("contr.treatment", "contr.poly"))
library(car)
design<-model.matrix(~0+stage/male, contrasts.arg=list(stage="contr.Sum",male="contr.treatment"))
# check if data matrix is correct
as.data.frame(cbind(nasdevEList_exon$targets, as.character(stage), as.character(male), design))

# fit models
fit_exon<-lmFit(nasdevEList_exon, design)
fit2_exon<-eBayes(fit_exon)

# Save exon models and EList
save(file="./Output/Exon_Basestats_fit2_LIMMAonly.RData", list=c("fit2_exon"))

# we now apply FDR correction to each factor
library(fdrtool)
## This version applies FDR correction to every factor separately
#FDRfits<-apply(fit2$p.value, 2, function(x){fdrtool(x, statistic="pvalue")$qval})

# This version applies FDR correction to all p-values
FDRfitall_exon<-fdrtool(as.vector(fit2_exon$p.value), statistic="pvalue", verbose=T)
Fdrfitall_exon<-matrix(FDRfitall_exon$qval, ncol=ncol(fit2_exon$p.value), nrow=nrow(fit2_exon$p.value), dimnames=list(fit2_exon$genes, colnames(fit2_exon$p.value))) # add global (tail area based)
fdrfitall_exon<-matrix(FDRfitall_exon$lfdr, ncol=ncol(fit2_exon$p.value), nrow=nrow(fit2_exon$p.value), dimnames=list(fit2_exon$genes, colnames(fit2_exon$p.value))) # local fdr (density based)

# Save Fdr and fdr to base model
fit2_exon$Fdr<-Fdrfitall_exon
fit2_exon$fdr<-fdrfitall_exon

# and save as csv
write.csv(Fdrfitall_exon, file="./Output/NasoniaDev_Exon_limma_qval_unsigned.csv")

# Add sign of sex bias
fitq_exon<-Fdrfitall_exon
fitq_exon[,6:10]<-fitq_exon[,6:10]*sign(fit2_exon$coefficients[,6:10])
# and save as csv
write.csv(fitq_exon, file="./Output/NasoniaDev_Exon_limma_qval_signed.csv")

# and save updated model
save(file="./Output/Exon_BaseStats_fit2_fdrtool.RData", list=c("fit2_exon"))

## Update gene models with information on which genes contain DE exons
### Save list of GENES with DE exons in any treatment

# # Save lower Fdr per exon per condition for each gene to original models
library(plyr)
Fdr_exon<-as.data.frame(fit2_exon$Fdr)
# Gene names are defined as the portion of the exon name that goes from the beginning to the last non "t" character
Fdr_exon$geneID<-unlist(regmatches(row.names(Fdr_exon), gregexpr("^[^t]*",row.names(Fdr_exon))))
Fdr_exon<-ddply(Fdr_exon, .variable=.(geneID), numcolwise(min))
row.names(Fdr_exon)<-Fdr_exon$geneID
Fdr_exon<-as.matrix(subset(Fdr_exon, select=-geneID))
# add to limma results
fit2_exon$Fdr_exon<-Fdr_exon


# Select lowest EXON fdr for each GENE (significant differentially expressed exon in that treatment)
fdr_exon<-as.data.frame(fit2_exon$fdr)
# Gene names are defined as the portion of the exon name that goes from the beginning to the last non "t" character
fdr_exon$geneID<-unlist(regmatches(row.names(fdr_exon), gregexpr("^[^t]*",row.names(fdr_exon))))
fdr_exon<-ddply(fdr_exon, .variable=.(geneID), numcolwise(min))
row.names(fdr_exon)<-fdr_exon$geneID
fdr_exon<-as.matrix(subset(fdr_exon, select=-geneID))
# add to limma results
fit2_exon$fdr_exon<-fdr_exon
# and save as csv
write.csv(fdr_exon, file="./Output/NasoniaDev_Gene_exon_fdr.csv")


# Store position of minimum q-value exon per treatment
exon_T<-as.data.frame(fit2_exon$fdr)
exon_T$geneID<-unlist(regmatches(row.names(exon_T), gregexpr("^[^t]*",row.names(exon_T))))
exon_T<-ddply(exon_T, .variable=.(geneID), numcolwise(which.min))
# Create dataframe of coefficients
exon_coefs<-as.data.frame(fit2_exon$coefficients)
exon_coefs$geneID<-unlist(regmatches(row.names(exon_coefs), gregexpr("^[^t]*",row.names(exon_coefs))))
# Create data.frame with all geneIDs on array, set gene counter to zero
# Main loop: For all gene IDs in qvalue dataframe, add gene position in coef to min qval exon position to get row position in coef dataset
# Internal loop: For each minimum qvalue estimate value add one column to dataframe
exon_coefs_2<-data.frame(row.names=unique(exon_coefs$geneID))
exon_coefs_2<-data.frame(matrix(nrow=length(unique(exon_coefs$geneID)), ncol=10, dimnames=list(row=unique(exon_coefs$geneID), column=colnames(exon_coefs)[1:10])))
i=0
for (gID in row.names(exon_coefs_2)){
  i=i+1 # row to be used for storing coefficients
  j=0 # column counter
  rows<-unlist(which(exon_coefs$geneID%in%gID)[1]+exon_T[which(exon_T$geneID==gID),2:11]-1)
  for (r in rows){
    j<-j+1 # Column counter
    exon_coefs_2[i, j]<-exon_coefs[r, j] # add to the gene row the coefficient in the appropriate exon row, column by column
  }
}
colnames(exon_coefs_2)<-colnames(fit2_exon$fdr_exon)
  # gsub(".",":",colnames(exon_coefs_2), fixed=T)
# add to limma results
fit2_exon$coefficients_exon<-as.matrix(exon_coefs_2)

# Save workspace
save.image(file="./Output/Splcing_Basestats.RData")
# Save updated gene models
save(file="./Output/Exon_Basestats_fit2_gene.RData", list=c("fit2_exon"))