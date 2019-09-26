## Method for detection of splicing without assumptions about sample identity/modelling:
## Workflow: start from exon intensity files, calculate median and SD of each gene in each sample, subtract means to raw exons
## Caveats: in the null case (no splicing), sd goes down with each exon in gene model, single exon genes have zero degrees of freedom, does not consider information from the same gene in other samples
## Pros: No assumption about sample identity required, mean is estimated directly from data (no problems with different normalization)

## Import data
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON

# filter only rows with samples (exclude means)
nasdevexon<-nasonia_devtesova.exon.r99[,grep("[.][123]",names(nasonia_devtesova.exon.r99))]
# Exclude gonads 
nasdevexon<-nasdevexon[,-grep("(ovaries|testes)", names(nasdevexon))]
# add gene IDs
nasdevexon$geneID<-unlist(regmatches(row.names(nasdevexon), gregexpr("^[^t]*",row.names(nasdevexon))))

## Calculate medians and MAD (Median Absolute Deviation (from total median), from http://www.r-bloggers.com/absolute-deviation-around-the-median/)
library(plyr)

nasdevintercepts<-ddply(nasdevexon, .(geneID), numcolwise(median))

# Calculate mads and number of exons
nasdevmads<-ddply(nasdevexon, .(geneID), numcolwise(mad))
nasdevNexons<-ddply(nasdevexon, .(geneID), nrow)
# merge and plot
library(reshape)
nasdevmads2<-merge(melt(nasdevmads), melt(nasdevNexons), all.x=T, by="geneID")
names(nasdevmads2)<-c("geneID", "sample", "mad", "V1", "exons")
nasdevmads2<-nasdevmads2[,-4]
library(ggplot2)
tiff(filename="./Graphics/Splicing_intercept_MADvsExons.tiff", width=7.43, height=8.72, units="in", res=500, compression="lzw")
ggplot(nasdevmads2, aes(y=mad, x=log10(exons), col=sample))+geom_smooth()+geom_point()
dev.off()
## MAD increases until genes reach 10 exons, genes with <10 exons will be more conservative in splicing detection
## Store list of genes with zero variance
zeromadgenes<-nasdevmads2[which(nasdevmads2$mad==0),"geneID"]

#### Subtract intercepts from exons
nasdevintercepts2<-merge(data.frame(geneID=nasdevexon$geneID), nasdevintercepts, by="geneID", all.x=T)
nasdevintercepts2<-nasdevintercepts2[,-grep("geneID", names(nasdevintercepts2))]
nasdevexon2<-nasdevexon[,-grep("geneID", names(nasdevexon))]
## Check if column names match between two sets (should return all 0)
match(names(nasdevexon2), names(nasdevintercepts2))-(1:30)
## subtract
nasdevsplicing<-nasdevexon2-nasdevintercepts2
# ## Remove zero variance genes (list too long for grep, using loop)
# zeromadexons<-c()
# for (geneID in unique(zeromadgenes)){
#   a<-grep(geneID, row.names(nasdevsplicing))
#   zeromadexons<-c(a,zeromadexons)
# }

# Save dateset as csv file
write.csv(nasdevsplicing, file="./Output/Nasonia_Dev_Splicing.csv")

## LIMMA analyses
## convert into EList object (list first, then create as EList)
library(limma)
nasoniadevlist_splicing<-list(E=as.matrix(nasdevsplicing), genes=as.character(row.names(nasdevsplicing)), targets=names(nasdevsplicing))
nasdevEList_splicing<-new("EList", nasoniadevlist_splicing)
rm(list=c("nasoniadevlist_splicing"))

# Create complete design matrix
stage<-factor(as.character(regmatches(nasdevEList_splicing$targets, gregexpr("^[[:alnum:]]+",nasdevEList_splicing$targets))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
stage_ordered<-ordered(stage, levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
# male<-grepl("[^[:alnum:]]male|testes", nasdevEList$targets)
# male<-factor(ifelse(male==T, "Male", "Female"),levels=c("Male", "Female"))
male<-factor(ifelse(grepl("[^[:alnum:]]male|testes", nasdevEList_splicing$targets)==T, "Male", "Female"),levels=c("Female", "Male"))

# options(contrasts=c("contr.treatment", "contr.poly"))
library(car)
design<-model.matrix(~0+stage/male, contrasts.arg=list(stage="contr.Sum",male="contr.treatment"))
# check if data matrix is correct
as.data.frame(cbind(nasdevEList_splicing$targets, as.character(stage), as.character(male), design))

# Fit model
fit_splicing<-lmFit(nasdevEList_splicing, design)
fit2_splicing<-eBayes(fit_splicing)

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
# and R object
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
colnames(splicing_coefs_2)<-colnames(fit2_splicing$fdr_splicing) # relies on same order of columns in matrix fdr_splicing, must check
#gsub(".",":",colnames(splicing_coefs_2), fixed=T)
# add to limma results
fit2_splicing$coefficients_splicing<-as.matrix(splicing_coefs_2)
# Save workspace
save.image(file="./Output/Splicing_Basestats.RData")
# Save updated gene models
save(file="./Output/Splicing_Basestats_fit2_fdrtool.RData", list=c("fit2_splicing"))
