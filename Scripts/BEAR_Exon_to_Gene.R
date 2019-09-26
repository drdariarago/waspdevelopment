#### OUTDATED now part of Exon_Basestats

## Update gene models with information on which genes contain DE exons
## BEAR version of conversion from exon dataset to gene-based values (will be incorporated into Exon_BaseStats)

# load model in BEAR
load(file="./Output/Exon_BaseStats_fit2_fdrtool.RData")

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
exon_coefs$geneID<-<-unlist(regmatches(row.names(exon_coefs), gregexpr("^[^t]*",row.names(exon_coefs))))
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
colnames(exon_coefs_2)<-gsub(".",":",colnames(exon_coefs_2), fixed=T)
# add to limma results
fit2_exon$coefficients_exon<-as.matrix(exon_coefs_2)

# Save workspace
save.image(file="./Output/Splcing_Basestats.RData")
# Save updated gene models
save(file="./Output/Exon_Basestats_fit2_gene.RData", list=c("fit2_exon"))