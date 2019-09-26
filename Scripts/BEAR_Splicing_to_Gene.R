### OUTDATED: now part of Splicing Basestats

## BEAR version of conversion from splicing dataset to gene-based values (will be incorporated into splicing_BaseStats)

# load model in BEAR
load(file="./Output/Splicing_BaseStats_fit2_fdrtool.RData")

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
splicing_coefs$geneID<-<-unlist(regmatches(row.names(splicing_coefs), gregexpr("^[^t]*",row.names(splicing_coefs))))
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
save.image(file="./Output/Splcing_Basestats.RData")
# Save updated gene models
save(file="./Output/Splicing_Basestats_fit2_gene.RData", list=c("fit2_splicing"))