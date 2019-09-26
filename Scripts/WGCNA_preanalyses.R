##### WGCNA preanalyses
rm(list=ls())
library(plyr)
### Complete version of genes plus net exons (exon minus gene value for the same sample)

## import gene matrix
nasonia_devtesova.gene.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.gene.r99.score")
row.names(nasonia_devtesova.gene.r99)<-nasonia_devtesova.gene.r99$GENE
# filter only rows with samples (exclude means)
nasdevgene<-nasonia_devtesova.gene.r99[,grep("[.][123]",names(nasonia_devtesova.gene.r99))]
# Exclude gonads 
nasdevgene<-nasdevgene[,-grep("(ovaries|testes)", names(nasdevgene))]

## import exon matrix
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON
# filter only rows with samples (exclude means)
nasdevexon<-nasonia_devtesova.exon.r99[,grep("[.][123]",names(nasonia_devtesova.exon.r99))]
# Exclude gonads 
nasdevexon<-nasdevexon[,-grep("(ovaries|testes)", names(nasdevexon))]
# add gene IDs
nasdevexon$geneID<-unlist(regmatches(row.names(nasdevexon), gregexpr("^[^t]*",row.names(nasdevexon))))
# calculate gene medians
nasdevintercepts<-ddply(nasdevexon, .(geneID), numcolwise(median))
# subtract intercepts from exons
nasdevintercepts2<-merge(data.frame(geneID=nasdevexon$geneID), nasdevintercepts, by="geneID", all.x=T)
nasdevintercepts2<-nasdevintercepts2[,-grep("geneID", names(nasdevintercepts2))]
nasdevexon2<-nasdevexon[,-grep("geneID", names(nasdevexon))]
## Check if column names match between two sets (should return all 0)
match(names(nasdevexon2), names(nasdevintercepts2))-(1:30)
## subtract
nasdevsplicing<-nasdevexon2-nasdevintercepts2

## merge gene and net exon datasets
nasdevgeneexon<-rbind(nasdevgene, nasdevsplicing)

##WGCNA based quality control
library(WGCNA)
allowWGCNAThreads()
# check genes with missing values or zero variance, then remove them
QC<-goodSamplesGenes(t(nasdevgeneexon), verbose=3)
t_nasdevgeneexon<-t(nasdevgeneexon)[,which(QC$goodGenes)]
# # cluster samples to detect anomalies (computationally intensive)
# h1<-flashClust(dist(t_nasdevgeneexon, method="euclidean"), method="average")
# pdf(file="./Graphics/WGCNA_Diagnostic_Sample_Clustering.pdf", paper = "a4")
# plot(as.phylo(h1), tip.color=grepl("female", h1$labels)+1, type="cladogram")
# dev.off()

### Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2), seq(from=25, to=50, by=5))
# Call the network topology analysis function
sft = pickSoftThreshold(t_nasdevgeneexon, powerVector = powers, verbose = 5, moreNetworkConcepts = T)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file="./Graphics/WGCNA_power_finder.pdf", paper = "a4")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

## Save datasets
save(nasdevgeneexon, file="./Output/nasdevgenexon_pre_WGCNA.RData")
save(t_nasdevgeneexon, file="./Output/t_nasdevgenexon_pre_WGCNA.RData")
write.csv(nasdevgeneexon, file="./Output/nasdevgenexon_pre_WGCNA.csv")
write.csv(t_nasdevgeneexon, file="./Output/t_nasdevgenexon_pre_WGCNA.csv")
