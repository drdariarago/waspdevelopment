## Create modules blockwise
# load core packages
date()
rm(list=ls())
library(flashClust)
library(WGCNA)
library(ggplot2)
library(reshape2)
allowWGCNAThreads()
# initialize output path
newdir<-file.path(getwd(), "Output/WGCNA_clustering_biweight_singleblock")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/WGCNA_clustering_biweight_singleblock")
dir.create(graphdir)

# import dataset
nasdevgeneexon<-read.csv(file="./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexon_evalues.csv")
row.names(nasdevgeneexon)<-nasdevgeneexon$eigenexonID
# nasdevgeneexon<-nasdevgeneexon[1:5000,] # reduce size for testing
nasdevgeneexon<-nasdevgeneexon[,which(colnames(nasdevgeneexon)%in%setdiff(names(nasdevgeneexon),"eigenexonID"))]
t_nasdevgeneexon<-transposeBigData(nasdevgeneexon, blocksize = 50000)

# QC 
# Remove nodes expressed in less than 3 samples
goodCols<-apply(t_nasdevgeneexon, 2, function(x){sum(na.exclude(x)>0)>2})
summary(goodCols)
t_nasdevgeneexon<-t_nasdevgeneexon[,which(goodCols==T)]
# Check for zero variance nodes
goodCols<-goodSamplesGenes(datExpr = t_nasdevgeneexon, minNSamples = 3, minFraction = 1/10)
lapply(goodCols, function(x) {table(x)})
lapply(goodCols, function(x) {prop.table(table(x))})
t_nasdevgeneexon<-t_nasdevgeneexon[,which(goodCols$goodGenes==T)]
str(t_nasdevgeneexon)
# How many gene and how many transcript values were removed?
table(grepl("con", row.names(nasdevgeneexon)))-table(grepl("con", colnames(t_nasdevgeneexon)))


# calculate block size (check RAM availability trough unix command vmstat, RAM expressed in bytes)
bs<-blockSize(matrixSize = ncol(t_nasdevgeneexon), maxMemoryAllocation = 509435*1048576, overheadFactor = 4, rectangularBlocks = T)
bs
# Block size is hard limited to a maximum of 46340
bs<-ifelse(bs>46340, 46340, bs)
bs

# pick optimal power threshold with biweight midcorrelations
threshold_biweight<-pickSoftThreshold(
  corFnc = bicor,
#   blockSize = bs,
  t_nasdevgeneexon, 
  corOptions = list(quick=0, use = "pairwise.complete.obs"),
  powerVector = c(seq(from = 1,to = 21,by = 10),seq(from = 22, to = 31, by = 3),32:35),
  networkType="signed",
  moreNetworkConcepts=T,
  verbose=1)
# save threshold selection
write.csv(threshold_biweight, file=file.path(newdir, "WGCNA_biweight_power_selection.csv"))

# Power selection plot
pdf(file=file.path(graphdir, "WGCNA_biweight_Power_diagnostics.pdf"))
ggplot(threshold_biweight$fitIndices, aes(x=Power, y=-sign(slope)*SFT.R.sq))+geom_text(aes(label=Power))+geom_hline(yintercept=0.9)+theme_bw()
ggplot(melt(threshold_biweight$fitIndices[,c(1,grep("k", names(threshold_biweight$fitIndices)))], id.vars = "Power"), aes(x=Power, col=variable, y=value, label=Power))+geom_text()+geom_line()+scale_y_log10()+theme_bw()
dev.off()

# Final selected power
threshold_biweight$powerEstimate

#Calculate adjacency matrix with chosen power level
# alternatively, use GTOMdist() (generalized topological overlap matrix, accounts for higher order interactions)
adj<-adjacency(datExpr = t_nasdevgeneexon, 
               power = threshold_biweight$powerEstimate, 
               corFnc = "bicor",
               corOptions = "quick=0, use = 'pairwise.complete.obs'",
               type = "signed")
str(adj)
save(adj, file = file.path(newdir, "/adj"), ascii = F)

#Convert to topological dissimilarity matrix
dTOM<-TOMdist(adjMat = adj, TOMType = "signed")
# Transfer transcript names from adj to dTOM
row.names(dTOM)<-row.names(adj)
colnames(dTOM)<-colnames(adj)
str(dTOM)

# Save TOM matrix
save(dTOM, file = file.path(newdir, "/dTOM"), ascii = F)

# Record gene names
transcriptIDs<-row.names(dTOM)
#Calculate tree for clustering
geneTree<-flashClust(d = as.dist(dTOM), method = "average")
#Create clusters
dynamicMods<-cutreeDynamic(dendro = geneTree, distM = dTOM, deepSplit = 4, pamRespectsDendro = F, minClusterSize = 20)
# Convert numeric labels to colours
colclusterassignments<-labels2colors(dynamicMods)

# reorder expression data based on cluster IDs
any(match(transcriptIDs,colnames(t_nasdevgeneexon))-1:length(t_nasdevgeneexon)!=0)
t_nasdevgeneexon<-t_nasdevgeneexon[,match(transcriptIDs,colnames(t_nasdevgeneexon))]

# Calculate module eigengenes
MEList<-moduleEigengenes(expr = t_nasdevgeneexon, colors = colclusterassignments)

# Cluster by eigengene dissimilarity
MEdist<-as.dist(1-bicor(MEList$eigengenes, quick=0, use = "pairwise.complete.obs"))
MEtree<-flashClust(MEdist, method = "average")

# Plot correlations
pdf(file = file.path(graphdir, "WGCNA_Cross_Cluster_Correlation_beforemerging.pdf"), paper = "a4")
plot(MEtree, main="Cross-Cluster Correlations", xlab="", ylab="Proportion of Correlation Divergence")
abline(a = 0.1, b=0, col="green")
dev.off()

# merge highly correlated modules
MElist2<-mergeCloseModules(exprData = t_nasdevgeneexon, useAbs = F,
                           colors = colclusterassignments, unassdColor = "grey",
                           corFnc = bicor, corOptions = list(quick=0, use="pairwise.complete.obs"), 
                           cutHeight = 0.1, relabel = T)
colclusterassignments2<-MElist2$colors

#Save workspace
save.image(file=file.path(newdir, "WGCNA_clustering_biweight.RData"), ascii = F, compress = T, safe = T)

# create and save data.frame with gene/cluster assignments
clusterassignments<-data.frame(transcriptID=colnames(t_nasdevgeneexon), clusterID=colclusterassignments2)
write.csv(clusterassignments, file=file.path(newdir, "clusterassignments_simple.csv"))

# create and save data.frame with gene/cluster assigments and expression data
clusterassignments2<-data.frame(transposeBigData(t_nasdevgeneexon), clusterID=colclusterassignments2)
write.csv(clusterassignments2, file=file.path(newdir, "clusterassignments_full.csv"))

# plot clustering trees for each module
pdf(file=file.path(graphdir, "WGCNA_biweight_Cluster_trees.pdf"), paper = "a4r")
plot(geneTree)
dev.off()

# plot the clusters before and after merging
pdf(file = file.path(graphdir, "WGCNA_Cross_Cluster_Correlation_aftermerging.pdf"), paper = "a4r")
plotDendroAndColors(flashClust(as.dist(1-TOMsimilarityFromExpr(t_nasdevgeneexon, power = threshold_biweight$powerEstimate, TOMType = "signed", corType = "bicor"))), 
                    cbind(colclusterassignments, colclusterassignments2),
                    dendroLabels = F)
dev.off()

#### outdated code

# ## re stack TOM matrices
# mainTOMmatrix<-matrix(nrow=ncol(t_nasdevgeneexon), ncol=ncol(t_nasdevgeneexon))
# row.names(mainTOMmatrix)<-unlist(sapply(adj, row.names))
# colnames(mainTOMmatrix)<-unlist(sapply(adj, colnames))
# 
# # match(row.names(dTOM[[1]]), row.names(mainTOMmatrix))
# 
# # calculate starting position of each cluster on the main matrix
# mainMatPos<-cumsum(sapply(X = dTOM, function(x){cumsum(nrow(x))}))
# mainMatPos<-c(0,mainMatPos)
# # insert values into new dataset
# for (i in 1:length(dTOM)){
#   mainTOMmatrix[(mainMatPos[i]+1):mainMatPos[i+1],(mainMatPos[i]+1):mainMatPos[i+1]]<-dTOM[[i]]
# }