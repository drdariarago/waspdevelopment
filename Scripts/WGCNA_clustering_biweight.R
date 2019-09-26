## Create modules blockwise
# load core packages
date()
rm(list=ls())
library(ggplot2)
library(reshape2)
library(flashClust)
library(WGCNA)
library(plyr)
allowWGCNAThreads()
source(file = "./Scripts/multiplot.R")

logit<-function(x){
  log10(x/(1-x))
}

z<-function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
}

# initialize output path
newdir<-file.path(getwd(), "Output/WGCNA_clustering_biweight")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/WGCNA_clustering_biweight")
dir.create(graphdir)

# import dataset
nasdevgeneexon<-read.csv(file="./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexon_evalues.csv")
row.names(nasdevgeneexon)<-nasdevgeneexon$eigenexonID
# nasdevgeneexon<-nasdevgeneexon[1:1000,] # reduce size for testing
nasdevgeneexon<-nasdevgeneexon[,which(colnames(nasdevgeneexon)%in%setdiff(names(nasdevgeneexon),"eigenexonID"))]
t_nasdevgeneexon<-transposeBigData(nasdevgeneexon, blocksize = 50000)
# QC 
# Check for zero variance nodes and nodes with too many missing values
goodCols<-goodSamplesGenes(datExpr = t_nasdevgeneexon, minNSamples = 3, minFraction = 1/10)
lapply(goodCols, function(x) {table(x)})
lapply(goodCols, function(x) {prop.table(table(x))})
t_nasdevgeneexon[,which(goodCols$goodGenes==F)]
t_nasdevgeneexon<-t_nasdevgeneexon[,which(goodCols$goodGenes==T)]
# How many gene and how many transcript values were removed?
table(grepl("con", row.names(nasdevgeneexon)))-table(grepl("con", colnames(t_nasdevgeneexon)))
# Normalize transcription and splicing nodes separately
m_nasdevgeneexon<-data.frame(melt(t_nasdevgeneexon),
                             sampleID=row.names(t_nasdevgeneexon))
m_nasdevgeneexon$con<-grepl("_con", m_nasdevgeneexon$variable)

distr_pre<-ggplot(data=m_nasdevgeneexon, aes(x=value, col=con))+geom_density()+theme_bw()

ddply(.data = m_nasdevgeneexon, .variables = .(con), summarize, 
      mean=mean(value, na.rm = T),
      sd=sd(value, na.rm = T),
      zero=(0-mean(value, na.rm = T))/sd(value, na.rm = T))

m_nasdevgeneexon_z<-ddply(.data = m_nasdevgeneexon, .variables = .(con), summarize, 
                        variable=variable, 
                        value=z(value),
                        sampleID=sampleID)

distr_post<-ggplot(data=m_nasdevgeneexon_z, aes(x=value, col=con))+geom_density()+theme_bw()

t_nasdevgeneexon<-recast(m_nasdevgeneexon_z[,2:4], sampleID~variable)
row.names(t_nasdevgeneexon)<-t_nasdevgeneexon$sampleID
t_nasdevgeneexon<-t_nasdevgeneexon[,-1]

# Print distributions
pdf(file = file.path(graphdir, "distr_scaling.pdf"), paper = "a4r", useDingbats = F)
multiplot(distr_pre, distr_post, cols = 1)
dev.off()

# Save transposed and normalized expression dataset
save(t_nasdevgeneexon, file = file.path(newdir, "t_nasdevgeneexon.RData"))
write.csv(nasdevgeneexon, file = file.path(newdir, "nasdevgeneexon.csv"))

# calculate block size (check RAM availability trough unix command vmstat, RAM expressed in bytes)
bs<-blockSize(matrixSize = ncol(t_nasdevgeneexon), maxMemoryAllocation = 509435*1048576, overheadFactor = 4, rectangularBlocks = T)
bs
# Block size is hard limited to a maximum of 46340
bs<-ifelse(bs>46339, 46339, bs)
bs

# pick optimal power threshold with biweight midcorrelations
threshold_biweight<-pickSoftThreshold(
  corFnc = bicor,
#   blockSize = bs,
  t_nasdevgeneexon, 
  corOptions = list(quick=0, use = "pairwise.complete.obs"),
  powerVector = c(seq(from =1, to = 21, by = 5), 22:23, 30),
  networkType="signed",
  moreNetworkConcepts=T,
  verbose=0)
# save threshold selection
write.csv(threshold_biweight, file=file.path(newdir, "WGCNA_biweight_power_selection.csv"))

# Power selection plot
power_1<-ggplot(threshold_biweight$fitIndices, aes(x=Power, y=-sign(slope)*SFT.R.sq))+geom_text(aes(label=Power))+geom_hline(yintercept=0.9)+theme_bw()
power_2<-ggplot(melt(threshold_biweight$fitIndices[,c(1,grep("k", names(threshold_biweight$fitIndices)))], id.vars = "Power"), aes(x=Power, col=variable, y=value, label=Power))+geom_text()+geom_line()+scale_y_log10()+theme_bw()

pdf(file=file.path(graphdir, "WGCNA_biweight_Power_diagnostics.pdf"), paper = "a4r")
multiplot(power_1, power_2, cols = 2)
dev.off()

# Calculate block assignments
Kblocks<-projectiveKMeans(datExpr = t_nasdevgeneexon, preferredSize = bs, sizePenaltyPower = 2, networkType = "signed", checkData = T, verbose = 0)

# create list of sub-datasets for simplified calculations
Eval_blocks<-lapply(1:length(unique(Kblocks$clusters)), function(x){t_nasdevgeneexon[,which(Kblocks$clusters==x)]})
# validation of block assignments: subtracts indexes from blocks to indexes of matrices on main dataset (should return only zeroes)
# sapply(1:length(unique(Kblocks$clusters)), function(x){which(colnames(t_nasdevgeneexon)%in%colnames(Eval_blocks[[x]]))-which(Kblocks$clusters==x)})
all(unlist(sapply(1:length(unique(Kblocks$clusters)), function(x){which(colnames(t_nasdevgeneexon)%in%colnames(Eval_blocks[[x]]))-which(Kblocks$clusters==x)}))==0)

#Calculate adjacency matrices with chosen power level
# alternatively, use GTOMdist() (generalized topological overlap matrix, accounts for higher order interactions)
adj<-lapply(Eval_blocks, function(x){
  adjacency(x, 
            power = threshold_biweight$powerEstimate, 
            corFnc = "bicor",
            corOptions = "quick=0, use = 'pairwise.complete.obs'",
            type = "signed")
})
str(adj)
save(adj, file = file.path(newdir, "/adj"), ascii = F)

#Convert to topological dissimilarity matrix
dTOM<-lapply(adj, FUN = TOMdist, TOMType="signed")
#Add names from original blocks
for (block in 1:length(dTOM)){
  colnames(dTOM[[block]])<-colnames(adj[[block]])
  row.names(dTOM[[block]])<-row.names(adj[[block]])
}
str(dTOM)
# Save TOM matrix
save(dTOM, file = file.path(newdir, "/dTOM"), ascii = F)

# Record gene names for each cluster
transcriptIDs<-lapply(dTOM, row.names)
#Calculate tree for clustering
geneTree<-lapply(X = dTOM, FUN = function(x){flashClust(as.dist(x), method="average")})
#Create clusters
dynamicMods<-mapply(function(den, dis){cutreeDynamic(dendro = den, distM = dis, deepSplit = 4, pamRespectsDendro = F, minClusterSize = 20)}, 
                    den=geneTree, dis=dTOM)
str(dynamicMods)
# Generate unique labels for each cluster (final digits represent block number)
clusterassignments<-mapply(function(x,y){paste0(y,x)}, x=1:length(dynamicMods), y=dynamicMods)
# Fuse bloks
clusterassignments<-unlist(clusterassignments)
# Mark unassigned clusters with zeroes
clusterassignments<-as.numeric(ifelse(grepl("^0", clusterassignments)==T,0,clusterassignments))
# convert to numeric labels to colours
colclusterassignments<-labels2colors(clusterassignments, zeroIsGrey = T)
# fuse cluster IDs 
transcriptIDs<-unlist(transcriptIDs)
# reorder expression data based on cluster IDs
t_nasdevgeneexon<-t_nasdevgeneexon[,match(transcriptIDs,colnames(t_nasdevgeneexon))]

# Calculate module eigengenes
MEList<-moduleEigengenes(expr = t_nasdevgeneexon, colors = colclusterassignments, grey = "grey", excludeGrey = F, nPC = 10)
# Cluster by eigengene dissimilarity
MEtree<-flashClust(as.dist(1-bicor(MEList$eigengenes, quick=0, use = "pairwise.complete.obs")), method = "average")
# Plot correlations
pdf(file = file.path(graphdir, "WGCNA_Cross_Cluster_Correlation_beforemerging.pdf"), width = 46.8, height = 33.1)
plot(MEtree, main="Cross-Cluster Correlations", xlab="", ylab="Proportion of Correlation Divergence")
abline(a = 0.2, b=0, col="green")
dev.off()

# merge highly correlated modules
cut=0.1
MElist2<-mergeCloseModules(exprData = t_nasdevgeneexon, colors = colclusterassignments, 
                           corFnc = bicor, corOptions = list(quick=0, use="pairwise.complete.obs"), 
                           cutHeight = cut, relabel = T,
                           unassdColor = "grey", getNewMEs = T, getNewUnassdME = T)
colclusterassignments2<-MElist2$colors

# Save cluster information
save(MElist2, file = file.path(newdir, "MEdata_aftermerging"), ascii = F, compress = T)

# Save workspace
rm(adj)
rm(dTOM)
save.image(file=file.path(newdir, "WGCNA_clustering_biweight.RData"), ascii = F, compress = T, safe = T)

# create and save data.frame with gene/cluster assignments
clusterassignments<-data.frame(transcriptID=colnames(t_nasdevgeneexon), clusterID=colclusterassignments2)
write.csv(clusterassignments, file=file.path(newdir, "clusterassignments_simple.csv"))

# create and save data.frame with gene/cluster assigments and expression data
clusterassignments2<-data.frame(transposeBigData(t_nasdevgeneexon), clusterID=colclusterassignments2)
write.csv(clusterassignments2, file=file.path(newdir, "clusterassignments_full.csv"))

# plot correlation heatmap
pdf(file = file.path(graphdir, file="Cluster_heatmaps.pdf"), width = 66, height = 66)
ggplot(data=melt(bicor(MElist2$oldMEs)), aes(x=Var1, y=Var2, fill=value))+geom_tile()+scale_fill_gradient2()
ggplot(data=melt(bicor(MElist2$newMEs)), aes(x=Var1, y=Var2, fill=value))+geom_tile()+scale_fill_gradient2()
dev.off()

# plot clustering trees for each module
pdf(file=file.path(graphdir, "WGCNA_biweight_Cluster_trees.pdf"), width = 46.8, height = 33.1)
lapply(geneTree, plot)
dev.off()

# plot the clusters before and after merging
pdf(file = file.path(graphdir, "WGCNA_Cross_Cluster_Correlation_aftermerging.pdf"), width = 46.8, height = 33.1)
plot(MElist2$oldDendro)
abline(cut,0, col="red")
plot(MElist2$dendro)
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
# 
# # plot the clusters before and after merging
# pdf(file = file.path(graphdir, "WGCNA_Cross_Cluster_Correlation_aftermerging.pdf"), paper = "a4")
# plotDendroAndColors(flashClust(as.dist(1-TOMsimilarityFromExpr(t_nasdevgeneexon, power = threshold_biweight$powerEstimate, TOMType = "signed", corType = "bicor"))), 
#                     cbind(colclusterassignments, colclusterassignments2),
#                     dendroLabels = F)
# dev.off()
