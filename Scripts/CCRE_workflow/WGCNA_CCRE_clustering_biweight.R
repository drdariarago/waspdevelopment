#### WGCNA CCRE Clustering
#### Creates static network using CCRE data
date()
rm(list=ls())
library(ggplot2)
library(reshape2)
library(flashClust)
library(WGCNA)
library(plyr)
allowWGCNAThreads()
library(lattice)
source(file = "./Scripts/multiplot.R")
logit<-function(x){log10(x/(1-x))}
# initialize output path
newdir<-file.path(getwd(), "Output/WGCNA_CCRE_clustering_biweight")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/WGCNA_CCRE_clustering_biweight")
dir.create(graphdir)

# import dataset
CCRE_nasdevgeneexon <- read.csv(file = "./Output/CCREs/CCRE_nasdevgeneexon.csv")
row.names(CCRE_nasdevgeneexon) <- CCRE_nasdevgeneexon$X
CCRE_nasdevgeneexon <- CCRE_nasdevgeneexon[,-1]
# CCRE_nasdevgeneexon <- CCRE_nasdevgeneexon[sample(1:nrow(CCRE_nasdevgeneexon),1000),]
t_CCRE_nasdevgeneexon <- t(CCRE_nasdevgeneexon)

## pick optimal power threshold with biweight midcorrelations
threshold_biweight<-pickSoftThreshold(
  RsquaredCut = 0.85,
  corFnc = bicor,
  t_CCRE_nasdevgeneexon, 
  corOptions = list(quick=0, use = "pairwise.complete.obs", maxPOutliers = .05),
  powerVector = c(1,11,21,26,28:30),
  networkType="signed",
  moreNetworkConcepts=T,
  verbose=0)
# save threshold selection
write.csv(threshold_biweight, file=file.path(newdir, "WGCNA_CCRE_biweight_power_selection.csv"))
# Power selection plot
power_1<-ggplot(threshold_biweight$fitIndices, aes(x=Power, y=-sign(slope)*SFT.R.sq))+geom_text(aes(label=Power))+geom_hline(yintercept=0.85)+theme_bw()
power_2<-ggplot(melt(threshold_biweight$fitIndices[,c(1,grep("k", names(threshold_biweight$fitIndices)))], id.vars = "Power"), aes(x=Power, col=variable, y=value, label=Power))+geom_text()+geom_line()+scale_y_log10()+theme_bw()

pdf(file=file.path(graphdir, "WGCNA_CCRE_biweight_Power_diagnostics.pdf"), paper = "a4r")
multiplot(power_1, power_2, cols = 2)
dev.off()

# Which threshold are we using?
powerthreshold <- ifelse(is.na(threshold_biweight$powerEstimate), 18, threshold_biweight$powerEstimate)
powerthreshold
# Calculate adjacency matrix
adj <- adjacency(t_CCRE_nasdevgeneexon, 
          power = powerthreshold, 
          corFnc = "bicor",
          corOptions = "quick=0, use = 'pairwise.complete.obs'",
          type = "signed")
save(adj, file = file.path(newdir, "/adj"), ascii = F)

# Convert to topological dissimilarity
dTOM <- TOMdist(adjMat = adj, TOMType = "signed")
row.names(dTOM) <- row.names(adj)
colnames(dTOM) <- colnames(adj)
# Save TOM matrix
save(dTOM, file = file.path(newdir, "/dTOM"), ascii = F)

#Calculate tree for clustering
geneTree <- flashClust(d = as.dist(dTOM), method = "average")
#Create clusters
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dTOM, 
                             deepSplit = 1, pamRespectsDendro = F, 
                             # maxCoreScatter = .7, minGap = .4,
                             minClusterSize = 20, minSplitHeight = .2)
# convert to numeric labels to colours
colclusterassignments<-labels2colors(dynamicMods, zeroIsGrey = T)
table(colclusterassignments)["grey"]
length(unique(colclusterassignments))
# Calculate module eigengenes
MEList<-moduleEigengenes(expr = t_CCRE_nasdevgeneexon, colors = colclusterassignments, 
                         grey = "grey", excludeGrey = F, nPC = 10)

# Savev module eigengenes
save(MEList, file = file.path(newdir, "CCRE_MEList.RData"))

# create and save data.frame with gene/cluster assignments
clusterassignments<-data.frame(transcriptID=colnames(t_CCRE_nasdevgeneexon), clusterID=colclusterassignments)
write.csv(clusterassignments, file=file.path(newdir, "CCRE_clusterassignments_simple.csv"))

# create and save data.frame with gene/cluster assigments and expression data
clusterassignments2<-data.frame(transposeBigData(t_CCRE_nasdevgeneexon), clusterID=colclusterassignments)
write.csv(clusterassignments2, file=file.path(newdir, "CCRE_clusterassignments_full.csv"))

# Cluster by eigengene dissimilarity
MEtree<-flashClust(as.dist(1-(bicor(MEList$eigengenes, quick=0, use = "pairwise.complete.obs"))), method = "average")
# Plot correlations
pdf(file = file.path(graphdir, "WGCNA_Cross_Cluster_Correlation.pdf"), width = 46.8, height = 33.1)
plot(MEtree, main="Cross-Cluster Correlations", xlab="", ylab="Proportion of Topological Overlap Divergence")
abline(a = 0.15, b=0, col="green")
dev.off()

# plot correlation heatmap
pdf(file = file.path(graphdir, "CCRE_Cluster_heatmaps.pdf"), width = 33, height = 33)
meltedcor <- MEList$eigengenes[,order.dendrogram(as.dendrogram(MEtree))]
meltedcor <- melt(bicor(meltedcor))
ggplot(data=meltedcor, aes(x=Var1, y=Var2, fill=value))+geom_tile()+scale_fill_gradient2()
dev.off()

# plot size distribution of clusters, highlighting grey
clustersizes <- as.data.frame(table(clusterassignments$clusterID))
names(clustersizes) <- c("clusterID","Size")
pdf(file = file.path(graphdir, "CCRE_Clustersizes.pdf"))
ggplot(clustersizes, aes(x = Size, col = clusterID=="grey")) + geom_density() + scale_x_log10() + theme_bw()
dev.off()

# plot ME by sample heatmap, clustered by tree similarity
meltedME <- MEList$eigengenes[,order.dendrogram(as.dendrogram(MEtree))]
row.names(meltedME) <- row.names(t_CCRE_nasdevgeneexon)
meltedME$sample <- factor(row.names(meltedME), levels = c("emb10_female.1", "emb10_female.2", "emb10_female.3", "emb10_male.1", "emb10_male.2", "emb10_male.3", "emb18_female.1", "emb18_female.2", "emb18_female.3", "emb18_male.1", "emb18_male.2", "emb18_male.3", "lar51_female.1", "lar51_female.2", "lar51_female.3", "lar51_male.1", "lar51_male.2", "lar51_male.3", "pupyel_female.1", "pupyel_female.2", "pupyel_female.3", "pupyel_male.1", "pupyel_male.2", "pupyel_male.3", "adult_female.1", "adult_female.2", "adult_female.3", "adult_male.1", "adult_male.2", "adult_male.3"))
meltedME <- melt(meltedME)
ggplot(data=meltedME, aes(x=sample, y=variable, fill=value))+geom_raster()+scale_fill_distiller(type = "seq", palette = 6)

### Outdated code
## Check distribution of different node types
# m_nasdevgeneexon<-data.frame(melt(t_CCRE_nasdevgeneexon),
#                              sampleID=row.names(t_CCRE_nasdevgeneexon))
# m_nasdevgeneexon$nodetype <- grepl(pattern = "CCRE", m_nasdevgeneexon$Var2)+grepl(pattern = "con", x = m_nasdevgeneexon$Var2)*2+grepl(pattern = "fac", x = m_nasdevgeneexon$Var2)*4
# ggplot(data = m_nasdevgeneexon, aes(x=value, col=as.factor(nodetype)))+geom_density()
# prop.table(table(m_nasdevgeneexon$nodetype)/10)zzzz

## QC and preprocessing
# # Detect noise quantile based on plot
# pdf(file = file.path(graphdir, "expr_threshold.pdf"))
# densityplot(as.numeric(t_CCRE_nasdevgeneexon), xlab = "Normalized Expression", panel=function(x,...){
#   panel.densityplot(x,...)
#   panel.abline(v=quantile(x,c(.2,.3,.4,.5)), col.line="red") 
#   panel.abline(v=quantile(x,c(.25,.35,.45,.55)), col.line="brown") 
# })
# dev.off()
# # Set data below noise threshold to threshold value
# noise <- quantile(as.numeric(t_CCRE_nasdevgeneexon), probs = .3)
# t_CCRE_nasdevgeneexon <- apply(t_CCRE_nasdevgeneexon, c(1,2), function(x){
#   ifelse(x < noise, noise, x)
# })
# # Rescale to 0-1 range
# t_CCRE_nasdevgeneexon <- apply(t_CCRE_nasdevgeneexon, c(1,2), function(x){x+abs(noise)}) # add minimum observed value
# upper <- max(as.numeric(t_CCRE_nasdevgeneexon))
# t_CCRE_nasdevgeneexon <- apply(t_CCRE_nasdevgeneexon, c(1,2), function(x){x/abs(upper)}) # divide by largest observed value
# # Apply logit to rescale to -inf +inf range
# t_CCRE_nasdevgeneexon
# # Check for zero variance nodes and nodes with too many missing values
# goodCols <- goodSamplesGenes(datExpr = t_CCRE_nasdevgeneexon, minNSamples = 3, minFraction = 1/10)
# lapply(X = goodCols, FUN = summary)
# # Remove genes below noise threshold
# t_CCRE_nasdevgeneexon <- t_CCRE_nasdevgeneexon[,which(goodCols$goodGenes==T)]
# # Save normalized and cleaned datasets
# save(t_CCRE_nasdevgeneexon, file = file.path(newdir, "t_CCRE_nasdevgeneexon_allOK.RData"))
# save(t(t_CCRE_nasdevgeneexon), file = file.path(newdir, "CCRE_nasdevgeneexon_allOK.RData"))

# # merge highly correlated modules (unnecessary since we work only on one block, replaced by option minAbsSplitHeight in cutreeDynamic)
# cut=0.15
# MElist<-mergeCloseModules(exprData = t_CCRE_nasdevgeneexon, colors = colclusterassignments, 
#                            corFnc = bicor, corOptions = list(quick=0, use="pairwise.complete.obs"), 
#                            cutHeight = cut, relabel = T,
#                            unassdColor = "grey", getNewMEs = T, getNewUnassdME = T)
# 
# colclusterassignments_2 <- MElist$colors