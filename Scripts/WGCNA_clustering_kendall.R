## Create modules blockwise
# load core packages
date()
rm(list=ls())
library(WGCNA)
allowWGCNAThreads()
# initialize output path
newdir<-file.path(getwd(), "Output/WGCNA_clustering_kendall")
dir.create(newdir)
# import dataset
nasdevgeneexon<-read.csv(file="./Output/exon_to_transcript_clustering_rankingonly/eigenexon_evalues.csv")
row.names(nasdevgeneexon)<-nasdevgeneexon$eigenexonID
# nasdevgeneexon<-nasdevgeneexon[1:1000,] # reduce size for testing
nasdevgeneexon<-nasdevgeneexon[,which(colnames(nasdevgeneexon)%in%setdiff(names(nasdevgeneexon),"eigenexonID"))]
t_nasdevgeneexon<-transposeBigData(nasdevgeneexon, blocksize = 50000)

# QC 
goodCols<-goodSamplesGenes(datExpr = t_nasdevgeneexon, minNSamples = 3)
lapply(goodCols, function(x) {table(x)})
lapply(goodCols, function(x) {prop.table(table(x))})
t_nasdevgeneexon<-t_nasdevgeneexon[,which(goodCols$goodGenes==T)]
# How many gene and how many transcript values were removed?
table(grepl("con", row.names(nasdevgeneexon)))-table(grepl("con", colnames(t_nasdevgeneexon)))

# calculate block size (check RAM availability trough unix command vmstat, RAM expressed in bytes)
bs<-blockSize(matrixSize = ncol(t_nasdevgeneexon), maxMemoryAllocation = 509435*1048576, overheadFactor = 4, rectangularBlocks = T)

# pick optimal power threshold with kendall distances
threshold_kendall<-pickSoftThreshold(
  blockSize = bs,
  t_nasdevgeneexon, 
  corOptions = list(method = "kendall", use = "pairwise.complete.obs"),
  powerVector = c(seq(from = 1,to = 14,by = 2),seq(from = 15,to = 21,by = 1),seq(from = 25, to = 30, by = 5)),
  networkType="signed",
  moreNetworkConcepts=T,
  verbose=5)
# save threshold selection
write.csv(threshold_kendall, file=file.path(newdir, "WGCNA_Kendall_power_selection.csv"))

# Power selection plot
library(ggplot2)
library(reshape2)
pdf(file=file.path(newdir, "WGCNA_Kendall_Power_diagnostics.pdf"))
ggplot(threshold_kendall$fitIndices, aes(x=Power, y=-sign(slope)*SFT.R.sq))+geom_text(aes(label=Power))+geom_hline(yintercept=0.9)+theme_bw()
ggplot(melt(threshold_kendall$fitIndices[,c(1,grep("k", names(threshold_kendall$fitIndices)))], id.vars = "Power"), aes(x=Power, col=variable, y=value, label=Power))+geom_text()+geom_line()+scale_y_log10()+theme_bw()
dev.off()

# Calculate block assignments
Kblocks<-projectiveKMeans(datExpr = t_nasdevgeneexon, preferredSize = bs, sizePenaltyPower = 2, networkType = "signed", checkData = T, verbose = 5)

# create list of sub-datasets for simplified calculations
Eval_blocks<-lapply(1:length(unique(Kblocks$clusters)), function(x){t_nasdevgeneexon[,which(Kblocks$clusters==x)]})
# validation of block assignments: subtracts indexes from blocks to indexes of matrices on main dataset (should return only zeroes)
# sapply(1:length(unique(Kblocks$clusters)), function(x){which(colnames(t_nasdevgeneexon)%in%colnames(Eval_blocks[[x]]))-which(Kblocks$clusters==x)})
all(unlist(sapply(1:length(unique(Kblocks$clusters)), function(x){which(colnames(t_nasdevgeneexon)%in%colnames(Eval_blocks[[x]]))-which(Kblocks$clusters==x)}))==0)

#Calculate adjacency matrices with chosen power level
# alternatively, use GTOMdist() (generalized topological overlap matrix, accounts for higher order interactions)
adj<-lapply(Eval_blocks, function(x){
  adjacency(x, 
            power = threshold_kendall$powerEstimate, 
            corOptions = list("method='kendall', use = 'pairwise.complete.obs'"), 
            type = "signed")
})
str(adj)

#Convert to topological dissimilarity matrix
dTOM<-lapply(adj, FUN = TOMdist, TOMType="signed")
for (block in 1:length(dTOM)){
  colnames(dTOM[[block]])<-colnames(adj[[block]])
  row.names(dTOM[[block]])<-row.names(adj[[block]])
}
str(dTOM)
# Save TOM matrix
save(dTOM, file = file.path(newdir, "/dTOM"), ascii = F)

#Calculate tree for clustering
geneTree<-lapply(X = dTOM, FUN = function(x){flashClust(as.dist(x), method="average")})

# plot clustering trees for each module
pdf(file=file.path(newdir, "WGCNA_Kendall_Cluster_trees.pdf"), paper = "a4")
lapply(X = geneTree, FUN = plot)
dev.off()

# Record gene names for each cluster
transcriptIDs<-lapply(dTOM, row.names)
#Create clusters
dynamicMods<-mapply(function(den, dis){cutreeDynamic(dendro = den, distM = dis, deepSplit = 4, pamRespectsDendro = F, minClusterSize = 20)}, den=geneTree, dis=dTOM)
str(dynamicMods)
# Generate unique labels for each cluster (add block ID after underscore)
clusterassignments<-mapply(function(x,y){paste0(x,"_",y)}, x=1:length(dynamicMods), y=dynamicMods)
# Fuse bloks and convert to numeric labels to colours
clusterassignments<-unlist(clusterassignments)
colclusterassignments<-labels2colors(clusterassignments)
# fuse cluster IDs 
transcriptIDs<-unlist(transcriptIDs)
# reorder expression data based on cluster IDs
t_nasdevgeneexon<-t_nasdevgeneexon[,match(transcriptIDs,colnames(t_nasdevgeneexon))]

# Calculate module eigengenes
MEList<-moduleEigengenes(expr = t_nasdevgeneexon, colors = colclusterassignments)
# Cluster by eigengene dissimilarity
MEtree<-flashClust(as.dist(1-cor(MEList$eigengenes, method = "kendall", use = "pairwise.complete.obs")), method = "average")
# Plot correlations
pdf(file = file.path(newdir, "WGCNA_Cross_Cluster_Correlation_beforemerging.pdf"), paper = "a4")
plot(MEtree, main="Cross-Cluster Correlations", xlab="", ylab="Proportion of Correlation Divergence")
abline(a = 0.1, b=0, col="green")
dev.off()

# merge highly correlated modules
MElist2<-mergeCloseModules(exprData = t_nasdevgeneexon, colors = colclusterassignments, 
                           corFnc = cor, corOptions = list(method="kendall", use="pairwise.complete.obs"), 
                           cutHeight = 0.1, relabel = F)
colclusterassignments2<-MElist2$colors

#Save workspace
save.image(file=file.path(newdir, "WGCNA_clustering_Kendall.RData"), ascii = F, compress = T, safe = T)

# create and save data.frame with gene/cluster assignments
clusterassignments<-data.frame(transcriptID=colnames(t_nasdevgeneexon), clusterID=colclusterassignments2)
write.csv(clusterassignments, file=file.path(newdir, "clusterassignments_simple.csv"))

# create and save data.frame with gene/cluster assigments and expression data
clusterassignments2<-data.frame(transposeBigData(t_nasdevgeneexon), clusterID=colclusterassignments2)
write.csv(clusterassignments2, file=file.path(newdir, "clusterassignments_full.csv"))

# plot the clusters before and after merging
pdf(file = file.path(newdir, "WGCNA_Cross_Cluster_Correlation_aftermerging.pdf"), paper = "a4")
plotDendroAndColors(flashClust(as.dist(1-TOMsimilarityFromExpr(t_nasdevgeneexon, power = threshold_kendall$powerEstimate, TOMType = "signed", corType = "bicor"))), 
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
system.time()