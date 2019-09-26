### Error report when calculating adjacency matrix with biweight midcorrelation on BEAR

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
goodCols<-apply(t_nasdevgeneexon, 2, function(x){sum(x>0)>=3})
t_nasdevgeneexon<-t_nasdevgeneexon[,which(goodCols==T)]
# Check for zero variance nodes
goodCols<-goodSamplesGenes(datExpr = t_nasdevgeneexon, minNSamples = 3, minFraction = 1/10)
lapply(goodCols, function(x) {table(x)})
lapply(goodCols, function(x) {prop.table(table(x))})
t_nasdevgeneexon<-t_nasdevgeneexon[,which(goodCols$goodGenes==T)]
# How many gene and how many transcript values were removed?
table(grepl("con", row.names(nasdevgeneexon)))-table(grepl("con", colnames(t_nasdevgeneexon)))

str(t_nasdevgeneexon)
summary(apply(X = t_nasdevgeneexon, MARGIN = c(1,2), FUN = function(x)(x!=0)))

# # calculate block size (check RAM availability trough unix command vmstat, RAM expressed in bytes)
# bs<-blockSize(matrixSize = ncol(t_nasdevgeneexon), maxMemoryAllocation = 509435*1048576, overheadFactor = 4, rectangularBlocks = T)
# bs
# 
# # pick optimal power threshold with biweight midcorrelations
# threshold_biweight<-pickSoftThreshold(
#   corFnc = bicor,
#   blockSize = bs,
#   t_nasdevgeneexon, 
#   corOptions = list(quick=0.3, use = "pairwise.complete.obs"),
#   powerVector = c(seq(from = 1,to = 21,by = 10),seq(from = 22, to = 31, by = 3),31:35),
#   networkType="signed",
#   moreNetworkConcepts=T,
#   verbose=1)
# # save threshold selection
# write.csv(threshold_biweight, file=file.path(newdir, "WGCNA_biweight_power_selection.csv"))
# 
# # Power selection plot
# pdf(file=file.path(graphdir, "WGCNA_biweight_Power_diagnostics.pdf"))
# ggplot(threshold_biweight$fitIndices, aes(x=Power, y=-sign(slope)*SFT.R.sq))+geom_text(aes(label=Power))+geom_hline(yintercept=0.9)+theme_bw()
# ggplot(melt(threshold_biweight$fitIndices[,c(1,grep("k", names(threshold_biweight$fitIndices)))], id.vars = "Power"), aes(x=Power, col=variable, y=value, label=Power))+geom_text()+geom_line()+scale_y_log10()+theme_bw()
# dev.off()

#Calculate adjacency matrix with chosen power level
# alternatively, use GTOMdist() (generalized topological overlap matrix, accounts for higher order interactions)
adj<-adjacency(datExpr = t_nasdevgeneexon, 
               power = 33, 
               corFnc = "bicor",
               corOptions = "quick=0, use = 'pairwise.complete.obs'",
               type = "signed")
str(adj)
summary(apply(X = adj, MARGIN = c(1,2), FUN = function(x)(x!=0)))