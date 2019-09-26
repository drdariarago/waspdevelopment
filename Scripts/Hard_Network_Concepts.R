## Calculate Fundamental Network concepts based on binary network (hard threshold)
## Initialize script
date()
rm(list=ls())
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(WGCNA)
allowWGCNAThreads()

source(file = "./Scripts/multiplot.R")

logit<-function(x){
  log10(x/(1-x))
}

z<-function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
}


# initialize output path
newdir<-file.path(getwd(), "Output/HardNetworkConcepts")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/HardNetworkConcepts")
dir.create(graphdir)

## Load data
nasdevgeneexon<-read.csv(file="./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexon_evalues.csv")
nasdevclusters<-read.csv(file="./Output/WGCNA_clustering_biweight/clusterassignments_simple.csv")
nasdevgeneexon<-nasdevgeneexon[which(nasdevgeneexon$eigenexonID%in%nasdevclusters$transcriptID),]

## Prepare dataset for network calculations
# filter only transcripts present in the final network
nasdevgeneexon<-nasdevgeneexon[which(nasdevgeneexon$eigenexonID%in%nasdevclusters$transcriptID),]
# store transcript IDs in row.names
row.names(nasdevgeneexon)<-nasdevgeneexon$eigenexonID
nasdevgeneexon<-nasdevgeneexon[,-grep("eigenexonID", colnames(nasdevgeneexon))]
# reorder cluster IDs according to expression data
nasdevclusters<-nasdevclusters[match(row.names(nasdevgeneexon),nasdevclusters$transcriptID),]
# transpose expression data
t_nasdevgeneexon<-transposeBigData(nasdevgeneexon, blocksize = 50000)

# QC 
# Check for zero variance nodes and nodes with too many missing values
goodCols<-goodSamplesGenes(datExpr = t_nasdevgeneexon, minNSamples = 3, minFraction = 1/10)
lapply(goodCols, function(x) {table(x)})
lapply(goodCols, function(x) {prop.table(table(x))})
t_nasdevgeneexon[,which(goodCols$goodGenes==F)]
t_nasdevgeneexon<-t_nasdevgeneexon[,which(goodCols$goodGenes==T)]

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

### Split into blocks for correlation calculations
# calculate block size (check RAM availability trough unix command vmstat, RAM expressed in bytes)
bs<-blockSize(matrixSize = ncol(t_nasdevgeneexon), maxMemoryAllocation = 509435*1048576, overheadFactor = 4, rectangularBlocks = T)
# Block size is hard limited to a maximum of 46340
bs<-ifelse(bs>46339, 46339, bs)
bs
# Calculate block assignments
Kblocks<-projectiveKMeans(datExpr = t_nasdevgeneexon, preferredSize = bs, sizePenaltyPower = 2, networkType = "signed", checkData = T, verbose = 0)
# create list of sub-datasets for simplified calculations
Eval_blocks<-lapply(1:length(unique(Kblocks$clusters)), function(x){t_nasdevgeneexon[,which(Kblocks$clusters==x)]})

### Crete correlation matrix with bicor function 
cormat <- lapply(X = Eval_blocks, FUN = bicor,  quick = 0, use = "pairwise.complete.obs")
str(cormat)

## Assemble blocks (assumes 2 blocks only)
cormat2 <- list(NA,NA)
cormat2[[1]] <- rbind(
  cormat[[1]], matrix(NA, nrow = nrow(cormat[[2]]), ncol = ncol(cormat[[1]]))
)
cormat2[[2]] <- rbind(
  matrix(NA, nrow = nrow(cormat[[1]]), ncol = ncol(cormat[[2]])), cormat[[2]]
)
## Then merge matrices as columns
cormat2 <- cbind(cormat2[[1]], cormat2[[2]])
row.names(cormat2)<-c(row.names(cormat[[1]]),row.names(cormat[[2]]))
## Save correlation matrix
save(cormat2, file=file.path(newdir,"cormat2.RData"))

### pick optimal threshold with biweight midcorrelations
# Rescale correlations to 0-1 range
cormat01 <- (cormat+1)/2
cormat01[which(is.na(cormat01))] <- 0.5
# Calculate threshold
threshold_biweight <- pickHardThreshold.fromSimilarity(
  cormat2,
  RsquaredCut = 0.85,
  cutVector = seq( from = 0.1, to = 0.9, by = 0.05),
  moreNetworkConcepts = TRUE
)
# save threshold selection
write.csv(threshold_biweight, file=file.path(newdir, "WGCNA_biweight_threshold_selection.csv"))

###Calculate adjacency matrices with chosen threshold
cormat2[which(is.na(cormat2))] <- 0

signedAdjacencyFunction <- function(corMat, threshold){
  adjmat <- as.matrix(abs(corMat) >= threshold)
  adjmat <- adjmat*sign(corMat)
  dimnames(adjmat) <- dimnames(corMat)
  diag(adjmat) <- 0
  adjmat
}

adj <- signedAdjacencyFunction(
  corMat = cormat2,
  threshold = threshold_biweight$cutEstimate
)

str(adj)
save(adj, file = file.path(newdir, "/hard_adj"), ascii = F)

### Calculate network concepts