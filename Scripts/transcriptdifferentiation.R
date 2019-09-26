## Plot relationship between of exons/transcripts/clusters
# Initialize script
library(ggplot2)
library(plyr)
rm(list=ls())
newdir<-file.path(getwd(),"Output/transcriptdifferentiation")
dir.create(newdir)

# Load data
eigenexon_data<-read.csv(file="./Output/Results_compiler/eigenexon_data.csv")

# Plot n exons vs n transcripts

# Plot n transcripts vs n clusters
trans_clusters<-ddply(.data = eigenexon_data, .variables = .(geneID), summarize, nEigenexons=length(unique(eigenexonID)), nClusters=length(unique(clusterID)), .progress="text")
ggplot(data = eigenexon_data, aes(x=))

# Plot 