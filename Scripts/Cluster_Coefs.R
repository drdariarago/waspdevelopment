# Calculate clustering coefficient distribution for global network and modules
# then compare with networks using degree-preserving redistribution of links

## Initialize script
date()
rm(list=ls())
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(WGCNA)
library(fdrtool)
library(tnet)
# initialize output path
newdir<-file.path(getwd(), "Output/Cluster_Coefs")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/Cluster_Coefs")
dir.create(graphdir)

## Load data
powers<-read.csv(file="./Output/WGCNA_clustering_biweight/WGCNA_biweight_power_selection.csv")
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

## Calcuate global and intramodular connectivities
emp_modularconnectivities<- intramodularConnectivity.fromExpr(datExpr = t_nasdevgeneexon, colors = nasdevclusters$clusterID,
                                    ignoreColors = NULL, scaleByMax = F,
                                    corFnc = "bicor", corOptions= "use='p', verbose='4', robustX=T, quick=0",
                                    distFnc = "dist", distOptions = "method='kendall'", 
                                    networkType = "signed", power = powers$powerEstimate[1],
                                    getWholeNetworkConnectivity = T)
# Save as csv
write.csv(emp_modularconnectivities, file = file.path(newdir, "emp_modularconnectivities.csv"))

# Calculate node cluster-specificity (module connectivity vs global connectivity)
emp_modularconnectivities$kWithin*(1/emp_modularconnectivities$kTotal)

# load adjacency matrix
load(file = "./Output/WGCNA_clustering_biweight/adj")

# Calculate empirical network concepts (esp clusterCoef)
emp_networkconcepts <- fundamentalNetworkConcepts(adjMat = adj)
# Save cluster coefficients
save(list = "emp_networkconcepts", file = file.path(newdir, "emp_networkconcepts.RData"))





# Plot raw distributions


# Randomize network
# set seeds for reproducibility
sample(x = 1:10000,size = 500, replace = F)
print(sample)
# convert adjacency matrix to weighted edgelist
# creating adj from sample data to test
adj <- adjacency(t_nasdevgeneexon[,1:2000])


# Reshuffle
rg_reshuffling_w(net = vectorizeMatrix(adj), option = "links")
