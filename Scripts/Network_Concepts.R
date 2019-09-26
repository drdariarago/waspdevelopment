## Calculate Fundamental Network concepts
## Initialize script
date()
rm(list=ls())
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(WGCNA)
allowWGCNAThreads()

# initialize output path
newdir<-file.path(getwd(), "Output/NetworkConcepts")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/NetworkConcepts")
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


## Calculate network concepts
emp_network_concepts <- networkConcepts(datExpr = t_nasdevgeneexon, power = powers$powerEstimate[1], networkType = "signed")

## Save as list
save(list = emp_network_concepts, file = file.path(newdir, "EMP_Network_Concepts.RData"))