## Analyze network structure of whole network and clusters (divide between NDR/DE/DI)

## Initialize script
date()
rm(list=ls())
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(WGCNA)
allowWGCNAThreads()
## initialize output path
newdir<-file.path(getwd(), "Output/NetworkConcepts_analysis")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/NetworkConcepts_analysis")
dir.create(graphdir)

## Load dataset
load("./Output/NetworkConcepts/EMP_Network_Concepts.RData")

## Compile node data (add DE/DI/Cluster/Transcript information)
clusterdata <- read.csv(file = "./Output/Results_compiler/clusterdata_full.csv")[,-1]

## Plot whole network distributions of Connectivity vs Clustering vs Connection Strength

## Split by Transcription/Splicing

## Split by Cluster

## Find hubs in each cluster

## Annotate clusters 