## Detect hubs based on connectivity, waypoints based on betweenness and bottlenecks based on betweenness/connectivity
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
newdir<-file.path(getwd(), "Output/hub_characterization")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/hub_characterization")
dir.create(graphdir)

## Load data
nodedata <- read.csv("./Output/Results_compiler_CCRE/nodedata_full.csv")[,-1]
nodedata <- droplevels(nodedata[-which(is.na(nodedata$geneID)),])
nodedata$Grey <- as.factor(ifelse(nodedata$clusterID=="grey","grey","cluster"))
load(file = "./Output/WGCNA_CCRE_clustering_biweight/CCRE_MEList.RData")
eigenexons <- MEList$eigengenes

### Calculate cluster hubs
# Restrict to core and exclude grey
cores <- nodedata[which(nodedata$KME>0.5&nodedata$clusterID!="grey"),]
# normalize by KME, 
ScaledConnLM <- lm(formula = log10(ScaledConnectivity_WithinCluster) ~ KME, data = cores)
# plot(ScaledConnLM)
summary(ScaledConnLM)
cores$ScaledConnectivity_WithinResiduals <- resid(ScaledConnLM)
# top 5 nodes per cluster
cores <- dlply(.data = cores, .variables = .(clusterID))
cores <- lapply(X = cores, FUN = function(x){x[order(x$ScaledConnectivity_WithinResiduals, decreasing = T, na.last = T),][1:5,]})
cores <- ldply(cores)
head(cores)

