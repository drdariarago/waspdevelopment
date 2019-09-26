# Calculate weighted betweenness 
### Calculate network concepts based on CCRE module network and CCRE network 
## Initialize script
date()
rm(list=ls())
library(plyr)
library(WGCNA)
allowWGCNAThreads()
library(igraph)
library(tnet)

# initialize output path
newdir<-file.path(getwd(), "Output/w_betweenness_CCREs")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/w_betweenness_CCREs")
dir.create(graphdir)

## Load data
load(file = "./Output/NetworkConcepts_CCRE/module_adj.RData")

# Convert to edgelists
edgelists <- lapply(l_adj, as.tnet, type = "weighted one-mode tnet")
# Calculate within block betweenness
module_emp_betweenness <- lapply(edgelists, betweenness_w, directed = F)
# Convert to data.frame and add node IDs
module_emp_betweenness <- ldply(module_emp_betweenness)
names(module_emp_betweenness) <- c("clusterID","nodeID","Betweenness")
module_emp_betweenness$nodeID <- unlist(sapply(X = l_adj, FUN = row.names))
# Save as csv
write.csv(module_emp_betweenness, file = file.path(newdir, "module_emp_betweenness.csv"))

## Calculate global betweenness 
# Load adjacency matrix and convert to edgelist
load(file = "./Output/WGCNA_CCRE_clustering_biweight/adj")
global_emp_betweenness <- as.tnet(net = adj, type = "weighted one-mode tnet")
# Calculate global betweenness
global_emp_betweenness <- betweenness_w(net = global_emp_betweenness, directed = F)
# Add node IDs
names(global_emp_betweenness) <- c("nodeID", "Betweenness")
global_emp_betweenness$betweenness <- row.names(adj)
# Save as csv
write.csv(global_emp_betweenness, file = file.path(newdir, "global_emp_betweenness.csv"))

#### Outdated code
## Load datasets
# load(file = "./Output/WGCNA_CCRE_clustering_biweight/adj")
# clusterassignments <- read.csv(file = "./Output/WGCNA_CCRE_clustering_biweight/CCRE_clusterassignments_simple.csv")[,-1]

## Edgelist conversion per module
# l_adj <- lapply(l_adj, function(x){
#   diag(x)<-0
# })
# edgelists <- lapply(l_adj, graph.adjacency, weighted = T)
# edgelists <- lapply(edgelists, get.data.frame)