### Calculate network concepts based on CCRE module network and CCRE network 
## Initialize script
date()
rm(list=ls())
library(plyr)
library(WGCNA)
allowWGCNAThreads()

# initialize output path
newdir<-file.path(getwd(), "Output/NetworkConcepts_CCRE")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/NetworkConcepts_CCRE")
dir.create(graphdir)

## Load data
load(file = "./Output/WGCNA_CCRE_clustering_biweight/adj")
clusterassignments <- read.csv(file = "./Output/WGCNA_CCRE_clustering_biweight/CCRE_clusterassignments_simple.csv")[,-1]

## Replace missing values with zeroes
adj[which(is.na(adj))]<-0

## Split into subnetworks (one for each module)
l_adj <- split(x = clusterassignments$transcriptID, f = clusterassignments$clusterID, drop = T)
l_adj <- lapply(X = l_adj, FUN = function(nodeID){
  adj[which(row.names(adj)%in%nodeID), which(colnames(adj)%in%nodeID)]
})

# # Remove grey
# l_adj <- l_adj[-which(names(l_adj)=="grey")]

# Save list of matrices
save(l_adj, file = file.path(newdir, "module_adj.RData"))

# Calculate network concepts within each
module_emp_network_concepts <- lapply(X = l_adj, FUN = fundamentalNetworkConcepts)

# Save network concept list
save(module_emp_network_concepts, file = file.path(newdir, "cluster_emp_networkconcepts.RData"))

# Convert node parameters to data.frame and save as csv
nodeparams <- lapply(module_emp_network_concepts, function(moduledata){
  data.frame(
    eigenexonID = names(moduledata[[1]]),
    moduledata[1:4])
})
nodeparams <- ldply(nodeparams, .id = "clusterID")
write.csv(nodeparams, file = file.path(newdir, "CCRE_nodeparams.csv"))

# Convert cluster parameters to data.frame and save as csv
clusterparams <- lapply(module_emp_network_concepts, function(moduledata){
  as.data.frame(moduledata[5:7])
})
clusterparams <- ldply(clusterparams)
names(clusterparams)[1] <- "clusterID"
write.csv(clusterparams, file = file.path(newdir, "CCRE_clusterparams.csv"))

## Calculate global network concepts
emp_network_concepts <- fundamentalNetworkConcepts(adj = adj, GS = NULL)
str(emp_network_concepts)
dim(emp_network_concepts)
## Save as list
save(emp_network_concepts, file = file.path(newdir, "EMP_Network_Concepts_adj_CCRE.RData"))

# Convert and save as data.frames
global_nodeparams <- data.frame(
  # nodeID = names(emp_network_concepts[[1]]),
  Connectivity = emp_network_concepts[[1]],
  ScaledConnectivity = emp_network_concepts[[2]],
  ClusterCoef = emp_network_concepts[[3]],
  MAR = emp_network_concepts[[4]]
)
  
write.csv(global_nodeparams, file = file.path(newdir, "CCRE_nodeparams_global.csv"))

# Convert cluster parameters to data.frame and save as csv
global_networkparams <- as.data.frame(emp_network_concepts[5:7])
write.csv(global_networkparams, file = file.path(newdir, "CCRE_networkparams.csv"))