## Calculate Fundamental Network concepts using adjacency matrix

## Initialize script
date()
rm(list=ls())
library(WGCNA)
allowWGCNAThreads()

# initialize output path
newdir<-file.path(getwd(), "Output/NetworkConcepts")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/NetworkConcepts")
dir.create(graphdir)

## Load data
load(file = "./Output/WGCNA_clustering_biweight/adj")

## Merge the two matrices
# Add the extra rows to each matrix

adj2 <- list(NA,NA)

adj2[[1]] <- rbind(
  adj[[1]], matrix(0, nrow = nrow(adj[[2]]), ncol = ncol(adj[[1]]))
)

adj2[[2]] <- rbind(
  matrix(0, nrow = nrow(adj[[1]]), ncol = ncol(adj[[2]])), adj[[2]]
)
# Then merge matrices as columns
adj2 <- cbind(adj2[[1]], adj2[[2]])
row.names(adj2)<-c(row.names(adj[[1]]),row.names(adj[[2]]))

# Save merged adjacency matrix
save(adj2, file = file.path(newdir, "adj2_merged"))

## Calculate network concepts
emp_network_concepts <- fundamentalNetworkConcepts(adj = adj2, GS = NULL)

ls()
warnings()

## Save as list
save(list = emp_network_concepts, file = file.path(newdir, "EMP_Network_Concepts_adj.RData"))