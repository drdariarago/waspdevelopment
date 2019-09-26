### Create module metanetwork using CCRE data
# Initialize script
date()
rm(list=ls())
library(ggplot2)
library(reshape2)
library(flashClust)
library(WGCNA)
allowWGCNAThreads()
library(plyr)
library(stringr)
library(tnet)
source(file = "./Scripts/multiplot.R")

# initialize output path
newdir<-file.path(getwd(), "Output/Metanetwork_CCRE")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/Metanetwork_CCRE")
dir.create(graphdir)

# Load data
load(file = "./Output/WGCNA_CCRE_clustering_biweight/CCRE_MEList.RData")
# Create clustering tree from eigengene correlation matrix
MEtree<-flashClust(as.dist(
  1-(bicor(MEList$eigengenes, quick=0, use = "pairwise.complete.obs"))),
  method = "average")
OrderedMEs <- str_extract(string = names(MEList$eigengenes[,order.dendrogram(as.dendrogram(MEtree))]), pattern = "[^E]*$")

# generate correlation matrix 
cormat <- MEList$eigengenes[,order.dendrogram(as.dendrogram(MEtree))]
cormat <- bicor(MEList$eigengenes, quick = 0, use = "pairwise", maxPOutliers = "0.1")
colnames(cormat) <- substring(text = colnames(cormat), first = 3)
row.names(cormat) <- substring(text = row.names(cormat), first = 3)

# Convert to adjacency matrix
adjmat <- adjacency.fromSimilarity(cormat, type="signed")
# adjmat <- matrix(adjmat, nrow=153)
# diag(adjmat) <- 0
# row.names(adjmat) <- row.names(cormat)
# colnames(adjmat) <- colnames(cormat)

## Calculate cluster metanetwork attributes
emp_metanetworkconcepts <- fundamentalNetworkConcepts(adj = adjmat, GS = NULL)
## Metanetwork betweenness
adj_edgelist <- as.tnet(net = adjmat, type = "weighted one-mode tnet")
emp_metanetworkbetweenness <- betweenness_w(net = adj_edgelist, directed = F)
row.names(emp_metanetworkbetweenness) <- row.names(adjmat)
emp_metanetworkbetweenness <- as.data.frame(emp_metanetworkbetweenness)
## Scale betweenness by maximal betweenness possible
emp_metanetworkbetweenness$ScaledBetweenness <- with(emp_metanetworkbetweenness,
                                                     betweenness/((length(betweenness)-1)*(length(betweenness)-2)/2)
)
# Add betweenness to other concepts and save
emp_metanetworkconcepts$Betweenness <- emp_metanetworkbetweenness$betweenness
emp_metanetworkconcepts$ScaledBetweenness <- emp_metanetworkbetweenness$ScaledBetweenness
save(emp_metanetworkconcepts, file = file.path(newdir, "emp_metanetworkconcepts.RData"))

# Reshape as data.frame
data_emp_metanetwork <- data.frame(
  Connectivity = emp_metanetworkconcepts$Connectivity,
  ScaledConnectivity = emp_metanetworkconcepts$ScaledConnectivity, 
  ClusterCoef = emp_metanetworkconcepts$ClusterCoef, 
  MAR = emp_metanetworkconcepts$MAR, 
  Betweenness = emp_metanetworkconcepts$Betweenness,
  ScaledBetweenness = emp_metanetworkconcepts$ScaledBetweenness
)
write.csv(data_emp_metanetwork, file = file.path(newdir, "emp_metanetwork.csv"))

AbsMetanetwork <- data.frame(
  clusterID = row.names(data_emp_metanetwork),
  Connectivity = data_emp_metanetwork$Connectivity,
  Betweenness = data_emp_metanetwork$Betweenness,
  ClusterCoef = data_emp_metanetwork$ClusterCoef,
  MAR = data_emp_metanetwork$MAR,
  AbsConnectivity = data_emp_metanetwork$Connectivity/(length(data_emp_metanetwork$Connectivity)-1),
  AbsBetweenness = data_emp_metanetwork$Betweenness/((length(data_emp_metanetwork$Betweenness)-1)*(length(data_emp_metanetwork$Betweenness)-2)/2)
  #   MaxBetweenness = rep(length(data_emp_metanetwork$Betweenness)*(length(data_emp_metanetwork$Betweenness)-2)/2,length(data_emp_metanetwork$Betweenness)),
  #   MaxConnectivity = rep(length(data_emp_metanetwork$Connectivity),length(data_emp_metanetwork$Connectivity))
)

RelMetanetwork = data.frame(
  clusterID = AbsMetanetwork$clusterID,
  RelConnectivity = AbsMetanetwork$AbsConnectivity/max(AbsMetanetwork$AbsConnectivity),
  RelBetweenness = AbsMetanetwork$AbsBetweenness/max(AbsMetanetwork$AbsBetweenness),
  RelClusterCoef = data_emp_metanetwork$ClusterCoef/max(AbsMetanetwork$ClusterCoef),
  RelMAR = data_emp_metanetwork$MAR/max(AbsMetanetwork$MAR)
)
AllMetanetworks <- merge(AbsMetanetwork, RelMetanetwork, by = "clusterID")
summary(AllMetanetworks)
# write as .csv
write.csv(x = AllMetanetworks, file = file.path(newdir, "AllMetanetworks.csv"))

# Plot correlations between factors
pdf(file = file.path(graphdir, "ClusterProperties_correlations.pdf"))
pairs(data_emp_metanetwork) 
dev.off()

ggplot(data = data_emp_metanetwork, aes(x = ScaledConnectivity, y = ClusterCoef, col = MAR, size = log10(ScaledBetweenness))) + geom_point() + theme_bw() + scale_color_distiller(type = "div")

## Plot metanetwork
library(qgraph)

# Convert to adjacency matrix
adjmat <- matrix(adjmat, nrow=nrow(adjmat))
diag(adjmat) <- 0
row.names(adjmat) <- row.names(cormat)
colnames(adjmat) <- colnames(cormat)

# Plot metanetwork parameters (uses qgraph implementation, check how it handles weights)
pdf(file = file.path(newdir, "metanetwork_attributes.csv"))

meta_CtT <- centralityTable(adjmat, standardized = T)
meta_CtT$node <- factor(meta_CtT$node, levels = OrderedMEs)
meta_CtT <- meta_CtT[order(meta_CtT$node, meta_CtT$measure),]
ggplot(data = meta_CtT, aes(x = value, y = node, group = measure)) + geom_path() + geom_point() + 
  xlab("") + ylab("") + facet_grid( ~ measure, scales = "free")

meta_ClT <- clusteringTable(adjmat, standardized = T)
meta_ClT$node <- factor(meta_ClT$node, levels = OrderedMEs)
meta_ClT <- meta_ClT[order(meta_ClT$node, meta_ClT$measure),]
ggplot(data = meta_ClT, aes(x = value, y = node, group = measure)) + geom_path() + geom_point() + 
  xlab("") + ylab("") + facet_grid( ~ measure, scales = "free")

dev.off()

# Plot graph
qgraph(adjmat, layout = "spring")
# qgraph(adjmat, layout = "circle")
qgraph(adjmat, layout = layout.auto)

# Create significance graph (add sex bias later)
qgraph(input = cormat, diag = F, 
       minimum = 0.1, cut = 0.3, 
       vsize = 2, title = "Developmental Metanetwork", esize = 2,
       legend = T, borders = T,
       layout = "spring", graph = "sig2")

## Graph idea: superimpose groups of developmental stage bias and color nodes by sexbias
clusterdata <- read.csv(file = "./Output/Results_compiler_CCRE/clusterdata_full.csv")[,-1]

clustersexbias <- grepl(pattern = "m", x = clusterdata$cluster_devsexbias)+2*grepl(pattern = "f", x = clusterdata$cluster_devsexbias)
clustersexbias <- factor(clustersexbias, labels = c("unbiased", "malebiased", "femalebiased", "conflicting"))
sexshape <- factor(clustersexbias, labels = c("circle", "square", "triangle", "diamond"))

# Spring loaded layout
qgraph(input = cormat, diag = F, groups = clustersexbias,
       threshold = 0.1, cut = 0.3, maximum = 0.8,
       vsize = log10(clusterdata$clustersize), title = "Developmental Metanetwork", esize = 2,
       legend = T, borders = T, overlay = F,
       layout = "spring", graph = "sig2")
# Circular layouts
qgraph(input = cormat, diag = F, groups = clustersexbias,
       minimum = 0.1, cut = 0.3, maximum = 0.8,
       vsize = log10(clusterdata$clustersize), title = "Developmental Metanetwork", esize = 2,
       legend = T, borders = T, overlay = F,
       layout = "circular", graph = "sig2")

# Separate by stage-specific bias
qgraph(input = cormat, diag = F, groups = as.factor(clusterdata$cluster_devsexbias),
       threshold = 0.1, minimum = 0.3, cut = 0.7, maximum = 0.8, 
       shape = sexshape, vsize = log10(clusterdata$clustersize), 
       title = "Developmental Metanetwork", esize = 2,
       legend = T, borders = T, overlay = F,
       layout = "circular", graph = "sig2")

# Restrict to adult bias
adultbias <- factor(grepl(pattern = "m$", clusterdata$cluster_devsexbias)+2*grepl(pattern = "f$", clusterdata$cluster_devsexbias), 
                    labels = c("Unbiased", "Male_adult", "Female_adult"))
qgraph(input = cormat, diag = F, groups = adultbias,
       threshold = 0.1, minimum = 0.3, cut = 0.5, maximum = 0.8, 
       shape = sexshape, vsize = log10(clusterdata$clustersize), 
       title = "Developmental Metanetwork", esize = 2,
       legend = T, borders = T, overlay = F,
       layout = "circular", graph = "sig2")
# Create factor for : adult (m/f/u) and preadult (m/f/u)
preadultbias <- factor(grepl(pattern = "m.", clusterdata$cluster_devsexbias)+2*grepl(pattern = "f.", clusterdata$cluster_devsexbias), 
                       labels = c("Unbiased", "Male_preadult", "Female_preadult", "Confl_preadult"))
devbias <- ifelse(adultbias == "Unbiased", as.character(preadultbias), as.character(adultbias))
devbias <- ifelse(adultbias == "Female_adult" & preadultbias == "Male_preadult", "Confl_preadult", devbias)
devbias <- ifelse(adultbias == "Male_adult" & preadultbias == "Female_preadult", "Confl_preadult", devbias)
devbias <- factor(x = devbias, levels = c("Confl_preadult", "Female_preadult", "Female_adult", "Unbiased", "Male_adult", "Male_preadult"))

corgraph <- qgraph(input = cormat, diag = F, groups = as.factor(devbias),
                   threshold = 0.1, minimum = 0.3, cut = 0.5, maximum = 0.8, 
                   shape = sexshape, vsize = log10(clusterdata$clustersize), 
                   title = "Developmental Metanetwork", esize = 2,
                   legend = T, borders = T, overlay = F,
                   layout = "circular", graph = "sig2")

## IDEA: erase all connections that don't include at least one biased node and re-plot
coredges <- corgraph$Edgelist
sexcoredges <- which(coredges$from%in%which(devbias!="Unbiased")|coredges$to%in%which(devbias!="Unbiased"))
sexcoredges <- lapply(X = coredges, FUN = function(x){
  x[sexcoredges]
})
sexcoredges <- cbind(sexcoredges$from,
                     sexcoredges$to, 
                     sexcoredges$weight)
# Replace numbers with node IDs in the from and to columns
sexcoredges[,1:2] <- apply(X = sexcoredges[,1:2], MARGIN = c(1,2), FUN = function(x){
  clusterdata$clusterID[x]
})

sexcorenodes <- which(union(sexcoredges[,1], sexcoredges[,2])%in%clusterdata$clusterID)
sexparams <- list(
  devbias = devbias[sexcorenodes],
  devsexbias = clusterdata$cluster_devsexbias[sexcorenodes],
  sexshape = sexshape[sexcorenodes],
  # DIDC = paste(clusterdata$DI_sex, clusterdata$DI_stage, sep = "")[sexcorenodes],
  clustersize = clusterdata$clustersize[sexcorenodes])


qgraph(input = sexcoredges, groups = sexparams$devsexbias, 
       threshold = 0.5, minimum = 0.6, cut = 0.8, maximum = 1, 
       vsize = log10(sexparams$clustersize), shape = sexparams$sexshape,
       title = "Developmental Metanetwork", esize = 2,
       legend = T, borders = T, overlay = F,
       layout = "circular", directed = F, curve = .2, curveAll = T)

qgraph(input = sexcoredges, groups = sexparams$devsexbias, 
       threshold = 0.2, minimum = 0.3, cut = 0.6, maximum = 1, 
       vsize = log10(sexparams$clustersize)+1, shape = sexparams$sexshape,
       title = "Developmental Metanetwork", esize = 2,
       legend = T, borders = T, overlay = F,
       layout = "spring", directed = F, curve = .2, curveAll = T)

# # Sort clusters by stagebiased similarity
# stagebiasmatrix <- str_extract_all(string = clusterdata$cluster_stagebiased, pattern = ".", simplify = T)
# stagebiasmatrix <- matrix(as.numeric(as.factor(stagebiasmatrix)), ncol = 5)
# row.names(stagebiasmatrix) <- clusterdata$clusterID
# cmdscale(d = dist(stagebiasmatrix), k = 1)
# 
# library(psych)
# stagecormat <- cor(t(stagebiasmatrix))
# stagecormat <- apply(stagecormat, c(1,2), function(x) ifelse(is.na(x),0,x))
# stagepca <- principal(r = cormat, nfactors = 5, rotate = "promax")
# 
# qgraph(input = stagepca, groups = as.factor(clusterdata$cluster_devsexbias),
#        minimum = 0.4, cut = 0.6, maximum = 0.9, threshold = 0.3,
#        shape = sexshape, vsize = log10(clusterdata$clustersize),
#        title = "Developmental Metanetwork", vTrans = 200,
#        legend = T, borders = T, overlay = F,
#        layout = "circular", rotation = "promax")
# 
# stageK <- kmeans(x = stagebiasmatrix, centers = 10)
# a <-data.frame(clusterdata$cluster_stagebiased, stageK$cluster)
# a[order(a$stageK.cluster),]
# 
# corgraph <- qgraph(input = cormat, diag = F, groups = as.factor(stageK$cluster),
#                    threshold = 0.1, minimum = 0.3, cut = 0.7, maximum = 0.8, 
#                    shape = sexshape, vsize = log10(clusterdata$clustersize), 
#                    title = "Developmental Metanetwork", esize = 2,
#                    legend = T, borders = T, overlay = F,
#                    layout = "circular", graph = "sig2")
# 
# treepart <- cutree(tree = MEtree, k = 10)
# a <-data.frame(clusterdata$cluster_stagebiased, clusterdata$cluster_devsexbias, treepart)
# a[order(a$treepart),]
# corgraph <- qgraph(input = cormat, diag = F, groups = as.factor(treepart),
#                    threshold = 0.1, minimum = 0.3, cut = 0.7, maximum = 0.8, 
#                    shape = sexshape, vsize = log10(clusterdata$clustersize), 
#                    title = "Developmental Metanetwork", esize = 2,
#                    legend = T, borders = T, overlay = F,
#                    layout = "circular", graph = "sig2")


##### Outdated code snippets

# library(igraph)
# library(sna)
# library(ndtv)
# library(edgebundleR)
# library(HiveR)

# qrage(links = as.data.frame(as.edgelist.sna(adjmat)))

# Generate graph
# adjgraph <- graph.adjacency(adjmatrix = adjmat, weighted = T, diag = 0, mode = "undirected")
# qgraph(adjmat, layout = layout.graphopt, charge = edge.attributes(adjgraph)[[1]], threshold = 0.001)

# ########## HiveR section
# library(HiveR)
# diag(cormat) <- 0
# hiveadj <- adj2HPD(M = abs(cormat), type = "2D")
# hiveadj <- mineHPD(HPD = hiveadj, option = "rad <- tot.edge.count")
# hiveadj <- mineHPD(HPD = hiveadj, option = "remove orphans")
# hiveadj$edges$color <- ifelse( sign(cormat[which(lower.tri(cormat, diag = F))]) > 0, "blue", "orange")
# hiveadj$nodes$axis<-c(1:2)
# 
# sumHPD(hiveadj)
# plotHive(hiveadj)