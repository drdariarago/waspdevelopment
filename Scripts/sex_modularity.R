## Calculate integration and parcellation coefficients for sexes
## Load libraries and clean space
date()
rm(list=ls())
library(plyr)
library(stringr)
library(WGCNA)
library(ggplot2)
allowWGCNAThreads()
# initialize output path
newdir<-file.path(getwd(), "Output/sex_modularity")
dir.create(newdir)

# import expression and cluster values
nasdevgeneexon<-read.csv(file="./Output/exon_to_transcript_clustering_rankingonly/eigenexon_evalues.csv")
nasdevclusters<-read.csv(file="./Output/WGCNA_clustering_kendall/clusterassignments_simple.csv")
# nasdevclusters<-nasdevclusters[1:500,]  # limit to first transcripts for testing
# nasdevclusters<-nasdevclusters[which(nasdevclusters$clusterID%in%(levels(nasdevclusters$clusterID)[1:9])),] # limit to first 3 clusters only
# filter only transcripts present in the final network
nasdevgeneexon<-nasdevgeneexon[which(nasdevgeneexon$eigenexonID%in%nasdevclusters$transcriptID),]
# store transcript IDs in row.names
row.names(nasdevgeneexon)<-nasdevgeneexon$eigenexonID
nasdevgeneexon<-nasdevgeneexon[,-grep("eigenexonID", colnames(nasdevgeneexon))]
# reorder cluster IDs according to expression data
nasdevclusters<-nasdevclusters[match(row.names(nasdevgeneexon),nasdevclusters$transcriptID),]
# transpose expression data
t_nasdevgeneexon<-transposeBigData(nasdevgeneexon, blocksize = 50000)


# # TEMPORARY GENERATION OF RANDOM CLUSTERS, for testing only
# clusters<-data.frame(rpois(n=nrow(nasdevgeneexon),lambda = 4),row.names=row.names(nasdevgeneexon))
# colnames(clusters)<-"cluster"
# clusters$cluster[runif(n=nrow(clusters)/10, min=0, max=nrow(clusters))]<-NA
# clusters$geneID<-row.names(clusters)

## remove random samples (three samples from a single stage)
# annotate stage in main environment
stage<-factor(as.character(str_extract(row.names(t_nasdevgeneexon),"^[[:alnum:]]+"),levels=c("emb10", "emb18", "lar51", "pupyel", "adult")))
# Create possible combinations for all stages and rename them
combinations<-lapply(levels(stage), function(x){combn(grep(x,row.names(t_nasdevgeneexon)),3)})
names(combinations)<-levels(stage)
# store combinations into data.frame
combinations<-as.data.frame(t(as.data.frame(combinations)))
# calculate degree of sex bias as number of female samples removed
combinations$sexbias<-apply(combinations, 1, function(x){sum(grepl("female", row.names(t_nasdevgeneexon[x,])))})
# annotate stage within combinations
combinations$stage<-factor(as.character(str_extract(row.names(combinations), "^[[:alnum:]]+"),levels=c("emb10", "emb18", "lar51", "pupyel", "adult")))


# calculate connectivities for each dataset

# load power used for network construction
powers<-read.csv(file="./Output/WGCNA_clustering_kendall/WGCNA_Kendall_power_selection.csv")

modularconnectivities<-alply(as.matrix(combinations[,1:3]), 1, function(expr){
  print(expr)
  intramodularConnectivity.fromExpr(datExpr = t_nasdevgeneexon[-c(expr),], colors = nasdevclusters$clusterID, 
                                    corFnc = "cor", corOptions= "use='p', method='kendall', verbose='4'",
                                    distFnc = "dist", distOptions = "method='kendall'", 
                                    networkType = "signed", power = powers$powerEstimate[1],
                                    getWholeNetworkConnectivity = T)
})
names(modularconnectivities)<-row.names(combinations)
modularconnectivities<-ldply(modularconnectivities)
modularconnectivities$transcriptID<-colnames(t_nasdevgeneexon)

# add cluster ID
modularconnectivities<-merge(modularconnectivities,nasdevclusters[,-1], by="transcriptID")
# add information on permutation
modularconnectivities<-merge(modularconnectivities,combinations,by.x="X1", by.y="row.names")
# add total network connectivity per permutation
kMain<-ddply(modularconnectivities, .(X1), numcolwise(median))[,c("X1","kTotal")]
names(kMain)<-c("X1", "kMain")
modularconnectivities<-merge(modularconnectivities, kMain, by.x="X1", by.y="X1", all.x=T)

# add cluster size
transcriptsWithin<-ddply(modularconnectivities, .(clusterID), nrow)
names(transcriptsWithin)<-c("clusterID", "clusterSize")
modularconnectivities<-merge(modularconnectivities, transcriptsWithin, by="clusterID", all.x=T)

## calculate density values for our variables for each gene*treatment
modularconnectivities$dWithin<-modularconnectivities$kWithin/(modularconnectivities$clusterSize-1)
modularconnectivities$dOut<-modularconnectivities$kOut/(nrow(modularconnectivities)-modularconnectivities$clusterSize)
modularconnectivities$dMain<-modularconnectivities$kMain/(nrow(modularconnectivities)-1)

# save dataset
write.csv(modularconnectivities, file=file.path(newdir, "modularconnectivities.csv"))

# summarize by module
modulardensities<-ddply(modularconnectivities, .(clusterID, stage, sexbias, X1), numcolwise(median))

# Save as csv
write.csv(modulardensities, file.path(newdir, "modulardensities.csv"))

# plot changes in modularity parameters
# reorder stages
modulardensities$stage<-factor(modulardensities$stage, c("emb10","emb18", "lar51", "pupyel", "adult"))

ggplot(modulardensities, aes(x=sexbias, y=dWithin, col=stage))+geom_point()+theme_bw()+geom_smooth(method="lm")+facet_grid(clusterID~.)+scale_color_brewer(type = "seq", palette = 4)
ggplot(modulardensities, aes(x=sexbias, y=dOut, col=stage))+geom_point()+theme_bw()+geom_smooth(method="lm")+facet_grid(clusterID~.)+scale_color_brewer(type = "seq", palette = 4)
ggplot(modulardensities, aes(x=dWithin, y=dOut, col=stage))+geom_point()+theme_bw()+geom_smooth(method="lm")+facet_grid(clusterID~.)+scale_color_brewer(type = "seq", palette = 4)+scale_x_log10()+scale_y_log10()

# #### OUTDATED code snippets

## Biweight multi-network density estimation
# modularconnectivities<-alply(as.matrix(combinations[,1:3]), 1, function(expr){
#   print(expr)
#   intramodularConnectivity.fromExpr(datExpr = t_nasdevgeneexon[-c(expr),], colors = nasdevclusters$clusterID, 
#                                     corFnc = "bicor", corOptions= "use='p', robustX=T, quick=0, pearsonFallback='individual', verbose=2",
#                                     distFnc = "dist", distOptions = "method='kendall'", 
#                                     networkType = "signed", power = powers$powerEstimate[1],
#                                     getWholeNetworkConnectivity = T)
# })


# # single dataset version for calculation of network values
# intramodularConnectivity.fromExpr(datExpr = nasdevgeneexon, colors = nasdevclusters$clusterID, 
#                                   corFnc = "cor", corOptions = list(method="kendall", use="pairwise.complete.obs"), 
#                                   distFnc = "dist", distOptions = list(method="kendall"), 
#                                   networkType = "signed", power = 16,
#                                   getWholeNetworkConnectivity = T)

# 
# 
# # Calculate correlation matrices for each combination 
# # WARNING restricted to first 1000 genes (rows) for testing
# cor_matrices<-alply(as.matrix(combinations[,1:3]),1,function(x){
#   cor(t(nasdevgeneexon[,-c(x)]), method = "kendall", use = "pairwise.complete.obs")}
#   ,.dims = T) 
# 
# ## Store position of genes for each cluster in data.frame
# # ## WARGING re-generating random labels for clusters and removing genes missing from the top 1000 used for matrices, testing ONLY
# # clusters$cluster<-rpois(n=nrow(clusters),lambda = 4)
# # clusters$cluster[runif(n=nrow(clusters)/10, min=0, max=nrow(clusters))]<-NA
# # clusters<-clusters[which(clusters$geneID%in%row.names(cor_matrices[[1]])),]
# 
# # initializing variables for row-col position and index
# clusters$x<-NA
# clusters$y<-NA
# index=0
# # recording position of clusters (pasting end-of-line after geneID to avoid confusion with exons)
# for(i in clusters$geneID) {
#   index<-index+1
#   clusters$x[index]<-grep(paste0(i,"$"), row.names(cor_matrices[[1]]))
#   clusters$y[index]<-grep(paste0(i,"$"), row.names(cor_matrices[[1]]))
# }
# rm(index)
# 
# ### Calculate integration and parcellation coefficients for each cluster
# modularity=list()
# # WARNING: reduce loop length for testing
# for (i in 1:nrow(combinations)){
#   print(i)
#   z<-i
#   tempdata<-ddply(.data = clusters, .variables = .(cluster), .fun = summarize, .inform = T,
#                   combination=names(get("cor_matrices", env=.GlobalEnv)[z]),
#                   integration = median(get("cor_matrices", env=.GlobalEnv)[[z]][x,y],na.rm = T),
#                   parcellation = median(get("cor_matrices", env=.GlobalEnv)[[z]][x,-y],na.rm = T),
#                   median_conectivity = median(get("cor_matrices", env=.GlobalEnv)[[z]],na.rm = T))
#   modularity[[names(get("cor_matrices", env=.GlobalEnv)[z])]]<-tempdata
#   rm(tempdata)
# }
# rm(z)
# # stack combinations in single data.frame
# modularity<-ldply(modularity)
# # add sex bias and stage data
# modularity<-merge(modularity, combinations[,4:5], by.y="row.names", by.x="combination")
# # and save as .csv
# write.csv(list = modularity, file="./Output/sex_modularity.csv")
# write.csv(list = modularity, file="./Output/sex_modularity/sex_modularity.csv")
# 
# 
# # # merge gene and exons with cluster identities (originally used to order before indexing extraction)
# # nasclusters<-merge(nasdevgeneexon, clusters, by = "row.names", all.x=T)
# # row.names(nasclusters)<-nasclusters$Row.names
# # nasclusters<-nasclusters[,-1]
# # # sort by cluster ID and gene ID, then remove cluster ID
# # nasclusters<-nasclusters[order(nasclusters$cluster,row.names(nasclusters)),]
# # nasclusters<-nasclusters[,-grep("cluster", names(nasclusters))]
# 
# 
# # # original version, before serialization to list
# # ddply(.data = clusters, .variables = .(cluster), .fun = summarize, 
# #       integration = median(get("cor_matrices", env=.GlobalEnv)[[1]][x,y],na.rm = T),
# #       parcellation = median(get("cor_matrices", env=.GlobalEnv)[[1]][x,-y],na.rm = T),
# #       median_conectivity = median(get("cor_matrices", env=.GlobalEnv)[[1]],na.rm = T)
# #       )
# 
# 
# ## calculate density values for our variables for each cluster*treatment
# # function to calculate number of pairwise interactions in a set
# handshakes<-function(x){
#   y<-x-1
#   (y*(y+1))/2
# }
# modulardensities$dWithin<-modulardensities$kWithin/handshakes(modulardensities$clusterSize)
# modulardensities$dOut<-modulardensities$kOut/(modulardensities$clusterSize*(nrow(modularconnectivities)-modulardensities$clusterSize))
# # modulardensities$dTotal<-modulardensities$kTotal/handshakes(nrow(modularconnectivities)) # this version calculates total density of each cluster with every node of the network
# modulardensities$dMain<-modulardensities$kMain/handshakes(nrow(modularconnectivities))
# 
# 
# ## idea for generalized analysis of modularity in stage and stage:sex interaction
# ## generate random paired and unpaired permutations of samples 
# ## label with proportion of each stage (stage-bias) and within-stage sex-bias (stage-specific sex bias, or interaction term)
# ## use GLM to assess effect of stage bias and sex bias in module integration and parcellation
# ## problems: stage-specific networks remove six samples, stage:sex-specific ones remove three
# ## solution: add number of removed samples as additional nuisance factor (still collinear with sex-bias and stage-bias values)