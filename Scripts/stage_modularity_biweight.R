## Calculate integration and parcellation coefficients for stages
## Load libraries and clean space
date()
rm(list=ls())
library(plyr)
library(stringr)
library(WGCNA)
library(ggplot2)
allowWGCNAThreads()
# initialize output path
newdir<-file.path(getwd(), "Output/stage_modularity_biweight")
dir.create(newdir)

# import expression and cluster values
nasdevgeneexon<-read.csv(file="./Output/exon_to_transcript_clustering_rankingonly/eigenexon_evalues.csv")
nasdevclusters<-read.csv(file="./Output/WGCNA_clustering_kendall/clusterassignments_simple.csv")
# filter only transcripts present in the final network
nasdevgeneexon<-nasdevgeneexon[which(nasdevgeneexon$eigenexonID%in%nasdevclusters$transcriptID),]
# store transcript IDs in row.names
row.names(nasdevgeneexon)<-nasdevgeneexon$eigenexonID
nasdevgeneexon<-nasdevgeneexon[,-grep("eigenexonID", colnames(nasdevgeneexon))]
# reorder cluster IDs according to expression data
nasdevclusters<-nasdevclusters[match(row.names(nasdevgeneexon),nasdevclusters$transcriptID),]
# transpose expression data
t_nasdevgeneexon<-transposeBigData(nasdevgeneexon, blocksize = 50000)

# # limit to 3 clusters only
# nasdevclusters<-nasdevclusters[which(nasdevclusters$clusterID%in%(levels(nasdevclusters$clusterID)[7:9])),] 
# t_nasdevgeneexon<-t_nasdevgeneexon[,which(colnames(t_nasdevgeneexon)%in%nasdevclusters$transcriptID)]

# Create data.frame annotating each sample to respective stage
sampleIDframe<-data.frame(sampleID=names(nasdevgeneexon), stage=str_extract(string = names(nasdevgeneexon), pattern = "^[[:alnum:]]*"))
# Create balanced pairs of removed samples for each stage
combinations1<-ddply(.data = sampleIDframe, .variables = .(stage), .fun = function(x){
  expand.grid(grep("female", x$sampleID, value = T), grep("female", x$sampleID, value = T, invert = T))
})
## Combine stage-specific balanced pairs in all possible groups of 3 (creates matrix with one row per removed sample, one column per permutation)
# create all possible combinations of balanced samples
combinations2<-t(combn(1:nrow(combinations1), m = 3))
# retrieve sample names for each combination
combinations3<-t(apply(X = combinations2, MARGIN = 1, function(x){unlist(combinations1[x,2:3])}))
# remove combinations with duplicated samples (possible due to presence of same sample in different matched pairs)
combinations<-combinations3[which(apply(combinations3, 1, anyDuplicated)==0),]
# remove permutations of the same combination (sort alphabetically and then filter out duplicates)
combinations<-combinations[which(duplicated(apply(combinations, 1, function(x){paste(x[sort.list(x)], collapse="")}))==F),]

# add sexbias information
sexbias<-apply(combinations, 1, function(x){sum(grepl("female", x))})
# add stagebias information, 
stage<-factor(as.character(regmatches(names(nasdevgeneexon), gregexpr("^[[:alnum:]]+",names(nasdevgeneexon)))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
stagebias<-adply(combinations,1, function(z){
  aaply(levels(stage), 1, function(x) {
    sum(grepl(x,z))
  })
})
names(stagebias)<-c("X1", levels(stage))
combinations<-cbind(combinations, stagebias[,-1], sexbias)


# convert samples removed to character and store in single column
combinations[sapply(combinations, is.factor)]<-lapply(combinations[sapply(combinations, is.factor)], as.character)
combinations$samplesremoved<-paste0(combinations[1,1:6], collapse = ":")
# then grep each sample to remove against names(nasdevgeneexon) to get indices in data.frame
combinations[,1:6]<-apply(combinations[,1:6],c(1,2), function(x){grep(x, names(nasdevgeneexon))})
names(combinations)<-c(paste(rep("x",6), 1:6, sep=""), names(combinations[,-c(1:6)]))
# the final dataset contains the column number of samples removed in each combination (1:6), the stage-bias of samples removed (7:11), the sex bias of each combination (12:15) and the identity of the samples removed

# calculate connectivities for each dataset

# load power used for network construction
powers<-read.csv(file="./Output/WGCNA_clustering_biweight/WGCNA_biweight_power_selection.csv")

modularconnectivities<-alply(as.matrix(combinations[,1:6]), 1, function(expr){
  print(expr)
  intramodularConnectivity.fromExpr(datExpr = t_nasdevgeneexon[-c(expr),], colors = nasdevclusters$clusterID, 
                                    corFnc = "bicor", corOptions= "use='p', verbose='0', robustX=T, quick=5",
                                    distFnc = "dist", distOptions = "method='kendall'", 
                                    networkType = "signed", power = powers$powerEstimate[1],
                                    getWholeNetworkConnectivity = T)
})

names(modularconnectivities)<-row.names(combinations)
modularconnectivities<-ldply(modularconnectivities)
modularconnectivities$transcriptID<-colnames(t_nasdevgeneexon)

# save basic output
write.csv(modularconnectivities, file="./Output/stage_modularity_biweight/modularconnectivities_basic.csv")

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

# # plot changes in modularity parameters
# # reorder stages
# modulardensities$stage<-factor(modulardensities$stage, c("emb10","emb18", "lar51", "pupyel", "adult"))
# 
# # reshape dataset
# 
# library(reshape2)
# modulardensities2<-modulardensities[,c("clusterID","dWithin", "dOut", "emb10", "emb18", "lar51", "pupyel", "adult")]
# modulardensities2<-recast(modulardensities, formula = clusterID ~ dWithin + dOut)
# 
# # plot
# 
# pdf(file = file.path(newdir, "modularity_plots.pdf"), paper = "a4")
# 
# ggplot(modulardensities, aes(x=emb10, y=rank(dWithin)))+geom_point()+theme_bw()+geom_smooth(method="lm")+facet_wrap(~clusterID)+scale_color_brewer(type = "seq", palette = 4)+ggtitle("Integration vs Number of emb10 Samples Removed")
# ggplot(modulardensities, aes(x=emb10, y=rank(dOut)))+geom_point()+theme_bw()+geom_smooth(method="lm")+facet_wrap(~clusterID)+scale_color_brewer(type = "seq", palette = 4)+ggtitle("Pleiotropy vs Number of emb10 Samples Removed")
# ggplot(modulardensities, aes(x=rank(dWithin), y=rank(dOut)))+geom_point()+theme_bw()+geom_smooth(method="lm")+facet_wrap(~clusterID)+scale_color_brewer(type = "seq", palette = 4)+scale_x_log10()+scale_y_log10()+ggtitle("Integration vs Pleiotropy")
# 
# dev.off()


#### Outdated code snippets

# # TEMPORARY GENERATION OF RANDOM CLUSTERS, for testing only
# clusters<-data.frame(rpois(n=nrow(nasdevgeneexon),lambda = 4),row.names=row.names(nasdevgeneexon))
# colnames(clusters)<-"cluster"
# clusters$cluster[runif(n=nrow(clusters)/10, min=0, max=nrow(clusters))]<-NA
# clusters$geneID<-row.names(clusters)

## remove all possible combinations of 6 samples, then retain only those that have equal number of males and females for each stage
# # calculate all combinations
# combinations<-t(combn(names(nasdevgeneexon),6))
# 
# # get all male and female stage IDs
# malestages<-as.character(regmatches(names(nasdevgeneexon), gregexpr("^[[:alnum:]]+_male",names(nasdevgeneexon))))
# malestages<-malestages[-grep("char",malestages)]
# femalestages<-as.character(regmatches(names(nasdevgeneexon), gregexpr("^[[:alnum:]]+_female",names(nasdevgeneexon))))
# femalestages<-femalestages[-grep("char",femalestages)]
# sexstages<-unique(as.data.frame(cbind(malestages, femalestages)))
# 

# for every line of sexstages dataframe grep the first and the second column against every line of the combination data.frame and store the value of the stage-specific sex-unbalance
# 
# stagesexbias<-adply(combinations,1, function(z){
#   daply(sexstages, 1, function(x) {
#     sum(grepl(x$malestages,z))-sum(grepl(x$femalestages,z))
#   })
# })
# # remove all combinations where at least one stage is sex-unbalanced from stagesexbias table
# stagesexbias$toremove<-apply(stagesexbias[-1], 1, function(x){
#   any(x!=0)}
# )
# stagesexbias<-stagesexbias[-which(stagesexbias$toremove==T),-grep("toremove", names(stagesexbias))]
# and from table of sample IDs to remove in each combination
# combinations<-combinations[stagesexbias$X1,]


# # merge gene and exons with cluster identities (originally used to order before indexing extraction)
# nasclusters<-merge(nasdevgeneexon, clusters, by = "row.names", all.x=T)
# row.names(nasclusters)<-nasclusters$Row.names
# nasclusters<-nasclusters[,-1]
# # sort by cluster ID and gene ID, then remove cluster ID
# nasclusters<-nasclusters[order(nasclusters$cluster,row.names(nasclusters)),]
# nasclusters<-nasclusters[,-grep("cluster", names(nasclusters))]


# # original version, before serialization to list
# ddply(.data = clusters, .variables = .(cluster), .fun = summarize, 
#       integration = median(get("cor_matrices", env=.GlobalEnv)[[1]][x,y],na.rm = T),
#       parcellation = median(get("cor_matrices", env=.GlobalEnv)[[1]][x,-y],na.rm = T),
#       median_conectivity = median(get("cor_matrices", env=.GlobalEnv)[[1]],na.rm = T)
#       )


# # # Calculate correlation matrices for each combination 
# # # WARNING restricted to first 1000 genes (rows) for testing
# # cor_matrices<-alply(as.matrix(combinations[,1:6]),1,function(x){
# #   cor(t(nasdevgeneexon[,-c(x)]))}
# #   ,.dims = T) 
# 
# ## Store position of genes for each cluster in data.frame
# 
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
# write.csv(list = modularity, file="./Output/stage_modularity.csv")
# write.csv(list = modularity, file="./Output/stage_modularity/stage_modularity.csv")


## idea for generalized analysis of modularity in stage and stage:sex interaction
## generate random paired and unpaired permutations of samples 
## label with proportion of each stage (stage-bias) and within-stage sex-bias (stage-specific sex bias, or interaction term)
## use GLM to assess effect of stage bias and sex bias in module integration and parcellation
## problems: stage-specific networks remove six samples, stage:sex-specific ones remove three
## solution: add number of removed samples as additional nuisance factor (still collinear with sex-bias and stage-bias values)

# # troubleshooting version for permutation filtering
# # create all possible combinations of balanced samples
# combinations2<-t(combn(38:44, m = 3))
# # retrieve sample names for each combination
# combinations3<-t(apply(X = combinations2, MARGIN = 1, function(x){unlist(combinations1[x,2:3])}))
# # remove combinations with duplicated samples (possible due to presence of same sample in different matched pairs)
# combinations<-combinations3[which(apply(combinations3, 1, anyDuplicated)==0),]
# # remove permutations of the same combination (sort alphabetically and then filter out duplicates)
# combinations<-combinations[which(duplicated(apply(combinations, 1, function(x){paste(x[sort.list(x)], collapse="")}))==F),]
# 
