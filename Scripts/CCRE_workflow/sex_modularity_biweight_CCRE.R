## Calculate integration and parcellation coefficients for sexes in CCRE network
## Load libraries and clean space
date()
rm(list=ls())
library(plyr)
library(stringr)
library(WGCNA)
library(ggplot2)
allowWGCNAThreads()

# Function to calculate maximum pairwise interactions within a cluster
# Double the number of handshakes since we're adding up all nodes in the network
maxConn<-function(x){
  x*(x-1)
}

# initialize output path
newdir<-file.path(getwd(), "Output/sex_modularity_biweight_CCRE")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/sex_modularity_biweight_CCRE")
dir.create(graphdir)

# import expression and cluster values
nasdevclusters <- read.csv(file="./Output/WGCNA_CCRE_clustering_biweight/CCRE_clusterassignments_simple.csv")[,-1]
nasdevgeneexon <- read.csv(file="./Output/CCREs/CCRE_nasdevgeneexon.csv")
row.names(nasdevgeneexon) <- nasdevgeneexon[,1]
nasdevgeneexon <- nasdevgeneexon[,-1]
t_nasdevgeneexon <- t(nasdevgeneexon)
# t_nasdevgeneexon <- t_nasdevgeneexon[,1:1000]

# # store transcript IDs in row.names
# row.names(nasdevgeneexon)<-nasdevgeneexon$eigenexonID
# nasdevgeneexon<-nasdevgeneexon[,-grep("eigenexonID", colnames(nasdevgeneexon))]
# # reorder cluster IDs according to expression data
# nasdevclusters<-nasdevclusters[match(row.names(nasdevgeneexon),nasdevclusters$transcriptID),]
# # transpose expression data
# t_nasdevgeneexon<-transposeBigData(nasdevgeneexon, blocksize = 50000)


# reorder cluster IDs according to expression data
nasdevclusters<-nasdevclusters[match(colnames(t_nasdevgeneexon),nasdevclusters$transcriptID),]
#####
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


## load power used for network construction
powers<-read.csv(file="./Output/WGCNA_CCRE_clustering_biweight/WGCNA_CCRE_biweight_power_selection.csv")
## calculate connectivities for each dataset
modularconnectivities<-alply(as.matrix(combinations[,1:3]), 1, function(expr){
  #   print(expr)
  intramodularConnectivity.fromExpr(datExpr = t_nasdevgeneexon[-c(expr),], colors = nasdevclusters$clusterID,
                                    ignoreColors = NULL, scaleByMax = F,
                                    corFnc = "bicor", corOptions= "use='p', verbose='4', robustX=T, quick=0",
                                    distFnc = "dist", distOptions = "method='kendall'", 
                                    networkType = "signed", power = powers$powerEstimate[1],
                                    getWholeNetworkConnectivity = T)
})
names(modularconnectivities)<-row.names(combinations)
modularconnectivities<-ldply(modularconnectivities)
modularconnectivities$transcriptID<-colnames(t_nasdevgeneexon)

# add cluster ID
modularconnectivities<-merge(modularconnectivities,nasdevclusters, by="transcriptID") # line not working
# Replace connectivities lesser than zero with zero (they arise either from genes with too few expression values)
modularconnectivities$kTotal<-ifelse(modularconnectivities$kTotal<0,0,modularconnectivities$kTotal)
modularconnectivities$kWithin<-ifelse(modularconnectivities$kWithin<0,0,modularconnectivities$kWithin)
modularconnectivities$kOut<-ifelse(modularconnectivities$kOut<0,0,modularconnectivities$kOut)

# Save node-wise connectivities
write.csv(modularconnectivities, file=file.path(newdir, "nodewise_modularconnectivities.csv"))

# calculate network size by dividing total genes by number of permuatations
networksize<-nrow(modularconnectivities)/length(levels(droplevels(modularconnectivities$X1)))
# add information on permutation
modularconnectivities<-merge(modularconnectivities,combinations,by.x="X1", by.y="row.names")

## calculate modulewise densities for every gene*treatment
#Calculate observed connectivities
modulardensities_2<-ddply(modularconnectivities, .(clusterID, stage, sexbias, X1), summarise, 
                          kMain=sum(kTotal, na.rm = T),
                          kWithin=sum(kWithin, na.rm = T),
                          kOut=sum(kOut, na.rm = T),
                          clusterSize=length(unique(transcriptID)))
# Calculate maximum possible connectivities
modulardensities_2$MaxkWithin<-maxConn(modulardensities_2$clusterSize)
modulardensities_2$MaxkOut<-(modulardensities_2$clusterSize*(networksize-modulardensities_2$clusterSize))/2
modulardensities_2$MaxkMain<-maxConn(networksize)
# Divide observed by maximum possible connectivities
modulardensities_2$dWithin<-modulardensities_2$kWithin/modulardensities_2$MaxkWithin
modulardensities_2$dOut<-modulardensities_2$kOut/modulardensities_2$MaxkOut
modulardensities_2$dMain<-modulardensities_2$kMain/modulardensities_2$MaxkMain

# Save as csv
write.csv(modulardensities_2, file.path(newdir, "sex_modulewise_modulardensities.csv"))

# plot changes in modularity parameters
# reorder stages
modulardensities_2$stage<-factor(modulardensities_2$stage, c("emb10","emb18", "lar51", "pupyel", "adult"))

pdf(file = file.path(graphdir, "modularity_plots_clusterdensities.pdf"), paper = "a4")

ggplot(modulardensities_2, aes(x=sexbias, y=dWithin, col=stage))+geom_point()+theme_bw()+geom_smooth(method="lm")+facet_wrap(~clusterID)+scale_color_brewer(type = "seq", palette = 4)+ggtitle("Integration vs Number of Female Samples Removed")
ggplot(modulardensities_2, aes(x=sexbias, y=dOut, col=stage))+geom_point()+theme_bw()+geom_smooth(method="lm")+facet_wrap(~clusterID)+scale_color_brewer(type = "seq", palette = 4)+ggtitle("Pleiotropy vs Number of Female Samples Removed")
ggplot(modulardensities_2, aes(x=dWithin, y=dOut, col=stage))+geom_point()+theme_bw()+geom_smooth(method="lm")+facet_wrap(~clusterID)+scale_color_brewer(type = "seq", palette = 4)+scale_x_log10()+scale_y_log10()+ggtitle("Integration vs Pleiotropy")

dev.off()

###### Outdated code snippets

# nasdevclusters<-droplevels(nasdevclusters[1:500,])  # limit to first transcripts for testing
# nasdevclusters<-nasdevclusters[which(nasdevclusters$clusterID%in%(levels(nasdevclusters$clusterID)[1:9])),] # limit to first 3 clusters only
# nasdevclusters<-nasdevclusters[grep("grey|antiquewhite", nasdevclusters$clusterID),] # limit to grey and red clusters
# nasdevclusters<-nasdevclusters[nasdevclusters$clusterID%in%problemclusters,] # limit to clusters with density>1 (solved, was using wrong maximum density)
# nasdevclusters<-nasdevclusters[nasdevclusters$clusterID%in%c("aquamarine3","snow4"),] # Filter two clusters containing some permutations with -1 connectivities (the -1 internal connections happen if all values = 0)
# nasdevclusters<-nasdevclusters[nasdevclusters$clusterID%in%c("darkorchid2","khaki","rosybrown3","paleturquoise4"),]# filter clusters containing non -1 negative connectivities, they are all near zero (E-15), likely rounding errors
# filter only transcripts present in the final network

# # Remove genes with sex-sepcifc expresion (they cannot have sex-specific interaction)
# sexspec<-read.csv(file = "./Output/sex_specific_nodes/sexspec.csv")[,-1]
# sexaspec<-sexspec[which(sexspec$spec=="Aspecific"),"eigenexonID"]
# nasdevgeneexon<-nasdevgeneexon[which(nasdevgeneexon$eigenexonID%in%sexaspec),]

# # TEMPORARY GENERATION OF RANDOM CLUSTERS, for testing only
# clusters<-data.frame(rpois(n=nrow(nasdevgeneexon),lambda = 4),row.names=row.names(nasdevgeneexon))
# colnames(clusters)<-"cluster"
# clusters$cluster[runif(n=nrow(clusters)/10, min=0, max=nrow(clusters))]<-NA
# clusters$geneID<-row.names(clusters)
