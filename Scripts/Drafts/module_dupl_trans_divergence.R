#### Check for gene duplication/splicing hypothesis
# Do genes within the same OG belong in the same clusters?
# Do transcripts within the same gene belong in the same clusters?
# Also correct for sequence divergence at hymenoptera

# load packages
library(plyr)
library(stringr)
library(lme4)
rm(list=ls())
date()
newdir<-file.path(getwd(), "Output/module_dupl_trans_divergence")
dir.create(newdir)

# load data
OGS2 <- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")
clusterassignments <- read.csv("./Output/WGCNA_clustering_biweight/clusterassignments_simple.csv")
# Filter only columns of interest
OGS2m<-OGS2[,c("geneID","ODB6_OG_ID","OG_Copynumber","Nresiduals","adult_female_meth_status")]
# Exclude methylation status (~1435 missing entries)
OGS2<-OGS2m[,-grep("adult_female_meth_status", names(OGS2m))]
# Remove data with missing entries or duplicates
OGS2<-na.exclude(OGS2)
OGS2<-unique(OGS2)
# and for methylation dataset
OGS2m<-na.exclude(OGS2m)
OGS2m<-unique(OGS2m)
# Merge with cluster assignments
clusterassignments$geneID<-as.factor(str_extract(pattern = "^[[:alnum:]]*",string = clusterassignments$transcriptID))
clusterdata<-merge(clusterassignments, OGS2, by="geneID")
clusterdata<-clusterdata[order(clusterdata$transcriptID),]
# and for methylation dataset
clusterdatam<-merge(clusterassignments, OGS2m)
clusterdatam<-clusterdatam[order(clusterdatam$transcriptID),]


# # remove OGs present in single copy in the final dataset
# duplOG<-unique(clusterdata[,c("ODB6_OG_ID","geneID")])
# duplOG<-names(table(duplOG$ODB6_OG_ID)[which(table(duplOG$ODB6_OG_ID)>1)])
# clusterdata<-clusterdata[which(clusterdata$ODB6_OG_ID%in%duplOG),]

# summarize data
# How many OGs are present in more than 1 cluster?
nOG_Clusters<-ddply(.data = clusterdata, .variables = .(ODB6_OG_ID), summarize, nClusters=length(levels(droplevels(clusterID))), .progress = "text")
table(nOG_Clusters$nClusters>1)
prop.table(table(nOG_Clusters$nClusters>1))
# How many genes are present in more than 1 cluste?r
nTr_Clusters<-ddply(.data = clusterdata, .variables = .(geneID), summarize, nClusters=length(levels(droplevels(clusterID))), .progress = "text")
table(nTr_Clusters$nClusters>1)
prop.table(table(nTr_Clusters$nClusters>1))

# nested: og/gene/transcipt belong to the same module (crate series of models, then compare via AIC)

glmerDOGT<-glmer(formula = clusterID~Nresiduals+(1|ODB6_OG_ID/geneID), data = clusterdata, family = binomial)
glmerDOG<-glmer(formula = clusterID~Nresiduals+(1|ODB6_OG_ID), data = clusterdata, family = binomial)
glmerDT<-glmer(formula = clusterID~Nresiduals+(1|geneID), data = clusterdata, family = binomial)
glmerOGT<-glmer(formula = clusterID~1+(1|ODB6_OG_ID/geneID), data = clusterdata, family = binomial)
glmerOG<-glmer(formula = clusterID~1+(1|ODB6_OG_ID), data = clusterdata, family = binomial)
glmerT<-glmer(formula = clusterID~1+(1|geneID), data = clusterdata, family = binomial)
glmD<-glm(formula = clusterID~Nresiduals, data = clusterdata, family = binomial)
glmnull<-glm(formula = clusterID~1, data = clusterdata, family = binomial)

AIC(glmerDOGT,glmerDOG,glmerDT,glmerOGT,glmerOG,glmerT,glmD,glmnull)

# modulewise: check within every module whether OG, gene and distance predict transcript assignment to that module (binomial, zero inflated?)
# annotate presence or absence from every module as a different variable (0/1)
clusterpresence<-sapply(levels(clusterdata$clusterID), FUN = function(x){clusterdata$clusterID%in%x})
clusterdata<-cbind(clusterpresence, clusterdata)

# ## solving problem as a nonrandom distribution of 2 entry matrix (OG by cluster)
# glm(formula = clusterID~ODB6_OG_ID, data=clusterdata[1:500,], family=binomial)

overlapTable # could be useful for cross-tabulation of ODB/genes vs modules
clusterdata<-clusterdata[1:100,]
a<-coClustering.permutationTest(clusters.ref = clusterdata$clusterID, clusters.test = clusterdata$ODB6_OG_ID, nPermutations = 100, verbose=3)

# Solve problem trough use of TOM matrix: is the overlap within OGs greater or smaller than the overlap between random groups of transcripts? same for genes

# load(file="./Output/WGCNA_clustering_biweight/dTOM")

# create small dTOM 
library(WGCNA)
allowWGCNAThreads()
nasdevgeneexon<-read.csv(file="./Output/exon_to_transcript_clustering_rankingonly/eigenexon_evalues.csv")
row.names(nasdevgeneexon)<-nasdevgeneexon$eigenexonID
nasdevgeneexon<-nasdevgeneexon[1:1000,] # reduce size for testing
nasdevgeneexon<-nasdevgeneexon[,which(colnames(nasdevgeneexon)%in%setdiff(names(nasdevgeneexon),"eigenexonID"))]
t_nasdevgeneexon<-transposeBigData(nasdevgeneexon, blocksize = 50000)
goodCols<-goodSamplesGenes(datExpr = t_nasdevgeneexon, minNSamples = 3)
lapply(goodCols, function(x) {table(x)})
lapply(goodCols, function(x) {prop.table(table(x))})
t_nasdevgeneexon<-t_nasdevgeneexon[,which(goodCols$goodGenes==T)]

threshold_biweight<-pickSoftThreshold(
  corFnc = bicor,
  t_nasdevgeneexon, 
  corOptions = list(quick=0.3, use = "pairwise.complete.obs"),
  powerVector = c(1:10,seq(from = 12,to = 21,by = 2),seq(from = 25, to = 30, by = 5)),
  networkType="signed",
  moreNetworkConcepts=T,
  verbose=5)

adj<-adjacency(
  datExpr = t_nasdevgeneexon,
  power = threshold_biweight$powerEstimate, 
  corFnc = "bicor",
  corOptions = "quick=0, use = 'pairwise.complete.obs'",
  type = "signed")

dTOM<-TOMdist(adj, TOMType = "signed")
