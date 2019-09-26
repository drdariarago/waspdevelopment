## Create modules blockwise
rm(list=ls())
library(WGCNA)
allowWGCNAThreads()
# import dataset
nasdevgeneexon<-read.csv(file="./Output/exon_to_transcript_clustering/eigenexons.csv")
row.names(nasdevgeneexon)<-nasdevgeneexon$eigenexonID
nasdevgeneexon<-nasdevgeneexon[,which(colnames(nasdevgeneexon)%in%setdiff(names(nasdevgeneexon),"eigenexonID"))]
t_nasdevgeneexon<-t(nasdevgeneexon)

# Check parameters of the dataset for debugging
dim(t_nasdevgeneexon)
str(t_nasdevgeneexon)
t_nasdevgeneexon


# calculate modules
# NOTES
# maxBlockSize as suggested for 16 GB RAM
# power validated by the WGCNA_preanalyses script
# deepSplit set for maximum number of clusters

nasdevgeneexonWGCNA<-blockwiseModules(
  t_nasdevgeneexon,
  maxBlockSize=20000,
  corType="pearson",
  power=20,
  networkType="signed",
  TOMType="signed",
  deepSplit=4,
  minModuleSize=10,
  saveTOMs=T,
  numericLabels=T,
  verbose=5,
  )

save(nasdevgeneexonWGCNA, file = "./Output/nasdevgeneexonWGCNA.RData", ascii = F)
save(nasdevgeneexonWGCNA$colors, file="./Output/nasoniadevgeneexonWGCNA_clusters.txt", ascii=T)

colors<-as.data.frame(nasdevgeneexonWGCNA$colors, row.names = colnames(t_nasdevgeneexon))
colors2<-as.data.frame(nasdevgeneexonWGCNA$colors, row.names = colnames(t_nasdevgeneexon)[nasdevgeneexonWGCNA$goodGenes])

save(colors, file="./Output/WGCNA_colors.RData")
save(colors2, file="./Output/WGCNA_colors2.RData")