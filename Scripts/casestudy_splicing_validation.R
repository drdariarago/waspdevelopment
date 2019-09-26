## Validate splicing by comparing with known genes
## Version using gene list plus apply loop

# Initialize script
rm(list=ls())
library(ggplot2)
library(plyr)
library(reshape2)
library(stringr)

## Import eigenexon Evalues, genome annotation, eigenexon assignments and raw score files
eigenexon_evalues <- read.csv("./Output/exon_to_transcript_clustering_rankingonly/eigenexon_evalues.csv")
OGS2 <- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")
eigenexons_assignments <- read.csv("./Output/exon_to_transcript_clustering_rankingonly/eigenexons_assignments.csv")
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")


## Create gene list

geneIDs<-as.character(
  OGS2[grep("Doublesex", OGS2$Name, value = F)[1],"geneID"])

## Create results list

splicing_data<-list(
  gID=geneIDs,
  eigenexons=NA
  )

## Fill results list
lapply(X = splicing_data, FUN = function(x){
  x$eigenexons<-grep(x["gID"], eigenexon_evalues$eigenexonID, value=T)
})