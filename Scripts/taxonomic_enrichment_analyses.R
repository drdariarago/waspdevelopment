# Analyze clusters for taxonomic patterns
# Initialize script
rm(list=ls())
library(ggplot2)
library(reshape2)
library(stringr)
newdir<-file.path(getwd(), "Output/taxonomic_enrichment_analyses")

# Load data
eigenexon_data <- read.csv("./Output/Results_compiler/eigenexon_data.csv")

# Plot general overview of taxonomic depth in the genome
ggplot(eigenexon_data, aes(x=Taxascan_max, fill=constitutive))+geom_histogram(position="dodge")+theme_bw()+scale_x_discrete(breaks=c(0:5),labels=c("Nasonia","Hymenoptera","Holometabola","Insecta","Arthropoda","Metazoa"))+theme_bw()

# Plot dWithin sexbias over developmental time of genes in each cluster, split by taxonomic depth and splicing
dWithin_Taxa<-eigenexon_data[,c("geneID","clusterID","eigenexonID","constitutive","Taxascan_max", grep("sexbias.*dWithin",names(eigenexon_data),value = T))]
dWithin_Taxa<-melt(dWithin_Taxa, id.vars = c("geneID","clusterID","eigenexonID","constitutive","Taxascan_max"))
dWithin_Taxa$stage<-str_extract(string = dWithin_Taxa$variable, pattern = "stage[[:alnum:]]*")
dWithin_Taxa$stage<-factor(dWithin_Taxa$stage, levels=c("stageemb10", "stageemb18", "stagelar51", "stagepupyel", "stageadult"))
dWithin_Taxa$value<-as.factor(sign(dWithin_Taxa$value))
# ggplot(dWithin_Taxa, aes(x=stage, y=sign(value), col=Taxascan_max))+geom_smooth()+theme_bw()

ggplot(dWithin_Taxa[which(is.na(dWithin_Taxa$value)==F),], aes(x=stage, fill=value))+geom_histogram(position="dodge")+theme_bw()
ggplot(dWithin_Taxa[which(is.na(dWithin_Taxa$value)==F),], aes(x=stage, fill=value))+geom_histogram(position="dodge")+theme_bw()+facet_grid(Taxascan_max~constitutive)

# Plot data for dOut sexbias over developmental time of genes in each cluster, split by taxonomic depth and splicing
dOut_Taxa<-eigenexon_data[,c("geneID","clusterID","eigenexonID","constitutive","Taxascan_max", grep("sexbias.*dOut",names(eigenexon_data),value = T))]
dOut_Taxa<-melt(dOut_Taxa, id.vars = c("geneID","clusterID","eigenexonID","constitutive","Taxascan_max"))
dOut_Taxa$stage<-str_extract(string = dOut_Taxa$variable, pattern = "stage[[:alnum:]]*")
dOut_Taxa$stage<-factor(dOut_Taxa$stage, levels=c("stageemb10", "stageemb18", "stagelar51", "stagepupyel", "stageadult"))
dOut_Taxa$value<-as.factor(sign(dOut_Taxa$value))

ggplot(dOut_Taxa[which(is.na(dOut_Taxa$value)==F),], aes(x=stage, fill=value))+geom_histogram(position="dodge")+theme_bw()
ggplot(dOut[which(is.na(dOut_Taxa$value)==F),], aes(x=stage, fill=value))+geom_histogram(position="dodge")+theme_bw()+facet_grid(Taxascan_max~constitutive)
