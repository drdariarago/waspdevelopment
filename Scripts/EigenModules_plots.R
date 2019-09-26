## Plot mean eigenmodule behavior over time
date()
rm(list=ls())
library(ggplot2)
library(reshape2)
library(flashClust)
library(WGCNA)
library(plyr)
library(stringr)
allowWGCNAThreads()
source(file = "./Scripts/multiplot.R")

# Load data
load(file = "./Output/WGCNA_clustering_biweight/WGCNA_clustering_biweight.RData")
clusterdata <- read.csv(file = "./Output/Results_compiler/clusterdata_full.csv")[,-1]

# initialize output path
newdir<-file.path(getwd(), "Output/EigenModules_plots")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/EigenModules_plots")
dir.create(graphdir)

EigenModuleE <- t(MElist2$newMEs)
colnames(EigenModuleE) <- row.names(t_nasdevgeneexon)

EigenModuleE<- as.data.frame(EigenModuleE)
row.names(EigenModuleE) <- str_extract(row.names(EigenModuleE), pattern = "[^[:upper:]]*$")
EigenModuleE$moduleID<- factor(row.names(EigenModuleE), row.names(EigenModuleE)[MElist2$dendro$order])

EigenModuleE<- melt(EigenModuleE)
EigenModuleE$Sex <- ifelse(grepl("female",EigenModuleE$variable),"Female","Male")
EigenModuleE$Stage <- factor(str_extract(string = EigenModuleE$variable, pattern = "^[^_]*"), levels = c("emb10","emb18","lar51","pupyel","adult"))
EigenModuleE$Replicate <- as.factor(str_extract(string = EigenModuleE$variable, pattern = ".$"))
EigenModuleE <- EigenModuleE[,c("moduleID","value","Sex","Stage","Replicate")]

# ggplot(data = EigenModuleE, aes(x=Stage, y=value, col=Sex))+geom_point()+facet_wrap( ~ moduleID)+theme_bw()
# ggsave(filename = file.path(graphdir,"MEplots.pdf"), width = 30, height = 30)

dWithin <- clusterdata[,c("clusterID","dWithinBias")]
names(dWithin) <- c("moduleID","dWithinBias")
EigenModuleE <- merge(EigenModuleE, dWithin, by = "moduleID")

ggplot(data = EigenModuleE, aes(x=Stage, y=value, pch=Sex, col=dWithinBias))+geom_point()+facet_wrap( ~ moduleID)+theme_bw()
ggsave(filename = file.path(graphdir,"MEplots_dWithin.pdf"), width = 30, height = 30)

ggplot(data = EigenModuleE, aes(x=Stage, y=value+1, pch=Sex, col=dWithinBias))+geom_point()+facet_wrap( ~ moduleID)+theme_bw()+scale_y_log10()
ggsave(filename = file.path(graphdir,"MEplots_dWithin_log10.pdf"), width = 30, height = 30)