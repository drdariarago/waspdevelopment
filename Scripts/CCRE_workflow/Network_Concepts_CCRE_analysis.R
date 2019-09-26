## Analyze network concepts generated via CCRE within cluster
## Initialize script
date()
rm(list=ls())
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(WGCNA)
allowWGCNAThreads()
## initialize output path
newdir<-file.path(getwd(), "Output/NetworkConcepts_CCRE_analysis")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/NetworkConcepts_CCRE_analysis")
dir.create(graphdir)

## Load dataset
load("./Output/NetworkConcepts_CCRE/cluster_emp_networkconcepts.RData")
nodeparams <- read.csv("./Output/NetworkConcepts_CCRE/CCRE_nodeparams.csv")[,-1]
nodeparams$grey <- as.factor(ifelse(nodeparams$clusterID=="grey","grey","cluster"))
clusterparams <- read.csv("./Output/NetworkConcepts_CCRE/CCRE_clusterparams.csv")[,-1]

# Plotting basic cluster metrics
ggplot(clusterparams, aes(x=Density, y=Centralization, size = Heterogeneity, col = clusterID=="grey"))+
  geom_point()+geom_smooth()+scale_x_log10()+scale_y_log10()+theme_bw()
ggplot(clusterparams, aes(x=Centralization, y=Heterogeneity, size = Density, col = clusterID=="grey"))+
  geom_point()+geom_smooth()+theme_bw()
pairs(log10(clusterparams[,2:4]))

# Denser clusters are more centralized, more centralized clusters have higher heterogeneity

# Plotting basic node metrics
pairs(log10(nodeparams[,4:6]))
ggplot(nodeparams, aes(x=ScaledConnectivity))+geom_density()+theme_bw()+scale_x_log10()
ggplot(nodeparams, aes(x=ScaledConnectivity, y=ClusterCoef))+geom_point()+geom_smooth()
ggplot(nodeparams, aes(x=MAR, y=ClusterCoef, col=clusterID))+geom_point()+geom_smooth(method="lm")
ggplot(nodeparams, aes(x=ScaledConnectivity, y=ClusterCoef, col=log10(MAR), alpha=0.2))+
  geom_point()+geom_smooth(method="lm")+facet_wrap(~grey)+scale_x_log10()+scale_y_log10()+scale_color_distiller(palette = 4)
ggplot(nodeparams, aes(x=rank(ScaledConnectivity), y=rank(ClusterCoef)))+geom_point()+facet_wrap(~grey)+geom_smooth()
