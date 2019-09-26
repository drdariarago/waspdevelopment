# Calculate network parameters on a module-by-module basis
# Load packages
date()
rm(list=ls())
library(plyr)
handshake<-function(x){(x*(x+1))/2}

# initialize output path
newdir<-file.path(getwd(), "Output/stage_modulewise_modularity")
dir.create(newdir)

# load network data
modularconnectivities <- read.csv("./Output/stage_modularity_biweight_pseudopaired/modularconnectivities.csv")
# calculate network size by dividing total genes by number of permuatations
networksize<-nrow(modularconnectivities)/length(levels(as.factor(modularconnectivities$X1)))
# remove extra columns
modularconnectivities<-modularconnectivities[,setdiff(names(modularconnectivities),c("X","V2","V2","V3","kDiff","dWithin","dOut","dMain","x1","x2","x3","x4","x5","x6"))]

# calculate modulewise densities for every gene*treatment
modulardensities_2<-ddply(modularconnectivities, .(clusterID, X1), summarise, dMain=sum(kTotal), dWithin=sum(kWithin), dOut=sum(kOut), clusterSize=length(X1), emb10=emb10[1], emb18=emb18[1], lar51=lar51[1], pupyel=pupyel[1], adult=adult[1])

modulardensities_2$dWithin<-modulardensities_2$dWithin/handshake(modulardensities_2$clusterSize)
modulardensities_2$dOut<-modulardensities_2$dOut/((modulardensities_2$clusterSize*(networksize-modulardensities_2$clusterSize))/2)
modulardensities_2$dMain<-modulardensities_2$dMain/handshake(networksize)

# Save as csv
write.csv(modulardensities_2, file.path(newdir, "stage_modulewise_modulardensities.csv"))