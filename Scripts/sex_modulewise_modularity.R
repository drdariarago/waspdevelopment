# Calculate network parameters on a module-by-module basis
# Load packages
date()
rm(list=ls())
library(plyr)
handshake<-function(x){(x*(x+1))/2}

# initialize output path
newdir<-file.path(getwd(), "Output/sex_modulewise_modularity")
dir.create(newdir)

# load network data
modularconnectivities <- read.csv("./Output/sex_modularity_biweight/modularconnectivities.csv")
# calculate network size by dividing total genes by number of permuatations
networksize<-nrow(modularconnectivities)/length(levels(droplevels(modularconnectivities$X1)))
# remove extra columns
modularconnectivities<-modularconnectivities[,setdiff(names(modularconnectivities),c("X","V2","V2","V3","kDiff","dWithin","dOut","dMain"))]

# calculate modulewise densities for every gene*treatment
modulardensities_2<-ddply(modularconnectivities, .(clusterID, stage, sexbias, X1), summarise, dMain=sum(kTotal), dWithin=sum(kWithin), dOut=sum(kOut), clusterSize=length(X1))

modulardensities_2$dWithin<-modulardensities_2$dWithin/handshake(modulardensities_2$clusterSize)
modulardensities_2$dOut<-modulardensities_2$dOut/((modulardensities_2$clusterSize*(networksize-modulardensities_2$clusterSize))/2)
modulardensities_2$dMain<-modulardensities_2$dMain/handshake(networksize)

# Save as csv
write.csv(modulardensities_2, file.path(newdir, "sex_modulewise_modulardensities.csv"))
