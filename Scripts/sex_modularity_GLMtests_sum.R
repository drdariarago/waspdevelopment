## test appropriate distribution for glm based analyses
rm(list=ls())
library(plyr)
# load connectivities
modularconnectivities <- read.csv("./Output/stage_modularity_biweight/modularconnectivities.csv")
# sum modular connectivities for each module
ddply(.data = modularconnectivities, .variables = .(clusterID, ))
