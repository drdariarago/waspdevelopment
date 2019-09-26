## Gonadbias vs devbias

# load data
clusterdata <- read.csv(file = "./Output/Results_compiler/clusterdata_full.csv")

# Check number of nodes with each devbias pattern
summary(clusterdata$gonadbias)
