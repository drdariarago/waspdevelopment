### Perform PCA as per DIDE_GLMs, but adding excess methylation proportion per cluster as factor

library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(MuMIn)
library(corrgram)
options(na.action = na.fail)
# Clean workspace
rm(list=ls())
z <- function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}
# Save text size (article KOMA)
textWidthHeightmm <- c(152.73, 216.00)
# Load multiplot function
source(file = "./Scripts/multiplot.R")
# Load colorblind friendly palette with grey
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Load colorblind palette for unbiased, female, male, both
cbPalette4 <- c("#999999", "#D55E00", "#0072B2", "#CC79A7")
# Create output path
newdir <- file.path(getwd(), "Output/PCA_clusters")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/PCA_clusters")
dir.create(graphdir)

## Load dataset
clusterdata <- read.csv(file = "./Output/Results_compiler_CCRE/clusterdata_full.csv")[,-1]
## Annotate DIDE factor
clusterdata$DIDE <- as.factor((1-is.na(clusterdata$DiffIntegrated))+2*grepl(pattern = "[m,f]", clusterdata$cluster_devsexbias))
clusterdata$DIDE <- factor(x = clusterdata$DIDE, labels = c("Unbiased","DI","DE","DIDE"))
clusterdata$DE <- grepl(pattern = "[m,f]", x = clusterdata$cluster_devsexbias)
clusterdata$DI <- clusterdata$DIDE=="DIDE"
clusterdata <- clusterdata[,c("clusterID", "cluster_devsexbias", "DiffIntegrated", "DIDE", "Density", "Centralization", "Heterogeneity", "nNodes", "excessSpl", "excessDup", "excessMet", "medianClusterCoef", "diameter", "DI", "DE")]

## Optional: remove grey cluster
clusterdata <- droplevels(clusterdata[-which(clusterdata$clusterID=="grey"),])

# Plot raw correlations
library(corrgram)
cbPalette4
cbPalette3 <- c("#999999","#D55E00","")

pdf(file = file.path(graphdir, "NetworkCor.pdf"), width = textWidthHeightmm[1]/1*.039, height = textWidthHeightmm[1]/1*.039)
names(clusterdata) <- c("clusterID", "cluster_devsexbias", "DiffIntegrated", "DIDE", "Density", "Centralization", "Heterogeneity", "Number\nof Nodes", "Excess\nSplicing\nProportion", "Excess\nDuplicate\nProportion", "Excess\nMethylated\nProportion", "Median\nClustering\nCoefficient","Diameter","DI","DE")
corrgram(clusterdata[,3:13], order=TRUE, lower.panel=panel.shade,upper.panel=panel.conf, text.panel=panel.txt, col.regions = colorRampPalette(cbPalette[c(3,1,5)])) # set labels parameter instead of changing factor names if necessary to run rest of the script
dev.off()

# Rescale numeric variables and reduce to relevant factors
clusterdata[,5:13] <- scale(x = clusterdata[,5:13], center = T, scale = T)

# Check PCs in our dataset
NetPCA <- prcomp(x = clusterdata[,5:13])
NetPCA
summary(NetPCA)
plot(NetPCA)
write.csv(x = NetPCA$rotation, file = file.path(newdir, "PCAScoreValues.csv"))

# Apply lm using all PCs as predictors
PCclusterdata <- cbind(clusterdata[,c(1:4,14:15)],NetPCA$x)
# DI
PCglmDI <- glm(formula = DI ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9, family = binomial, data = PCclusterdata)
summary(model.avg(dredge(PCglmDI)))
# -PC 1 (.93), -PC 5 (.85), PC 3 (.66)

# DE
PCglmDE <- glm(formula = DE ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9, family = binomial, data = PCclusterdata)
summary(model.avg(dredge(PCglmDE)))
# PC 7 (.90), -PC 5 (.85), PC 6 (.76)

# DE vs DI
PCglmDEvDI <- glm(formula = DI ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9, family = binomial, data = droplevels(PCclusterdata[which(PCclusterdata$DIDE!="Unbiased"),]))
summary(model.avg(dredge(PCglmDEvDI)))
# -PC1 (.96), -PC5 (.62)

# Interesting PCs 1,3,5,6,7
NetPCA2 <- data.frame(PCclusterdata[,c(1:6)],NetPCA$x)
# apply(X = NetPCA$rotation, MARGIN = 2, function(x){x-min(x)})
apply(X = NetPCA$rotation, MARGIN = 2, function(x){x/max(abs(x))})
ggplot(data = melt(apply(X = NetPCA$rotation, MARGIN = 2, abs)), mapping = aes(y=value, x=Var2, label=Var1))+geom_text()
# ggplot(data = melt(apply(X = NetPCA$rotation, MARGIN = 2, abs)), mapping = aes(y=value, x=Var2, col=Var1, shape=Var1))+geom_point()+scale_color_brewer(type = "qual")

### PC interpretation
## PC1: Density-dependent effects (Main cofactors ClusteringCoef and -NumberofNodes)
# PC2: Heterogeneity-dependent effects (Main cofactors Duplicates and -Methylated, positiveli correlated with centralization)
## PC3: Splicing effects (main cofactors Duplicates and methylated)
# PC4: Sparse-ness? (heterogeneity 2) (main cofactors -nNodes, -Diameter and -Density, negatively correlated with centralization)
## PC5: Splicing 2 (main cofactors pDuplicates and -Heterogeneity)
## PC6: Methylated (main cofactors duplicates)
## PC7: Hierarchy (Centralization, main cofactor -ClusteringCoefficient)
# PC8: Size-effects (nNodes, main cofactors -Diameter and ClusteringCoef)

library(ggbiplot)
ggbiplot(NetPCA, choices = c(1,2))

## Separate DE from rest (5 and 6 negative, 7 positive)
ggbiplot(NetPCA, groups = PCclusterdata$DE, choices = c(5,6), ellipse = T) + theme_bw() + 
  scale_color_brewer(type = "qual", palette = 3)
ggbiplot(NetPCA, groups = PCclusterdata$DE, choices = c(5,7), ellipse = T) + theme_bw() + 
  scale_color_brewer(type = "qual", palette = 3)
ggbiplot(NetPCA, groups = PCclusterdata$DE, choices = c(6,7), ellipse = T) + theme_bw() + 
  scale_color_brewer(type = "qual", palette = 3)

## Separate DI from rest (1,5, negative; 3, positive)
ggbiplot(NetPCA, groups = PCclusterdata$DI, choices = c(1,5), ellipse = T) + theme_bw() + 
  scale_color_brewer(type = "qual", palette = 3)
ggbiplot(NetPCA, groups = PCclusterdata$DI, choices = c(1,3), ellipse = T) + theme_bw() + 
  scale_color_brewer(type = "qual", palette = 3)
ggbiplot(NetPCA, groups = PCclusterdata$DI, choices = c(3,5), ellipse = T) + theme_bw() + 
  scale_color_brewer(type = "qual", palette = 3)

## Separate DE from DI (1 and 5, negative)
ggbiplot(NetPCA, groups = PCclusterdata$DIDE, choices = c(1,5), ellipse = T) + theme_bw() + 
  scale_color_brewer(type = "qual", palette = 3)


## -PC5 is the best predictor in DE vs others and significant in DI vs others but not DE vs DI
biplot(NetPCA, choices = c(1,5))
# Indicates they are more centralized and heterogeneous than expected given the excess of splicing and duplicated nodes
# Also potentially separates DI from DE (47% RI), hinting that DI might be even more hierarchical 

## -PC1 is characteristic of DI, indicates that they are more clustered and denser than expected per number of nodes
## Only separator between DE and DI

## -PC6 and -PC8 are characteristic of DE, 
## -PC6 indicates higher proportion of duplicates and methylated genes
## -PC8 indicates Higher centralization and lower median clustering coefficient (hierarchical index)

## PC3 is shared among DE and DI but differentiates DI from other clusters
biplot(NetPCA, choices = c(1,3))
# PC3 separates clusters enriched in duplicates from clusters rich in methylated genes and splicing. 
# DI clusters are more often duplicate rich and splicing/methylation poor

NetPCA2 <- data.frame(PCclusterdata[,c(1:6)],NetPCA$x)

apply(X = NetPCA$rotation, MARGIN = 2, function(x){x-min(x)})
apply(X = NetPCA$rotation, MARGIN = 2, function(x){x/max(abs(x))})

ggplot(data = NetPCA2, mapping = aes(x = DIDE, y = PC1)) + geom_boxplot(notch = T)
ggplot(data = NetPCA2, mapping = aes(x = DIDE, y = PC3)) + geom_boxplot(notch = T)

ggplot(data = NetPCA2, mapping = aes(x = DIDE, y = PC5)) + geom_boxplot(notch = T)
ggplot(data = NetPCA2, mapping = aes(x = DIDE, y = PC6)) + geom_boxplot(notch = T)
ggplot(data = NetPCA2, mapping = aes(x = DIDE, y = PC8)) + geom_boxplot(notch = T)


## Report for Jack

pdf(file = file.path(graphdir, "Met_report.pdf"))
corrgram(clusterdata[,3:13], order=TRUE, lower.panel=panel.shade,upper.panel=panel.conf, text.panel=panel.txt, col.regions = colorRampPalette(cbPalette[c(3,1,5)]))
ggbiplot(NetPCA, groups = PCclusterdata$DE, choices = c(3,6), ellipse = T) + theme_bw() + 
  scale_color_brewer(type = "qual", palette = 3, name = "Differentially Sex-Biased?")
dev.off()
