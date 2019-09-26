## Use GLM to detect unique features of DE and DIDE clusters
# Fit 3 models: nonDE vs DE, nonDE vs DIDE, DE vs DI
# Load packages
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
newdir <- file.path(getwd(), "Output/DIDE_GLMs")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/DIDE_GLMs")
dir.create(graphdir)

## Load dataset
clusterdata <- read.csv(file = "./Output/Results_compiler_CCRE/clusterdata_full.csv")[,-1]
## Annotate DIDE factor
clusterdata$DIDE <- as.factor((1-is.na(clusterdata$DiffIntegrated))+2*grepl(pattern = "[m,f]", clusterdata$cluster_devsexbias))
clusterdata$DIDE <- factor(x = clusterdata$DIDE, labels = c("Unbiased","DI","DE","DIDE"))
clusterdata$DE <- grepl(pattern = "[m,f]", x = clusterdata$cluster_devsexbias)
clusterdata$DI <- clusterdata$DIDE=="DIDE"
clusterdata <- clusterdata[,c("clusterID", "cluster_devsexbias", "DiffIntegrated", "DIDE", "Density", "Centralization", "Heterogeneity", "nNodes", "excessSpl", "excessDup", "medianClusterCoef", "diameter", "DI", "DE")]

# Do we see a small-world effect?
lm1 <- lm(formula = diameter ~ nNodes, data = clusterdata)
lm2 <- lm(formula = diameter ~ log(nNodes), data = clusterdata)
lm3 <- lm(formula = log(diameter) ~ nNodes, data = clusterdata)
lm4 <- lm(formula = log(diameter) ~ log(nNodes), data = clusterdata)
lm5 <- lm(formula = diameter ~ log(log(nNodes)), data = clusterdata)
lm6 <- lm(formula = log(diameter) ~ log(log(nNodes)), data = clusterdata)
AICc(lm1, lm2, lm3, lm4, lm5, lm6) #  best fit is log-log linear, coef>1.6, no small-world effect
summary(lm4)
ggplot(data = clusterdata, mapping = aes(x = log(nNodes), y = log(diameter), col = DIDE)) + geom_point() + geom_smooth(method = "lm", aes(group = grepl(".", cluster_devsexbias)))
 
# lm1 <- lm(formula = diameter ~ nTranscripts, data = clusterdata)
# lm2 <- lm(formula = diameter ~ log(nTranscripts), data = clusterdata)
# lm3 <- lm(formula = log(diameter) ~ nTranscripts, data = clusterdata)
# lm4 <- lm(formula = log(diameter) ~ log(nTranscripts), data = clusterdata)
# lm5 <- lm(formula = diameter ~ log(log(nTranscripts)), data = clusterdata)
# lm6 <- lm(formula = log(diameter) ~ log(log(nTranscripts)), data = clusterdata)
# AICc(lm1, lm2, lm3, lm4, lm5, lm6) #  some evidence for small-world when using nTranscripts
# ggplot(data = clusterdata, mapping = aes(x = log(nTranscripts), y = log(log(diameter)), col = DIDE)) + geom_point() + geom_smooth(method = "lm", aes(group = grepl(".", cluster_devsexbias)))


# Plot raw correlations
library(corrgram)
cbPalette4
cbPalette3 <- c("#999999","#D55E00","")

## Remove cluster GREY
clusterdata <- clusterdata[-which(clusterdata$clusterID=="grey"),]

# pdf(file = file.path(graphdir, "NetworkCor.pdf"), width = textWidthHeightmm[1]/1*.039, height = textWidthHeightmm[1]/1*.039)
# names(clusterdata) <- c("clusterID", "cluster_devsexbias", "DiffIntegrated", "DIDE", "Density", "Centralization", "Heterogeneity", "Number\nof Nodes", "Excess\nSplicing\nProportion", "Excess\nDuplicate\nProportion", "Median\nClustering\nCoefficient","Diameter","DI","DE")
# corrgram(clusterdata[,3:12], order=TRUE, lower.panel=panel.shade,upper.panel=panel.conf, text.panel=panel.txt, col.regions = colorRampPalette(cbPalette[c(3,1,5)])) # set labels parameter instead of changing factor names if necessary to run rest of the script
# dev.off()

# Rescale numeric variables and reduce to relevant factors
clusterdata[,5:12] <- scale(x = clusterdata[,5:12], center = T, scale = T)

# Check PCs in our dataset
NetPCA <- prcomp(x = clusterdata[,5:12])
NetPCA
summary(NetPCA)
plot(NetPCA)
write.csv(x = NetPCA$rotation, file = file.path(newdir, "PCAScoreValues.csv"))

par(mfrow = c(7,4))
combos <- data.frame(t(combn(x = 1:8,2, simplify = T)))
apply(combos, 1, function(x){
  biplot(NetPCA, choices = c(x))
})
par(mfrow = c(1,1))
# biplot(NetPCA, choices = c(1,2))
# biplot(NetPCA, choices = c(1,3))
# biplot(NetPCA, choices = c(1,4))
# biplot(NetPCA, choices = c(2,3))
# biplot(NetPCA, choices = c(2,4))
biplot(NetPCA, choices = c(3,5))

# First PC explains ~51% total variance, inv density/centralization/clustercoef vs nNodes and Diameter (size vs density)
# Second PC explains 17%, due to heterogeneity, centralization and diameter (inv scale-free ness?)
# Third PC explains 13%, duplicates vs splicing
# Fourth PC explains 10%, Heterogeneity/Splicing vs density/diameter/nNodes (random-ness?/sparse-ness?)
# Fifth PC explains 7%, duplication/splicing vs centralization/heterogeneity
# Sixth PC explains 1.1%, centralization vs clusterCoef (hierarchy)
# Seventh PC, >1%, nNodes/clusterCoef vs Diameter/density (small world)
# Eight PC, >0.1%, Density

# pairs(PCclusterdata[,5:12], col = PCclusterdata$DIDE)

# Apply lm using all PCs as predictors
PCclusterdata <- cbind(clusterdata[,c(1:4,13:14)],NetPCA$x)
PCclusterdata$DIDE <- factor(x = PCclusterdata$DIDE, levels = c("Unbiased", "DE", "DI", "DIDE"))
# DI
PCglmDI <- glm(formula = DI ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, family = binomial, data = PCclusterdata)
summary(model.avg(dredge(PCglmDI)))
# -PC 1 (.92), PC 3 (.71), -PC 5 (.72)

# DE
PCglmDE <- glm(formula = DE ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, family = binomial, data = PCclusterdata)
summary(model.avg(dredge(PCglmDE)))
# PC 6 (.86), -PC 4 (.56), -PC 5 (.53)

# DE vs DI
PCglmDEvDI <- glm(formula = DI ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, family = binomial, data = droplevels(PCclusterdata[which(PCclusterdata$DIDE!="Unbiased"),]))
summary(model.avg(dredge(PCglmDEvDI)))
# -PC 1 (.95), -PC5 (.54)

# # DIDE both at the same time
# PCglmDIDE <- glm(DIDE=="DIDE" ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, family = binomial, data = PCclusterdata)
# summary(model.avg(dredge(PCglmDIDE)))
# # -PC 1 (.9), PC 3 (.8), PC 5 (.6), PC 7 (.6)
# 
# # DIDE, at least one
# PCglmDIDE2 <- glm(DIDE!="Unbiased" ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, family = binomial, data = PCclusterdata)
# summary(model.avg(dredge(PCglmDIDE2)))
# # PC 7 (.7), PC 3 (.6), -PC 2 (.5), -PC6 (.5)


## Result summary (DE and DI)
## -PC5 moderate support in both (.55 and .65)
## -PC2 low support in both (.54, .54)
## PC7 hi support in DE, low in DI (.77 and .3)
## -PC1 hi support DI, no support DE (.92, .33)
## PC3 moderate support in DI, low in DE (.71, .49)

# Ordered by PC
# -1 density, important for DI
# -2 heterogeneity+diameter, low importance in DE and DI 
# 3 duplicates vs splicing, moderate support for DI, low in DE
# -5 centralization/heterogeneity vs duplicates/diameter/splicing, moderate in both
# 7 centralization vs clustercoef, important in DE, low in DI 

# Plot PC -1 and 3 (DI)
ggplot(data = PCclusterdata, mapping = aes(x = -PC1, y = PC3, col = DI, shape = DIDE)) + geom_point() + theme_bw() + geom_rug(alpha="0.3") + geom_density2d()
# Plot PC 6 and 4 (DE)
ggplot(data = PCclusterdata, mapping = aes(x = PC4, y = PC6, col = DE, shape = DIDE)) + geom_point() + theme_bw() + geom_rug(alpha="0.3") + geom_density_2d()

# As boxplots
PC1plot <- ggplot(data = PCclusterdata, mapping = aes(x = DIDE, y = PC1, col = DIDE)) + geom_boxplot(notch = T, varwidth = T) + geom_point(position = position_jitter(width = .2))
PC1plot
PC2plot <- ggplot(data = PCclusterdata, mapping = aes(x = DIDE, y = PC2, col = DIDE)) + geom_boxplot(notch = T, varwidth = T) + geom_point(position = position_jitter(width = .2))
PC2plot
PC3plot <- ggplot(data = PCclusterdata, mapping = aes(x = DIDE, y = PC3, col = DIDE)) + geom_boxplot(notch = T, varwidth = T) + geom_point(position = position_jitter(width = .2))
PC3plot
PC5plot <- ggplot(data = PCclusterdata, mapping = aes(x = DIDE, y = PC5, col = DIDE)) + geom_boxplot(notch = T, varwidth = T) + geom_point(position = position_jitter(width = .2))
PC5plot
PC7plot <- ggplot(data = PCclusterdata, mapping = aes(x = DIDE, y = PC7, col = DIDE)) + geom_boxplot(notch = T, varwidth = T) + geom_point(position = position_jitter(width = .2))
PC7plot

multiplot(PC1plot, PC2plot, PC3plot, PC5plot, PC7plot, cols = 3)

# In only 1 ggplot
mPCclusterdata <- PCclusterdata[,c(1,4,7:14)]
mPCclusterdata[,3:10] <- scale(x = mPCclusterdata[,3:10], center = F, scale = T)
# mPCclusterdata$PC1 <- -mPCclusterdata$PC1 # reverse PC1 to make it proportional to density/centralization/clustercoef and inv to nNodes
# mPCclusterdata$PC2 <- -mPCclusterdata$PC2 # reverse PC2 to make it proportional to smallness (less diameter, less heterogeneity)
# mPCclusterdata$PC5 <- -mPCclusterdata$PC5 # reverse PC2 to make it proportional to centralization vs dupl/diam/splic
mPCclusterdata <- melt(data = mPCclusterdata, id.vars = c("clusterID","DIDE"))
mPCclusterdata$DIDE <- factor(x = mPCclusterdata$DIDE, labels = c("No Bias", "DE Only", "DC Only", "DC and DE"))
# mPCclusterdata$DIDE <- factor(x = mPCclusterdata$DIDE, levels = c("No Bias", "DC Only", "DE Only", "DC and DE"))
mPCclusterdata <- droplevels(mPCclusterdata[which(mPCclusterdata$variable%in%c("PC1","PC3","PC6")),])
mPCclusterdata$variable <- factor(x = mPCclusterdata$variable, labels = c("PC1 50% var\nSize versus Density","PC3 13% var\nDuplicates versus Splicing","PC6 >1% var\nHubness"))

cbPalette[2:3] <- cbPalette[c(3,2)]

ggplot(data = mPCclusterdata, mapping = aes(x = DIDE, y = value, col = DIDE)) +
  geom_boxplot(notch = F, varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.3), alpha = "0.3") + facet_wrap(~variable) +
  theme_bw() + ylab(label = "Scaled Loading Values") + xlab(label = "") + scale_y_continuous(limits = c(-3,3)) + scale_color_manual(values = cbPalette, name = "Sex-Bias Type") +
  theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom") + guides(col = guide_legend(ncol = 2))
# legend.position = c(0.8,0.25)) 
ggsave(filename = file.path(graphdir, "NetworkBoxplots.pdf"), device = "pdf", width = textWidthHeightmm[1], height = textWidthHeightmm[2]/2, units = "mm")

# CP PC 3 and PC 5
ggplot(data = PCclusterdata, mapping = aes(x = PC3, y = PC5, col =  grepl(pattern = "DE", DIDE))) + geom_point() + geom_rug() + geom_density2d() + geom_smooth(method = "lm") + scale_x_continuous(limits = c(-1.5,1.5)) + scale_y_continuous(limits = c(-1.5,1.5))
ggplot(data = PCclusterdata, mapping = aes(x = PC3, y = PC5, col =  grepl(pattern = "DIDE", DIDE))) + geom_point() + geom_rug() + geom_density2d() + geom_smooth(method = "lm") + scale_x_continuous(limits = c(-1.5,1.5)) + scale_y_continuous(limits = c(-1.5,1.5))

####### Do different strata have different connectivities?
transcriptdata <- read.csv(file = "./Output/Results_compiler_CCRE/transcriptdata_full.csv")[,-1]
# ggplot(data = transcriptdata, mapping = aes(x = strata, y = AbsConnectivity_wc)) + geom_boxplot(notch = T, varwidth = T) #YES!

transcriptdata <- merge(transcriptdata, clusterdata[,c("clusterID","DIDE","DI","DE", "nNodes")])
transcriptdata$grouper <- paste(transcriptdata$strata, transcriptdata$DIDE)
transcriptdata <- droplevels(transcriptdata[-which(is.na(transcriptdata$strata)),])
transcriptdata$strata <- factor(x = transcriptdata$strata, levels = c("Wasp", "Hymenoptera", "Insect", "Arthropod", "Metazoa"))
transcriptdata$strata <- factor(x = transcriptdata$strata, labels = c("Nasonia", "Hymenoptera", "Insecta", "Arthropoda", "Metazoa"))
# Filter only representative nodes
transcriptdata <- transcriptdata[which(transcriptdata$Representative_Node==T|is.na(transcriptdata$Representative_Node)),]

ggplot(data = transcriptdata, mapping = aes(x = strata, y = AbsConnectivity_wc)) + geom_boxplot(notch = T, varwidth = T) + facet_grid(~DIDE)
ggplot(data = transcriptdata, mapping = aes(x = strata, y = AbsBetweenness_wc)) + geom_boxplot(notch = T, varwidth = T) + facet_grid(~DIDE) + scale_y_log10()
ggplot(data = transcriptdata, mapping = aes(x = strata, y = ClusterCoef_wc)) + geom_boxplot(notch = T, varwidth = T) + facet_grid(~DIDE) + scale_y_log10()
ggplot(data = transcriptdata, mapping = aes(x = strata, y = MAR_wc)) + geom_boxplot(notch = T, varwidth = T) + facet_grid(~DIDE) + scale_y_log10()

ggplot(data = transcriptdata, mapping = aes(x = strata, y = AbsConnectivity_wc)) + geom_boxplot(notch = T, varwidth = T) + facet_wrap(~DIDE)
ggplot(data = transcriptdata, mapping = aes(x = DIDE, y = AbsConnectivity_wc)) + geom_boxplot(notch = T, varwidth = T) + facet_wrap(~strata)
ggplot(data = transcriptdata, mapping = aes(x = DIDE, y = AbsBetweenness_wc)) + geom_boxplot(notch = T, varwidth = T) + facet_wrap(~strata) + scale_y_log10()

# Set metazoan as intercept stratum
transcriptdata$strata <- factor(x = transcriptdata$strata, levels = c("Metazoa", "Arthropoda", "Insecta", "Hymenoptera", "Nasonia"))

# Connectivity
lm1 <- glm(formula = AbsConnectivity_wc ~ strata*(DI+DE)+nNodes, data = transcriptdata, family = Gamma)
plot(lm1)
summary(lm1)
summary(model.avg(dredge(lm1)))
write.csv(x = summary(model.avg(dredge(lm1)))$coefmat.full, file = file.path(newdir, "AbsConn_coefs.csv"))
write.csv(x = summary(model.avg(dredge(lm1)))$importance, file = file.path(newdir, "AbsConn_importance.csv"))

# Connectivity, nested
lm1n <- glm(formula = AbsConnectivity_wc ~ strata*(DE/DI)+nNodes, data = transcriptdata, family = Gamma)
plot(lm1n)
summary(lm1n)
summary(model.avg(dredge(lm1n)))

# MAR
lm2 <- glm(formula = MAR_wc ~ nNodes + AbsConnectivity_wc + ClusterCoef_wc + strata*(DI+DE), data = transcriptdata, family = gaussian(link = "logit"))
plot(lm2)
summary(lm2)
summary(model.avg(dredge(lm2)))

# ClusterCoef
lm3 <- glm(formula = ClusterCoef_wc ~ strata*(DI+DE)+nNodes+AbsConnectivity_wc, data = transcriptdata, family = Gamma(link = "logit"))
plot(lm3)
summary(lm3)
summary(model.avg(dredge(lm3)))

# Hubness
lm4 <- glm(formula = Hubness_wc ~ strata*(DI+DE)+nNodes, data = transcriptdata, family = Gamma(link = "logit"))
plot(lm4)
summary(lm4)
summary(model.avg(dredge(lm4)))
write.csv(x = summary(model.avg(dredge(lm4)))$coefmat.full, file = file.path(newdir, "Hub_coefs.csv"))
write.csv(x = summary(model.avg(dredge(lm4)))$importance, file = file.path(newdir, "Hub_importance.csv"))

## Plot results from each density and hub scores per stratum in DI/DE comparisons

########## Density

ConnGLM <- as.data.frame(confint(model.avg(dredge(lm1)), level=.95))
ConnGLM$DEDI <- factor(x = c(NA, "DE","DC",NA, "Unbiased",  "Unbiased",  "Unbiased",  "Unbiased", "DE", "DE", "DE", "DE", "DC", "DC", "DC", "DC"), 
                      levels = c("Unbiased","DC","DE"))
ConnGLM$Estimate <- summary(model.avg(dredge(lm1)))$coefmat.full[,1]
names(ConnGLM)[1:2] <- c("min","max")
ConnGLM$stratum <- factor(x = c(rep(NA, 4), rep(c("Arthropoda", "Insecta", "Hymenoptera", "Nasonia"), 3)), levels = c("Arthropoda", "Insecta", "Hymenoptera", "Nasonia"))

## Erase bottom ticks, set background to white
ggplot(data = ConnGLM[-c(1:4),], mapping = aes(x = DEDI, col = DEDI, ymin = min, ymax = max, y = Estimate)) + 
  geom_pointrange() + facet_grid(.~stratum) + scale_color_manual(values = cbPalette, name = "Sex Bias Type") +
  scale_x_discrete(name = NULL, labels = NULL) + scale_y_continuous(name = "Relative Connectivity") + 
  ggtitle(label = "Connectivity of Nodes from Different Strata") + 
  theme(axis.ticks = element_blank(), panel.border = element_blank(), panel.background = element_rect(fill = "white", colour = NA), panel.grid.major.y = element_line(colour = "grey92"), panel.grid.minor.y = element_line(colour = "grey92", size = 0.25), legend.position = "bottom")
ggsave(filename = file.path(graphdir, "Stratum_AbsConn_pointrange.pdf"), device = "pdf", width = textWidthHeightmm[1], height = textWidthHeightmm[1], units = "mm")

########## Hub Scores

hubGLM <- as.data.frame(confint(model.avg(dredge(lm4)), level=.95))
hubGLM$DEDI <- factor(x = c(NA, "DE","DC",NA, "Unbiased",  "Unbiased",  "Unbiased",  "Unbiased", "DE", "DE", "DE", "DE", "DC", "DC", "DC", "DC"), 
                      levels = c("Unbiased","DC","DE"))
hubGLM$Estimate <- summary(model.avg(dredge(lm4)))$coefmat.full[,1]
names(hubGLM)[1:2] <- c("min","max")
hubGLM$stratum <- factor(x = c(rep(NA, 4), rep(c("Arthropoda", "Insecta", "Hymenoptera", "Nasonia"), 3)), levels = c("Arthropoda", "Insecta", "Hymenoptera", "Nasonia"))

## Erase bottom ticks, set background to white
ggplot(data = hubGLM[-c(1:4),], mapping = aes(x = DEDI, col = DEDI, ymin = min, ymax = max, y = Estimate)) + 
  geom_pointrange() + facet_grid(.~stratum) + scale_color_manual(values = cbPalette, name = "Sex Bias Type") +
  scale_x_discrete(name = NULL, labels = NULL) + scale_y_continuous(name = "Relative Hub Score") + ggtitle(label = "Hub Scores of Nodes from Different Strata") + 
  theme(axis.ticks = element_blank(), panel.border = element_blank(), panel.background = element_rect(fill = "white", colour = NA), panel.grid.major.y = element_line(colour = "grey92"), panel.grid.minor.y = element_line(colour = "grey92", size = 0.25), legend.position = "bottom")
ggsave(filename = file.path(graphdir, "Stratum_HubScores_pointrange.pdf"), device = "pdf", width = textWidthHeightmm[1], height = textWidthHeightmm[1], units = "mm")

ggplot(data = hubGLM[-c(1:4),], mapping = aes(x = stratum, fill = DEDI, ymin = min, ymax = max, y = Estimate, group = DEDI)) + 
  geom_ribbon(alpha = 0.3) + geom_line() +
  scale_fill_manual(values = cbPalette, name = "Sex Bias Type") +
  scale_x_discrete(name = "Stratum") + scale_y_continuous(name = "Relative Hub Score") + theme_bw() + ggtitle(label = "Hub Scores of Nodes from Different Strata")
ggsave(filename = file.path(graphdir, "Stratum_HubScores_ribbon.pdf"), device = "pdf", width = textWidthHeightmm[1], height = textWidthHeightmm[2]/2, units = "mm")

## Correct for main intercepts in each cluster class
hubGLM2 <- hubGLM[-c(1:4),]
hubGLM2[,c(1,2,4)] <- as.matrix(hubGLM2[,c(1,2,4)]) +
  matrix(ncol = 3, byrow = T, c(
    as.numeric(rep(0,12)),
    as.numeric(rep(hubGLM[2,c(1,2,4)],4)),
    as.numeric(rep(hubGLM[3,c(1,2,4)],4))
    ))
hubGLM2 <- cbind(hubGLM2, hubGLM[-c(1:4),c("DEDI","stratum")])

ggplot(data = hubGLM2, mapping = aes(x = DEDI, col = DEDI, ymin = min, ymax = max, y = Estimate)) + 
  geom_pointrange() + facet_grid(.~stratum) + scale_color_manual(values = cbPalette, name = "Sex Bias Type") +
  scale_x_discrete(name = NULL, labels = NULL) + scale_y_continuous(name = "Relative Hub Score") + ggtitle(label = "Hub Scores of Nodes from Different Strata")

ggplot(data = hubGLM2, mapping = aes(x = stratum, fill = DEDI, ymin = min, ymax = max, y = Estimate, group = DEDI)) +   geom_ribbon(alpha = 0.3) + geom_line() +
  scale_fill_manual(values = cbPalette, name = "Sex Bias Type") +
  scale_x_discrete(name = "Stratum") + scale_y_continuous(name = "Relative Hub Score") + theme_bw() + ggtitle(label = "Hub Scores of Nodes from Different Strata")



#####################

# Betweenness may require zero inflated modelling (or maybe poisson)
# using hurdle here
library(pscl)
lm5 <- hurdle(formula = as.integer(Betweenness_wc) ~ strata*(DE+DI)+nNodes, data = transcriptdata, dist = "negbin", zero.dist = "binomial")
plot(residuals(lm5),fitted(lm5))
qqmath(residuals(lm5))
summary(lm5)

lm6 <- zeroinfl(formula = as.integer(Betweenness_wc) ~ strata*(DE+DI)+nNodes|nNodes, data = transcriptdata, dist = "poisson")
qqmath(residuals(lm6))
summary(lm6)

plot(glm(formula = s(AbsBetweenness_wc) ~ strata*(DE+DI)+nNodes, family = Gamma, data = transcriptdata))

library(mgcv)
gam(formula = s(AbsBetweenness_wc) ~ strata*(DE+DI)+nNodes, data = transcriptdata)


## Outdated code snippets
# # DE vs non DE
# DEvnDE <- glm(formula = grepl(pattern = "[m,f]", x = cluster_devsexbias) ~ Density + Centralization + Heterogeneity + nNodes + splicingProp + medianClusterCoef + diameter, family = binomial, data = clusterdata)
# 
# dredge(DEvnDE)[1:20,]
# summary(model.avg(dredge(DEvnDE)))
# 
# bestDEvnDE <- glm(formula = grepl(pattern = "[m,f]", x = cluster_devsexbias) ~ Centralization + splicingProp + medianClusterCoef, family = binomial, data = clusterdata)
# 
# plot(bestDEvnDE)
# summary(bestDEvnDE)
# 
# # DIDE vs rest
# 
# DIDEvnDIDE <- glm(formula = DIDE=="DIDE" ~ Density + Centralization + Heterogeneity + nNodes + splicingProp + medianClusterCoef + diameter, family = binomial, data = clusterdata)
#  
# dredge(DIDEvnDIDE)
# summary(model.avg(dredge(DIDEvnDIDE)))
# 
# bestDIDEvnDIDE <- glm(formula = DIDE=="DIDE" ~ Centralization + splicingProp, family = binomial, data = clusterdata)
# summary(bestDIDEvnDIDE)
# plot(bestDIDEvnDIDE)
# 
# # DIDE vs DE
# DIDEvDE <- glm(formula = DIDE=="DIDE" ~ Density + Centralization + Heterogeneity + nNodes + splicingProp + medianClusterCoef + diameter, family = binomial, data = clusterdata, subset = which(clusterdata$DIDE!="Unbiased"))
# 
# dredge(DIDEvDE)
# summary(model.avg(dredge(DIDEvDE)))
# 
# bestDIDEvDE <- glm(formula = DIDE=="DIDE" ~ Centralization, family = binomial, data = clusterdata, subset = clusterdata$DIDE!="Unbiased")
# summary(bestDIDEvDE)
# plot(bestDIDEvDE)
# 
# # DI and DE vs nonDE
# BiasVsRest <- glm(formula = DIDE!="Unbiased" ~ Density + Centralization + Heterogeneity + nNodes + splicingProp + medianClusterCoef + diameter, family = binomial, data = clusterdata)
# 
# dredge(BiasVsRest)
# summary(model.avg(dredge(BiasVsRest)))
# 
# bestBiasVsRest <- glm(formula = DIDE!="Unbiased" ~ Centralization + splicingProp + medianClusterCoef, family = binomial, data = clusterdata)
# summary(bestBiasVsRest)
# plot(bestBiasVsRest)
# 
# # DI vs non DI
# DIvnDI <- glm(formula = DI ~ Density + Centralization + Heterogeneity + nNodes + splicingProp + medianClusterCoef + diameter, family = binomial, data = clusterdata2)
# 
# dredge(DIvnDI)[1:20,]
# summary(model.avg(dredge(DIvnDI)))
# 
# bestDIvnDI <- glm(formula = DI ~ Centralization + splicingProp, family = binomial, data = clusterdata2)
# 
# plot(bestDIvnDI)
# summary(bestDIvnDI)
# 
# ## Check whether DIDE has higher centralization than DE
# 
# CentraLM <- lm(formula = Centralization ~ DI + DE + nNodes + splicingProp + Density, data = clusterdata2)
# 
# plot(CentraLM)
# summary(CentraLM)
# 
# summary(model.avg(dredge(CentraLM)))
# dredge(CentraLM)
# 
# #### Use DA
# library(MASS)
# # linear
# DA_all <- lda(DIDE ~ Density + Centralization + Heterogeneity + nNodes + splicingProp + medianClusterCoef + diameter, data = clusterdata, CV = T)
# 
# # Assess the accuracy of the prediction
# # percent correct for each category of G
# ct <- table(clusterdata$DIDE, DA_all$class)
# diag(prop.table(ct, 1))
# # total percent correct
# sum(diag(prop.table(ct)))
# 
# DA_all <- lda(DIDE ~ Density + Centralization + Heterogeneity + nNodes + splicingProp + medianClusterCoef + diameter, data = clusterdata, CV = F)
# DA_all
# plot(DA_all)
# 
# partimat(DIDE ~ Density + Centralization + Heterogeneity + nNodes + splicingProp + medianClusterCoef + diameter, data = clusterdata)
# 
# # quadratic
# DA_all <- qda(DIDE!="Unbiased" ~ Density + Centralization + Heterogeneity + nNodes + splicingProp + medianClusterCoef + diameter, data = clusterdata)
