####### Do different strata have different connectivities?
library(ggplot2)
newdir<-file.path(getwd(),"Output/Wasp_Phylostratigraphy_HubScores")

transcriptdata <- read.csv(file = "./Output/Results_compiler_CCRE/transcriptdata_full.csv")[,-1]
transcriptdata$Stratum = factor(
  x = transcriptdata$Stratum, 
  levels = c("Hymenoptera", "Apocrita", "Chalcid", "Pteromalid", "Nasonia"))

# Do nodes from different strata have different connectivities?
ggplot(data = transcriptdata, mapping = aes(x = Stratum, y = AbsConnectivity_wc)) +
  geom_boxplot(notch = T, varwidth = T) #YES!
# Do nodes from different strata have different Hub Scores?
ggplot(data = transcriptdata, mapping = aes(x = Stratum, y = Hubness_wc)) + 
  geom_boxplot(notch = T, varwidth = T) #YES!

# Add information about sex-bias
clusterdata <- read.csv(file = "./Output/Results_compiler_CCRE/clusterdata_full.csv", row.names = 1)
clusterdata$DIDE <- as.factor((1-is.na(clusterdata$DiffIntegrated))+2*grepl(pattern = "[m,f]", clusterdata$cluster_devsexbias))
clusterdata$DIDE <- factor(x = clusterdata$DIDE, labels = c("Unbiased","DI","DE","DIDE"))
clusterdata$DE <- grepl(pattern = "[m,f]", x = clusterdata$cluster_devsexbias)
clusterdata$DI <- clusterdata$DIDE=="DIDE"
clusterdata <- clusterdata[,c("clusterID", "cluster_devsexbias", "DiffIntegrated", "DIDE", "Density", "Centralization", "Heterogeneity", "nNodes", "excessSpl", "excessDup", "medianClusterCoef", "diameter", "DI", "DE")]

# Merge sex-bias information with stratum information
transcriptdata <- merge(transcriptdata, clusterdata[,c("clusterID","DIDE","DI","DE", "nNodes")])
# transcriptdata$grouper <- paste(transcriptdata$Stratum, transcriptdata$DIDE)

# # Remove nodes with no assigned stratum
# transcriptdata <- droplevels(transcriptdata[-which(is.na(transcriptdata$Stratum)),])

# Filter only representative nodes
transcriptdata <- transcriptdata[which(transcriptdata$Representative_Node==T|is.na(transcriptdata$Representative_Node)),]
# Remove nodes that map to strata above Hymenoptera according to Saxton
transcriptdata <- transcriptdata[which(!transcriptdata$strata%in%c("Metazoa", "Arthropod", "Insect")),]
transcriptdata <- droplevels(transcriptdata)

# Count nodes left
summary(as.factor(transcriptdata$Stratum))

# Merge Apocrita in Hymenoptera, NAs in Nasonia
transcriptdata$Stratum <- as.character(transcriptdata$Stratum)
transcriptdata$Stratum <- gsub(pattern = "Apocrita", replacement = "Hymenoptera", x = transcriptdata$Stratum)
# transcriptdata$Stratum <- gsub(pattern = "Nasonia", replacement = "Pteromalid", x = transcriptdata$Stratum)
transcriptdata$Stratum <- ifelse(is.na(transcriptdata$Stratum)&transcriptdata$strata=="Wasp", "Nasonia", transcriptdata$Stratum) # Assign all genes matching the wasp stratum in Saxton and missing matches in Trichogramma to Nasonia-specific stratum
transcriptdata$Stratum <- factor(x = transcriptdata$Stratum, levels = c("Hymenoptera", "Chalcid", "Pteromalid", "Nasonia"))

table(transcriptdata$Stratum, useNA = 'ifany')

transcriptdata <- transcriptdata[-which(is.na(transcriptdata$Stratum)),]

# Plot connectivities per stratum vs sex-bias
ggplot(data = transcriptdata, mapping = aes(x = DIDE, y = AbsConnectivity_wc, col = DIDE)) +
  geom_boxplot(notch = T, varwidth = F) +
  facet_grid(.~Stratum)
# Plot Hub scores per stratum vs sex-bias
ggplot(data = transcriptdata, mapping = aes(x = DIDE, y = Hubness_wc, col = DIDE)) +
  geom_boxplot(notch = T, varwidth = F) +
  facet_grid(.~Stratum)

# DI and DE have higher connectivity and hubness across all strata BUT
# DIDE have lower Hub/con for Pteromalid-stratum genes (sample size 167!)


### Model as GLMs
library(lme4)
library(MuMIn)
options(na.action = "na.fail")


# Connectivity
lm1 <- glm(
  formula = AbsConnectivity_wc ~ Stratum*(DE+DI)+log(nNodes), 
  data = transcriptdata, 
  family = Gamma(link = 'logit')
)
# plot(lm1)
summary(lm1)
summary(model.avg(dredge(lm1)))
write.csv(x = summary(model.avg(dredge(lm1)))$coefmat.full, file = file.path(newdir, "AbsConn_coefs.csv"))
write.csv(x = summary(model.avg(dredge(lm1)))$importance, file = file.path(newdir, "AbsConn_importance.csv"))

# Hubness
lm4 <- glm(
  formula = Hubness_wc ~ Stratum*(DE+DI)+ log(nNodes), 
  data = transcriptdata, 
  family = Gamma(link = 'logit')
)
# plot(lm4)
summary(lm4)
summary(model.avg(dredge(lm4)))
write.csv(x = summary(model.avg(dredge(lm4)))$coefmat.full, file = file.path(newdir, "Hub_coefs.csv"))
write.csv(x = summary(model.avg(dredge(lm4)))$importance, file = file.path(newdir, "Hub_importance.csv"))

# No overall effect of DE/DC on connectivity/hubness
# Compared to Hymenoptera, Nasonia has significantly lower connectivity and hubness (same Effect Sizes)
# DE clusters counterbalance loss of hubness/connectivity for Nasonia stratum
# No support for DI
# Cannot test for DE:DI interactions since no clusters are DI and not DE