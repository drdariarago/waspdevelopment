### Integrate hymenopteran phylostratigraphic data with DE_DC status and network parameters
rm(list=ls())

library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
textWidthHeightmm <- c(152.73, 216.00)
source(file = "./Scripts/multiplot.R")
# The colorblind friendly palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

graphdir<-file.path(getwd(), "Graphics/Wasp_Phylostratigraphy_Networks")
dir.create(graphdir)

# Annotate DE and DI clusters
clusterdata = read.csv(file = "./Output/Results_compiler_CCRE/clusterdata_full.csv", row.names = 1)
sexDE_clusters = clusterdata$clusterID[which(clusterdata$cluster_devsexbias!=".....")]
sexDI_clusters = clusterdata$clusterID[which(!is.na(clusterdata$DiffIntegrated))]

# Import dataset, remove older strata, collapse sparse taxonomic strata and convert to counts
transcriptdata <- read.csv(file = "./Output/Results_compiler_CCRE/transcriptdata_full.csv", row.names = 1)

transcriptdata <- transcriptdata[which(!transcriptdata$strata%in%c("Metazoa", "Arthropod", "Insect")),]
transcriptdata <- droplevels(transcriptdata)

transcriptdata$Stratum <- as.character(transcriptdata$Stratum)
transcriptdata$Stratum <- gsub(pattern = "Apocrita", replacement = "Hymenoptera", x = transcriptdata$Stratum)
transcriptdata$Stratum <- gsub(pattern = "Nasonia", replacement = "Pteromalid", x = transcriptdata$Stratum)
transcriptdata$Stratum <- factor(x = transcriptdata$Stratum, levels = c("Hymenoptera", "Chalcid", "Pteromalid"))

transcriptdata <- transcriptdata[which(is.na(transcriptdata$Stratum)==F),]

taxprop_cluster <- as.data.frame(table(transcriptdata$clusterID, transcriptdata$Stratum, useNA = "ifany"))

names(taxprop_cluster) <- c("clusterID","Stratum", "Counts")

# Annotate DE and DI status
taxprop_cluster$sexDE <- taxprop_cluster$clusterID%in%sexDE_clusters
#taxprop_cluster$sexDE <- ifelse(taxprop_cluster$sexDE, "DE", "N")
taxprop_cluster$sexDI <- taxprop_cluster$clusterID%in%sexDI_clusters
#taxprop_cluster$sexDI <- ifelse(taxprop_cluster$sexDI, "DI", "N")
taxprop_cluster$DEDI <- paste(taxprop_cluster$sexDE, taxprop_cluster$sexDI, sep = "")
taxprop_cluster <- recast(taxprop_cluster, formula = clusterID + sexDE + sexDI + DEDI ~ Stratum)

# Convert to proportions
taxprop_cluster <- droplevels(na.exclude(taxprop_cluster))
taxprop_cluster$total <- apply(X = taxprop_cluster[5:7], MARGIN = 1, FUN = sum)
taxprop_cluster <- taxprop_cluster[,c("clusterID", "sexDE", "sexDI", "DEDI", "total", "Hymenoptera", "Chalcid", "Pteromalid")]
# Normalize by networkwide proportions
taxmeans <- apply(taxprop_cluster[,5:8], 2, sum)
taxmeans <- taxmeans[2:4]/taxmeans[1]

taxprop_cluster_2 <- taxprop_cluster
taxprop_cluster_2[,6:8] <- taxprop_cluster_2[,6:8]/taxprop_cluster_2$total
taxprop_cluster_2 <- taxprop_cluster_2[,-grep(pattern = "total", x = names(taxprop_cluster_2))]
taxprop_cluster_2[,5:7] <- apply(X = taxprop_cluster_2[,5:7], MARGIN = 1, FUN = function(x) {x/taxmeans})

taxprop_cluster_2 <- recast(data = taxprop_cluster_2, formula = clusterID + sexDE + sexDI + DEDI + variable ~ .)
names(taxprop_cluster_2)[5:6] <- c("Stratum", "Proportion")

taxprop_cluster_2$DEDI <- factor(
  x = taxprop_cluster$DEDI,
  levels = c("FALSEFALSE", "FALSETRUE", "TRUEFALSE", "TRUETRUE"),
  labels = c("Unbiased", "DI", "DE", "DEDI")
)
taxprop_cluster_2$Stratum <- factor(
  x = taxprop_cluster_2$Stratum, 
  levels = c("Hymenoptera", "Chalcid", "Pteromalid"),
  labels = c("Hymenoptera", "Chalcid", "Pteromalid")
)

# Plot as boxplot
ggplot(data = taxprop_cluster_2, 
  mapping = aes(x = DEDI!="Unbiased", y = Proportion, col = DEDI=="Unbiased")) +
  geom_boxplot(notch = T) + facet_grid(. ~ Stratum) +
  geom_hline(yintercept = 1) +
  ggtitle(label = "Taxonomic Stratum Enrichment in Clusters with \nDifferential Expression and Differential Correlation") +
  ylab(label = "Excess Proportion of Nodes from Stratum") + 
  xlab(label = "") + 
  theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
  #  scale_color_brewer(type = "qual", name = "Differentially\nExpressed") + 
  scale_linetype_discrete(name = "Differentially\nCorrelated") + 
  scale_y_continuous(limits = c(0.1,4.5)) + 
  theme_bw()

# Plot as violinplot
ggplot(data = taxprop_cluster_2, 
  mapping = aes(x = DEDI!="Unbiased", y = Proportion, col = DEDI!="Unbiased")) +
  geom_violin() + facet_grid(. ~ Stratum) +
  geom_point(position = 'jitter', alpha = '0.2') +
  geom_hline(yintercept = 1) +
  ggtitle(label = "Taxonomic Stratum Enrichment in Clusters with \nDifferential Expression and Differential Correlation") +
  ylab(label = "Excess Proportion of Nodes from Stratum") + 
  xlab(label = "") + 
  theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
  #  scale_color_brewer(type = "qual", name = "Differentially\nExpressed") + 
  scale_linetype_discrete(name = "Differentially\nCorrelated") + 
  scale_y_continuous(limits = c(0.1,4.5)) + 
  theme_bw()

# Plot only unbiased, DE and DEDI levels
ggplot(data = droplevels(taxprop_cluster_2[which(taxprop_cluster_2$DEDI!='DI'),]), 
  mapping = aes(x = DEDI, y = Proportion, col = DEDI)) +
  geom_boxplot(notch = T) + facet_grid(. ~ Stratum) +
  geom_hline(yintercept = 1) +
  ggtitle(label = "Taxonomic Stratum Enrichment in Clusters with \nDifferential Expression and Differential Correlation") +
  ylab(label = "Excess Proportion of Nodes from Stratum") + 
  xlab(label = "") + 
  theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
  scale_color_brewer(type = "qual", name = "Sex-Bias Type") + 
  scale_linetype_discrete(name = "Differentially\nCorrelated") + 
  scale_y_continuous(limits = c(0,4)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
    )
ggsave(filename = file.path(graphdir, "Wasp_Strata.pdf"), width = textWidthHeightmm[2], height = textWidthHeightmm[1], units = 'mm')


# Plot only unbiased and DE
ggplot(data = droplevels(taxprop_cluster_2), 
  mapping = aes(x = sexDE, y = Proportion, col = sexDE)) +
  geom_boxplot(notch = T) + facet_grid(. ~ Stratum) +
  geom_hline(yintercept = 1) +
  ggtitle(label = "Taxonomic Stratum Enrichment in Clusters with \nDifferential Expression and Differential Correlation") +
  ylab(label = "Excess Proportion of Nodes from Stratum") + 
  xlab(label = "") + 
  theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
  scale_color_brewer(type = "qual", name = "Sex-Bias Type") + 
  scale_linetype_discrete(name = "Differentially\nCorrelated") + 
  scale_y_continuous(limits = c(0,4)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave(filename = file.path(graphdir, "Wasp_Strata_DE.pdf"), width = textWidthHeightmm[2], height = textWidthHeightmm[1], units = 'mm')


# # Collapse DI and DEDI
# taxprop_cluster_2$DEDI = as.character(taxprop_cluster_2$DEDI)
# taxprop_cluster_2$DEDI = as.factor(ifelse(taxprop_cluster_2$DEDI=="DI", "DEDI", taxprop_cluster_2$DEDI))
# 
# # Visualize trend over time for each class of clusters
# ggplot(data = taxprop_cluster_2, mapping = aes(x = as.numeric(Stratum), y = Proportion, col = DEDI!="Unbiased")) +
#   geom_point(position = 'jitter', alpha = '.4') + geom_smooth() + 
#   scale_y_continuous(limits = c(0.1,4))
# 
# ggplot(data = taxprop_cluster_2, mapping = aes(x = as.numeric(Stratum), y = Proportion, col = DEDI)) +
#   geom_point(position = 'jitter', alpha = '.4') + geom_smooth() + 
#   scale_y_continuous(limits = c(0.1,4))
# 
# 
# # Looks like we have a very mild enrichment in sex DE during the Chalcid/Pteromalid strata
# # Must contrast with previous analyses so graph is systematic

# Analyze as GLMMs
head(taxprop_cluster)
taxprop_cluster_3 <- recast(
  data = taxprop_cluster[,c("clusterID", "sexDE", "sexDI", "total", "Hymenoptera", "Chalcid", "Pteromalid")], 
  formula = clusterID + sexDE + sexDI + total + variable ~ .,
  measure.var = 5:7
)
head(taxprop_cluster_3)
names(taxprop_cluster_3)[5:6] <- c("Stratum", "Counts")

library(lme4)
GLMMstrata <- glm(
  formula = Counts ~  Stratum * (sexDE + sexDI) + scale(total),
  data = taxprop_cluster_3,
  family = poisson
)
summary(GLMMstrata)

library(MuMIn)
options(na.action = na.fail)

GLMMstrataSet <- dredge(GLMMstrata)
model.sel(GLMMstrataSet)
summary(model.avg(GLMMstrataSet))
confint(model.avg(GLMMstrataSet))

# Moderate support for enrichment of chalcid genes in DE clusters

#### Unused code snippets

# taxprop_cluster$DEDI <- factor(
#   x = taxprop_cluster$DEDI, 
#   levels = c("UnbiasedUnbiased", "DEUnbiased", "UnbiasedDI", "DEDI"), 
#   labels = c("Unbiased", "DE", "DI", "DEDI")
# )
# taxprop_cluster <- taxprop_cluster[,c("clusterID", "Stratum", "Counts", "DEDI")]
# #taxprop_cluster <- recast(data = taxprop_cluster, formula = ClusterID + DEDI ~ Stratum)

# Reformat to Unbiased/DE/DEDI

# Convert to proportions and normalize to transcriptome-wide proportions
# ClusterSizes = ddply(taxprop_cluster, .variables = .(clusterID), .fun = summarize, 
#   TotalSize = sum(Counts))
# taxprop_cluster = merge(taxprop_cluster, ClusterSizes)
# taxprop_cluster$Proportions <- taxprop_cluster$Counts/taxprop_cluster$TotalSize
# 
# StratumTotals = ddply(taxprop_cluster, .variables = .(Stratum), .fun = summarize, 
#   TotalNodes = sum(Counts))
# StratumTotals$ExpectedProportions = StratumTotals$TotalNodes/sum(StratumTotals$TotalNodes)
# taxprop_cluster = merge(taxprop_cluster, StratumTotals)
# taxprop_cluster$NormalizedProportion = taxprop_cluster$Counts/taxprop_cluster$ExpectedProportions
