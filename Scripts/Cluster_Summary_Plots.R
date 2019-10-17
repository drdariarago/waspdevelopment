### Plots for final paper
### Summary statistics on sexbiased nodes
# Load packages
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
# Clean workspace
rm(list=ls())
# Save text size (article KOMA)
textWidthHeightmm <- c(152.73, 216.00)
# Load multiplot function
source(file = "./Scripts/multiplot.R")
# Load colorblind friendly palette with grey
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Load colorblind palette for unbiased, female, male, both
cbPalette4 <- c("#999999", "#D55E00", "#0072B2", "#CC79A7")
# Create output path
newdir <- file.path(getwd(), "Output/SummaryPlots")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/SummaryPlots")
dir.create(graphdir)

## Load dataset
transcriptdata <- read.csv(file = "./Output/Results_compiler_CCRE/transcriptdata_full.csv")[,-1]

# Return devconflict nodes
transcriptdata[grep(pattern = "f.*m|m.*f", x = transcriptdata$devsexbias),]


# Rapid lookup
transcriptdata[which(transcriptdata$clusterID=="navajowhite3"),c("Hubness_wc","RelHubness_wc","RelBetweenness_wc","AbsBetweenness_wc","eigenexonID","nodeID")]

## Shape data for plotting sexDE clusters over stages
clusterGLMsummaries <- read.csv(file = "./Output/DEcluster_analysis_CCRE/clusterGLMsummaries.csv")[,-1]
clusterGLMsummaries$stage <- factor(x = clusterGLMsummaries$stage, levels = c("emb10","emb18","lar51","pupyel","adult"))
clusterGLMsummaries$sex <- sign(clusterGLMsummaries$coeff*(clusterGLMsummaries$fdr<0.05))
clusterGLMsummaries$sex <- factor(ifelse(clusterGLMsummaries$contrast=="Sex", clusterGLMsummaries$sex, NA), labels = c("Female","Unbiased","Male"))
clusterGLMsummaries$clusterID <- substring(text = clusterGLMsummaries$clusterID, first = 3)
# and plot
ggplot(data = clusterGLMsummaries[which(clusterGLMsummaries$contrast=="Sex"&clusterGLMsummaries$sex!="Unbiased"),], mapping = aes(x = stage, fill = sex))+geom_bar(position = "dodge")+theme_bw()+ggtitle(label = "Number of Sex-Biased \nClusters across Development")+scale_fill_manual(values = cbPalette4[-1], name = "Sex Bias") + theme(legend.position = "bottom")
ggsave(filename = file.path(graphdir, "clusterBiasvsStage.pdf"), device = "pdf", width = textWidthHeightmm[1]*.7, height = textWidthHeightmm[2]/3, units = "mm")

# add cluster sizes
clusterdata <- read.csv(file = "./Output/Results_compiler_CCRE/clusterdata_full.csv")[,-1]
clusterGLMsummaries <- merge(x = clusterGLMsummaries, y = clusterdata[,c("clusterID","nTranscripts")], by = "clusterID", all.x = T)

# and plot
ggplot(data = clusterGLMsummaries[which(clusterGLMsummaries$contrast=="Sex"&clusterGLMsummaries$sex!="Unbiased"),], 
       mapping = aes(x = stage, fill = sex, alpha = log10(nTranscripts))) +
  geom_bar(position = "dodge")+theme_bw()+ggtitle(label = "Number of Sex-Biased \nClusters across Development")+
  scale_fill_manual(values = cbPalette4[-1], name = "Sex Bias") + 
  scale_alpha_continuous() + 
  theme(legend.position = "bottom")
ggsave(filename = file.path(graphdir, "clusterBiasvsStage.pdf"), device = "pdf", width = textWidthHeightmm[1]*.7, height = textWidthHeightmm[2]/3, units = "mm")


## Which stage DE patterns are more frequent in male and femalebiased clusters?
clusterdata <- read.csv(file = "./Output/Results_compiler_CCRE/clusterdata_full.csv")[,-1]
sexDE <- ifelse(grepl(pattern = "f", clusterdata$cluster_devsexbias), ifelse(grepl(pattern = "m", clusterdata$cluster_devsexbias), "Both", "Female"), ifelse(grepl(pattern = "m", clusterdata$cluster_devsexbias), "Male", "Unbiased"))
sexDE <- factor(x = sexDE, levels = c("Unbiased", "Female", "Male", "Both"))
sexDE <- data.frame(clusterID = clusterdata$clusterID, sexbias = sexDE)
# sexDE <- sexDE[which(sexDE$sexbias!="Both"),]
devbiascounts <- clusterGLMsummaries[which(clusterGLMsummaries$contrast=="Stage"&clusterGLMsummaries$fdr<.05&clusterGLMsummaries$coeff>0),]
devbiascounts <- droplevels(merge(devbiascounts, sexDE))

ggplot(data = devbiascounts, mapping = aes(x = stage, fill = sexbias, group = sexbias)) + geom_bar(stat = "count", position = "dodge") + theme_bw() + ggtitle(label = "Developmental Differential Expression of Clusters\n with Male and Female Differential Expression") + scale_y_continuous(name = "Number of Differentially\n Expressed Clusters") + scale_x_discrete(name = "Stage") + scale_fill_manual(values = cbPalette4, name = "Differential Expression") + theme(legend.position = "bottom")
ggsave(filename = file.path(graphdir, "devsexbias_counts.pdf"), device = "pdf", width = textWidthHeightmm[1], height = textWidthHeightmm[2]/2, units = "mm")

# Normalize by number of overall stageDE clusters per stage
devbiascounts_2 <- table(devbiascounts$stage, devbiascounts$sexbias)
devbiascounts_2 <- devbiascounts_2/apply(X = devbiascounts_2, MARGIN = 1, sum)
devbiascounts_2 <- melt(devbiascounts_2)
names(devbiascounts_2) <- c("stage", "sexbias", "proportion")
devbiascounts_2$sexbias <- factor(x = devbiascounts_2$sexbias, levels = c("Unbiased", "Female", "Male", "Both"))
devbiascounts_2$stage <- factor(x = devbiascounts_2$stage, levels = c("emb10", "emb18", "lar51", "pupyel", "adult"))
# And plot
ggplot(data = devbiascounts_2, mapping = aes(x = stage, y = proportion, fill = sexbias)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = cbPalette4, name = "Differential Expression") + theme(legend.position = "bottom") + theme_bw() + theme(legend.position = "bottom") + ggtitle(label = "Differential Expression Class \nof Stage Biased Clusters")
ggplot(data = devbiascounts_2, mapping = aes(x = stage, y = proportion, fill = sexbias)) + geom_bar(stat = "identity") + scale_fill_manual(values = cbPalette4, name = "Differential Expression") + theme(legend.position = "bottom") + theme_bw() + theme(legend.position = "bottom") + ggtitle(label = "Differential Expression Class \nof Stage Biased Clusters")

## Add total number of stageDE clusters as width?
## Check clusters which are "expressed"?
## Multiply by number of transcripts?
devbiascounts <- merge(devbiascounts, clusterdata[,c("clusterID", "nTranscripts")], by = "clusterID")

devbiascounts_3 <- ddply(.data = devbiascounts, .variables = .(stage, sexbias), summarize, nTranscripts = sum(nTranscripts))

ggplot(data = devbiascounts_3, mapping = aes(x = stage, fill = sexbias, y = nTranscripts)) + geom_bar(stat = "identity", position = "dodge") + theme_bw() + ggtitle(label = "Developmental Differential Expression of Clusters\n with Male and Female Differential Expression") + scale_y_continuous(name = "Number of Differentially\n Expressed Clusters") + scale_x_discrete(name = "Stage") + scale_fill_manual(values = cbPalette4, name = "Differential Expression") + theme(legend.position = "bottom")

devbiascounts_ref <- ddply(.data = devbiascounts, .variables = .(stage), summarize, nTranscripts = sum(nTranscripts))
devbiasprop <- merge(devbiascounts_3, devbiascounts_ref, by = "stage", suffixes = c("","_ref"))
devbiasprop$prop <- devbiasprop$nTranscripts/devbiasprop$nTranscripts_ref

ggplot(data = devbiasprop, mapping = aes(x = stage, fill = sexbias, y = prop)) + geom_bar(stat = "identity", position = "dodge") + theme_bw() + ggtitle(label = "Developmental Differential Expression of Nodes \nin Clusters with Male and Female Differential Expression") + scale_y_continuous(name = "Proportion of Differentially\n Expressed Nodes") + scale_x_discrete(name = "Stage") + scale_fill_manual(values = cbPalette4, name = "Differential Expression") + theme(legend.position = "bottom")
ggsave(filename = file.path(graphdir, "devsexbias_nodeprops.pdf"), device = "pdf", width = textWidthHeightmm[1], height = textWidthHeightmm[2]/2, units = "mm")
## Early DE genes are dominated by female-biased transcripts (possibly because of maternal transcript inheritance)
## Other stages show similar proportions between both sexes, possibly disproving the more "specific" pattern of male-biased genes

ggplot(data = devbiasprop, mapping = aes(x = stage, fill = sexbias, y = prop)) + geom_bar(stat = "identity") + theme_bw() + ggtitle(label = "Developmental Differential Expression of Nodes in Clusters\n with Male and Female Differential Expression") + scale_y_continuous(name = "Proportion of Differentially\n Expressed Clusters") + scale_x_discrete(name = "Stage") + scale_fill_manual(values = cbPalette4, name = "Differential Expression") + theme(legend.position = "bottom")



## Add DI coefficients and plot in series of 2D graphs
# Negative terms indicate female-bias in the parameter, positive terms indicate male-bias
# import and format data as cluster+stage~coeff+fdr
DI_data <- read.csv("./Output/sex_modularity_GLMs_CCRE/fdr_perm_dWithin.csv")[,-1]
DI_data <- DI_data[,c("coefWithinFactornames",".id","Estimate","lFDR_adjusted")]
DI_data$stage <- substring(text = str_extract(string = DI_data$coefWithinFactornames, pattern = "^[^:]*"), first = 6)
DI_data$stage <- factor(x = DI_data$stage, levels = c("emb10","emb18","lar51","pupyel","adult"))
DI_data <- DI_data[,c(".id","stage","Estimate","lFDR_adjusted")]
names(DI_data) <- c("clusterID","stage","coeff","fdr")
# merge with reduced clusterDE data
DIDE_data <- merge(x = clusterGLMsummaries[which(clusterGLMsummaries$contrast=="Sex"),c("clusterID","stage","fdr","coeff")], 
                   y = DI_data, 
                   by = c("clusterID","stage"),
                   suffixes = c("_DE","_DI"))
DIDE_data$sex_DE <- factor(sign((DIDE_data$fdr_DE<0.05)*DIDE_data$coeff_DE), labels = c("Female","Unbiased","Male"))
DIDE_data$sex_DI <- factor(sign((DIDE_data$fdr_DI<0.1)*DIDE_data$coeff_DI), labels = c("Female","Unbiased","Male"))
DIDE_data$minfdr <- pmin(DIDE_data$fdr_DE, DIDE_data$fdr_DI)
DIDE_data$DIDE <- paste(DIDE_data$fdr_DI<0.1, DIDE_data$fdr_DE<0.05)
DIDE_data$DIDE <- factor(x = DIDE_data$DIDE, labels = c("No Bias", "DE Only", "DC Only", "DE and DC"))
DIDE_data$stage <- factor(x = DIDE_data$stage, labels = c("Embryo, Early", "Embryo, Late", "Larva", "Pupa", "Adult"))
# Add cluster data
clusterdata <- read.csv(file = "./Output/Results_compiler_CCRE/clusterdata_full.csv")[,-1]
DIDE_data <- merge(DIDE_data, clusterdata[,c("clusterID","nTranscripts","Density","Centralization","Heterogeneity")])

## 2D plot of signed effect sizes for DE and DC
ggplot(data = DIDE_data, mapping = aes(x = coeff_DE, y = coeff_DI, col = DIDE, alpha = DIDE=="No Bias")) + geom_point() + facet_wrap(~stage) + theme_bw() + scale_color_manual(values = cbPalette, name = "Sex-Bias Type") + scale_alpha_manual(values = c(.6,.1), guide = F) + scale_x_continuous(name = "Differential Expression Coefficient") + scale_y_continuous(name = "Differential Correlation Coefficient") + theme(legend.position = c(.8,.2), legend.key.size = unit(3, "mm")) + ggtitle(label = "Differential Expression and Differential\n Correlation Coefficients within Stages")
  # guides(colour = guide_legend(nrow = 2, override.aes=list(size=5)), size = guide_legend(nrow = 2)) # Add to plot to change legend to two columns
ggsave(filename = file.path(graphdir, "effectsizes_DIDE.pdf"), device = "pdf", width = textWidthHeightmm[1], height = textWidthHeightmm[2]/2, units = "mm")

## Use cluster names as dots
ggplot(data = DIDE_data, mapping = aes(x = coeff_DE, y = coeff_DI, col = DIDE, alpha = DIDE=="No Bias", label = clusterID)) + geom_text() + facet_wrap(~stage) + theme_bw() + scale_color_manual(values = cbPalette, name = "Sex-Bias Type") + scale_alpha_manual(values = c(.6,.1), guide = F) + scale_x_continuous(name = "Differential Expression Coefficient") + scale_y_continuous(name = "Differential Correlation Coefficient") + theme(legend.position = c(.8,.2), legend.key.size = unit(3, "mm")) + ggtitle(label = "Differential Expression and Differential\n Correlation Coefficients within Stages")
ggsave(filename = file.path(graphdir, "effectsizes_DIDE_clusternames.pdf"), device = "pdf", width = textWidthHeightmm[1], height = textWidthHeightmm[2]/2, units = "mm")

## q-value 2d plot
ggplot(data = DIDE_data, mapping = aes(x = -log10(fdr_DE)*sign(coeff_DE), y = -log10(fdr_DI)*sign(coeff_DI), col = DIDE))+geom_point()+facet_wrap(~stage)+theme_bw()+geom_hline(yintercept=log10(.1))+geom_vline(xintercept=log10(.05))+geom_hline(yintercept=-log10(.1))+geom_vline(xintercept=-log10(.05))+scale_color_manual(values = cbPalette, name = "Sex-Bias Type") + scale_x_continuous(name = "Differential Expression fdr") + scale_y_continuous(name = "Differential Correlation fdr") + theme(legend.position = c(.8,.2), legend.key.size = unit(3, "mm")) + ggtitle(label = "Differential Expression and Differential\n Correlation False Discovery Rates within Stages")
ggsave(filename = file.path(graphdir, "qvalues_DIDE.pdf"), device = "pdf", width = textWidthHeightmm[1], height = textWidthHeightmm[2]/2, units = "mm")

## Histogram of cluster counts across stages for DIDE
ggplot(data = DIDE_data[which(DIDE_data$DIDE!="No Bias"),], mapping = aes(x = stage, fill = DIDE))+geom_bar(position = "dodge")+theme_bw()+ggtitle(label = "Number of Sex-Biased Clusters across Development")+scale_fill_manual(values = cbPalette[-1], name = "Sex Bias")

# ## Histogram of male vs female DI and DE clusters
# ggplot(data = DIDE_data[which(DIDE_data$sex_DE!="Unbiased"),], mapping = aes(x = stage, fill = sex_DE))+geom_bar(position = "dodge")+theme_bw()+ggtitle(label = "Number of Sex-Biased\n Clusters across Development")+scale_fill_manual(values = cbPalette[-1], name = "Sex Bias")
# ggplot(data = DIDE_data[which(DIDE_data$sex_DI!="Unbiased"),], mapping = aes(x = stage, fill = sex_DI))+geom_bar(position = "dodge")+theme_bw()+ggtitle(label = "Number of Sex-Biased\n Clusters across Development")+scale_fill_manual(values = cbPalette[-1], name = "Sex Bias")

### Check proportion of transcription vs splicing and dupl vs single-copy nodes in sexDE clusters
# Import data
transcriptdata <- read.csv(file = "./Output/Results_compiler_CCRE/transcriptdata_full.csv")[,-1]
clusterdata <- read.csv(file = "./Output/Results_compiler_CCRE/clusterdata_full.csv")[,-1]
# Add annotation
sexDE_clusters <- clusterdata$clusterID[grep(pattern = "[m,f]", clusterdata$cluster_devsexbias)]
sexDI_clusters <- clusterdata$clusterID[which(x = is.na(clusterdata$DiffIntegrated)==F)]

## Calculate counts of transcription vs splicing nodes
transprop_cluster <- as.data.frame(table(transcriptdata$clusterID, transcriptdata$event_type))
names(transprop_cluster) <- c("ClusterID","EventType","Counts")
transprop_cluster$sexDE <- transprop_cluster$ClusterID%in%sexDE_clusters
transprop_cluster$sexDI <- transprop_cluster$ClusterID%in%sexDI_clusters
transprop_cluster$DIDE <- paste(transprop_cluster$ClusterID%in%sexDI_clusters, transprop_cluster$ClusterID%in%sexDE_clusters)
transprop_cluster$DIDE <- factor(x = transprop_cluster$DIDE, labels = c("No Bias", "DE Only", "DC Only", "DE and DC"))
transprop_cluster <- recast(data = transprop_cluster, formula = ClusterID + sexDE + sexDI + DIDE ~ EventType)
# Plot as boxplot
transprop_boxplot <- ggplot(data = transprop_cluster, mapping = aes(x = DIDE, y = splicing/(splicing+transcription), group = DIDE, col = DIDE)) + geom_boxplot(varwidth = T) + theme_bw() + ylab(label = "Proportion of\nSplicing Nodes") + xlab(label = "") + theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom") + scale_color_manual(values = cbPalette, name = "Sex-Bias Type") + scale_linetype_discrete(name = "Sex-Bias Type") + geom_point(position = position_jitter(width = 0.1, height = 0)) + scale_y_continuous(limits = c(0,1)) + guides(col = guide_legend(ncol = 2))
transprop_boxplot

## Calculate counts of single vs multicopy OGs (orthodb8)
orthoprop_cluster <- as.data.frame(table(transcriptdata$clusterID, transcriptdata$Copynumber>1, useNA = "always"))
names(orthoprop_cluster) <- c("ClusterID","Orthologs","Counts")
orthoprop_cluster$sexDE <- orthoprop_cluster$ClusterID%in%sexDE_clusters
orthoprop_cluster$sexDI <- orthoprop_cluster$ClusterID%in%sexDI_clusters
orthoprop_cluster$DIDE <- paste(orthoprop_cluster$ClusterID%in%sexDI_clusters, orthoprop_cluster$ClusterID%in%sexDE_clusters)
orthoprop_cluster <- recast(data = orthoprop_cluster, formula = ClusterID + sexDE + sexDI + DIDE ~ Orthologs)
names(orthoprop_cluster)[5:7] <- c("Single","Multi","NoOG")
orthoprop_cluster$DIDE <- factor(x = orthoprop_cluster$DIDE, labels = c("No Bias", "DE Only", "DC Only", "DE and DC"))
orthoprop_cluster <- droplevels(na.exclude(orthoprop_cluster))
# Normalize by total number of nodes in each cluster
orthoprop_cluster[,5:7] <- orthoprop_cluster[,5:7]/(orthoprop_cluster[,5]+orthoprop_cluster[,6]+orthoprop_cluster[,7])
# Reshape as cluster+orthoclass~counts
orthoprop_cluster <- recast(data = orthoprop_cluster, formula = ClusterID + sexDE + sexDI + DIDE + variable ~ .)
names(orthoprop_cluster)[5:6] <- c("Orthologs","Counts")
orthoprop_cluster <- droplevels(na.exclude(orthoprop_cluster))


multiprop_boxplot <- ggplot(data = orthoprop_cluster[which(orthoprop_cluster$Orthologs=="Multi"),], mapping = aes(x = DIDE, y = Counts, col = DIDE)) + geom_boxplot(varwidth = T) + theme_bw() + ylab(label = "Proportion of Nodes\nfrom Paralog Genes") + xlab(label = "") + theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "bottom") + scale_color_manual(values = cbPalette, name = "Sex-Bias Type") + geom_point(position = position_jitter(width = 0.1, height = 0)) + scale_y_continuous(limits = c(0,1)) + guides(col = guide_legend(nrow = 2))
multiprop_boxplot

pdf(file = file.path(graphdir, "Trans_Dupl_boxplot.pdf"), width = textWidthHeightmm[1]*.039, height = textWidthHeightmm[2]*.039/2, useDingbats = F)
multiplot(transprop_boxplot, multiprop_boxplot, cols = 2)
dev.off()

# Proportions from no OG genes
ggplot(data = orthoprop_cluster[which(orthoprop_cluster$Orthologs=="NoOG"),], mapping = aes(x = DIDE, y = Counts, col = DIDE)) + geom_boxplot(varwidth = T) + theme_bw() + ggtitle(label = "Proportions of Nodes from Orphan Genes in \nClusters with Sex Biased Expression and Correlation") + ylab(label = "Proportion of Nodes from Orphan Genes") + xlab(label = "") + theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_color_manual(values = cbPalette, name = "Sex-Bias Type")
# Proportion from single copy OG
ggplot(data = orthoprop_cluster[which(orthoprop_cluster$Orthologs=="Single"),], mapping = aes(x = DIDE, y = Counts, col = DIDE)) + geom_boxplot(varwidth = T) + theme_bw() + ggtitle(label = "Proportions of Nodes from Singlecopy Genes in \nClusters with Sex Biased Expression and Correlation") + ylab(label = "Proportion of Nodes from Singlecopy Genes") + xlab(label = "") + theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_color_manual(values = cbPalette, name = "Sex-Bias Type")

### Check and plot proportions of different strata
taxprop_cluster <- as.data.frame(table(transcriptdata$clusterID, transcriptdata$strata, useNA = "always"))
names(taxprop_cluster) <- c("ClusterID","Stratum","Counts")
taxprop_cluster$sexDE <- taxprop_cluster$ClusterID%in%sexDE_clusters
taxprop_cluster$sexDI <- taxprop_cluster$ClusterID%in%sexDI_clusters
taxprop_cluster$DIDE <- paste(taxprop_cluster$ClusterID%in%sexDI_clusters, taxprop_cluster$ClusterID%in%sexDE_clusters)
taxprop_cluster$DIDE <- factor(x = taxprop_cluster$DIDE, labels = c("No Bias", "DE Only", "DC Only", "DE and DC"))
taxprop_cluster <- recast(data = taxprop_cluster, formula = ClusterID + sexDE + sexDI + DIDE ~ Stratum)
names(taxprop_cluster)[10] <- c("None")
taxprop_cluster <- droplevels(na.exclude(taxprop_cluster))
taxprop_cluster$total <- apply(X = taxprop_cluster[,5:10], MARGIN = 1, FUN = sum)
taxprop_cluster <- taxprop_cluster[,c("ClusterID", "sexDE", "sexDI", "DIDE", "total", "None", "Wasp", "Hymenoptera", "Insect", "Arthropod", "Metazoa")]
# Store global proportions of strata across network
taxmeans <- apply(taxprop_cluster[,5:11], 2, sum)
taxmeans <- taxmeans[2:7]/taxmeans[1]
# Normalize by number of nodes per cluster
taxprop_cluster[6:11] <- taxprop_cluster[,6:11]/taxprop_cluster$total
taxprop_cluster <- taxprop_cluster[,-grep(pattern = "total", x = names(taxprop_cluster))]
# Normalize by the global proportions of strata in network
taxprop_cluster[,5:10] <- apply(X = taxprop_cluster[,5:10], MARGIN = 1, FUN = function(x){x/taxmeans})
taxprop_cluster <- recast(data = taxprop_cluster, formula = ClusterID + sexDE + sexDI + DIDE + variable ~ .)
names(taxprop_cluster)[5:6] <- c("Stratum", "Proportion")

# Remove class none, consistently rename other classes
taxprop_cluster <- droplevels(taxprop_cluster[which(taxprop_cluster$Stratum!="None"),])
taxprop_cluster$Stratum <- factor(x = taxprop_cluster$Stratum, labels = c("Nasonia", "Hymenoptera", "Insecta", "Arthropoda", "Metazoa"))

# Plot
ggplot(data = taxprop_cluster, mapping = aes(x = DIDE, y = Proportion, col = DIDE)) +
  geom_boxplot(varwidth = T) + 
  theme_bw() + 
  facet_grid(. ~ Stratum) + 
  ggtitle(label = "Taxonomic Depth of Nodes in Sex-Biased Clusters") + 
  ylab(label = "Excess Proportion of Nodes from Stratum") + 
  xlab(label = "") + 
  theme(
    panel.grid.major.x=element_blank(), 
    axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  scale_color_manual(values = cbPalette, name = "Sex-Bias Type") + 
  scale_linetype_discrete(name = "Differentially\nCorrelated") + 
  theme(legend.position = "bottom")
# + geom_point(position = position_jitter(width = 0.1, height = 0), alpha = "0.5")
ggsave(filename = file.path(graphdir, "PhyloStrata_boxplot.pdf"), device = "pdf", width = textWidthHeightmm[1], height = textWidthHeightmm[2]/2, units = "mm")

## Compare recombination rate distribution across clusters with different DEDI
transcriptdata <- read.csv(file = "./Output/Results_compiler_CCRE/transcriptdata_full.csv")[,-1]
recombplot <- merge(x = transcriptdata[,c("clusterID","geneID","recombrate","LinkageCluster")], y = DIDE_data[,c("clusterID","DIDE","sex_DE","sex_DI")], all.x = T, all.y = F)
# Set to one node per gene per cluster
recombplot <- unique(recombplot)
# Plot
ggplot(data = recombplot, mapping = aes(x = DIDE, y = recombrate, col = DIDE)) + geom_point(position = "jitter") + geom_boxplot(notch = T) + scale_y_log10() + theme_bw() + scale_color_manual(values = cbPalette, name = "Sex-Bias Type")
# Use LM to detect changes
library(lme4)
# remove entries with recombination rate zero or NA and grey genes
recombplot$recombrate <- ifelse(recombplot$recombrate==0,NA,recombplot$recombrate)
recombplot <- na.exclude(recombplot)
recombplot$sex_DE <- ifelse(recombplot$sex_DE=="Unbiased","Unbiased","DE")
recombplot$sex_DI <- ifelse(recombplot$sex_DI=="Unbiased","Unbiased","DI")
# glmer(formula = recombrate ~ sex_DI * sex_DE + (1|LinkageCluster), data = na.exclude(recombplot), family = Gamma)

glm1 <- glm(formula = recombrate ~ sex_DI * sex_DE , data = na.exclude(recombplot), family = Gamma(link = "log"))
plot(glm1)
step(glm1)
summary(glm1)
# No significant differences in RecombRate

### Compare cluster parameters of DIDE clusters (rework using general sexbias trajectories)
# networkplots <- unique(DIDE_data[,c("clusterID","DIDE")])
# networkplots <- networkplots[rev(order(networkplots$clusterID, networkplots$DIDE)),]
# networkplots <- networkplots[which(duplicated(networkplots$clusterID)==F),]
networkplots <- droplevels(unique(transprop_cluster[,c("ClusterID", "DIDE"),]))
names(networkplots) <- c("clusterID", "DIDE")
networkplots <- merge(networkplots, clusterdata[,1:11])

densibox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = Density, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) +   theme_bw() + theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type") + scale_y_continuous(limits = c(0,.15))

centrabox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = Centralization, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) +  theme_bw()  + theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type") + scale_y_continuous(limits = c(0,.15))

heterobox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = Heterogeneity, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) +  theme_bw() + theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type") + scale_y_continuous(limits = c(0,1.5))

genesbox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = nGenes, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + scale_y_log10(breaks = c(25,50,100,250,500)) + geom_point(position = position_jitter(width = 0.2)) + theme_bw()+ theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type")

transcriptbox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = nTranscripts, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + scale_y_log10(breaks = c(25,50,100,250,500)) + geom_point(position = position_jitter(width = 0.2)) + theme_bw()+ theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type")

clusterbox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = medianClusterCoef, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) +   theme_bw() + theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type")# + scale_y_continuous(limits = c(0,.15))

diambox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = diameter, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) +   theme_bw() + theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type") + scale_y_log10()

KMEbox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = medianKME+1, col =DIDE, shape = medianKME<0.5))+geom_boxplot(varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) +   theme_bw() + theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type") + scale_y_log10()

# ggplot(data = networkplots, mapping = aes(x = DIDE, y = nTranscripts, col =DIDE))+geom_boxplot() + scale_y_log10()

multiplot(densibox, heterobox, centrabox, transcriptbox, cols = 2)

# CP with networkwide associations
library(WGCNA)
library(corrgram)
z <- function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}
plot(hclust(as.dist(bicor(clusterdata[,2:11]))))


corrgram(clusterdata[,2:11], order=TRUE, lower.panel=panel.shade,upper.panel=panel.pie, text.panel=panel.txt)
corrgram(apply(X = clusterdata[,c(2:6,8:11)], 2, function(x){x/clusterdata$nNodes}), 
         order=TRUE, lower.panel=panel.shade,upper.panel=panel.pie, text.panel=panel.txt)
corrgram(apply(X = clusterdata[,c(2:6,8:11)], 2, function(x){x-clusterdata$nNodes}), 
         order=TRUE, lower.panel=panel.shade,upper.panel=panel.pie, text.panel=panel.txt)

corrgram(cbind(apply(X = clusterdata[,c(2:6,8:11)], 2, function(x){x/clusterdata$nNodes}), clusterdata$nNodes), 
         order=TRUE, lower.panel=panel.shade,upper.panel=panel.pie, text.panel=panel.txt)
corrgram(cbind(apply(X = clusterdata[,c(2:6,8:11)], 2, function(x){z(x)/z(clusterdata$nNodes)}), clusterdata$nNodes), 
         order=TRUE, lower.panel=panel.shade,upper.panel=panel.pie, text.panel=panel.txt)

plot(hclust(as.dist(bicor(transcriptdata[which(transcriptdata$Representative_Node==T|is.na(transcriptdata$Representative_Node)),c("recombrate", "KME_wc", "Connectivity_wc", "Betweenness_wc", "MAR_wc", "ClusterCoef_wc")], use = "pairwise"))))

corrgram(apply(X = transcriptdata[which(transcriptdata$Representative_Node==T|is.na(transcriptdata$Representative_Node)),c("recombrate", "KME_wc", "Connectivity_wc", "Betweenness_wc", "MAR_wc", "ClusterCoef_wc")], MARGIN = 2, FUN = function(x){log10(x+1)}), order=TRUE, lower.panel=panel.shade,upper.panel=panel.pie, text.panel=panel.txt)
corrgram(apply(X = transcriptdata[which(transcriptdata$Representative_Node==T|is.na(transcriptdata$Representative_Node)),c("recombrate", "RelKME_wc", "RelConnectivity_wc", "RelBetweenness_wc", "RelMAR_wc", "RelClusterCoef_wc")], MARGIN = 2, FUN = function(x){log10(x+1)}), order=TRUE, lower.panel=panel.shade,upper.panel=panel.pie, text.panel=panel.txt)
corrgram(apply(X = transcriptdata[which(transcriptdata$Representative_Node==T|is.na(transcriptdata$Representative_Node)),c("recombrate", "KME_wc", "AbsConnectivity_wc", "AbsBetweenness_wc", "MAR_wc", "ClusterCoef_wc")], MARGIN = 2, FUN = function(x){log10(x+1)}), order=TRUE, lower.panel=panel.shade,upper.panel=panel.pie, text.panel=panel.txt)


## Apply normalization to cluster size and then repeat all plots
networkplots_norm <- cbind(networkplots[,1:2], apply(networkplots[,3:12], 2, function(x){z(x)/z(networkplots$nNodes)}))


densibox <- ggplot(data = networkplots_norm, mapping = aes(x = DIDE, y = Density+5, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) +   theme_bw() + theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type") + scale_y_log10()# + scale_y_continuous(limits = c(0,.15))
densibox

centrabox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = Centralization+5, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) +  theme_bw()  + theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type") + scale_y_log10()# + scale_y_continuous(limits = c(0,.15))
centrabox

heterobox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = Heterogeneity+5, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) +  theme_bw() + theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type") + scale_y_log10() #scale_y_continuous(limits = c(0,1.5))
heterobox

genesbox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = nGenes, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + scale_y_log10(breaks = c(25,50,100,250,500)) + geom_point(position = position_jitter(width = 0.2)) + theme_bw()+ theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type")

transcriptbox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = nTranscripts, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + scale_y_log10(breaks = c(25,50,100,250,500)) + geom_point(position = position_jitter(width = 0.2)) + theme_bw()+ theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type")

clusterbox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = medianClusterCoef, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) +   theme_bw() + theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type")# + scale_y_continuous(limits = c(0,.15))

diambox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = diameter, col =DIDE))+geom_boxplot(varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) +   theme_bw() + theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type") + scale_y_log10()

KMEbox <- ggplot(data = networkplots, mapping = aes(x = DIDE, y = medianKME+1, col =DIDE, shape = medianKME<0.5))+geom_boxplot(varwidth = T, outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) +   theme_bw() + theme(legend.position = "bottom", panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_x_discrete(name = "Sex-Bias Type") + scale_y_log10()



## Compare proportions of nodes with transcription/splicing after removing nodes which are already transcriptionally sexbiased in the same stage

## Load DE data for nodes
glms2_fdr_coef <- read.csv(file = "./Output/sex_biased_nodes_CCRE/glms2_fdr_coef.csv")[,-1]
threshold <- 5E-2
## Reduce only to coefficients for nodes with significant lfdr (5E-2)
nodebias <- as.data.frame(apply(glms2_fdr_coef[,7:11], MARGIN = c(1,2), function(x){x<threshold})*apply(glms2_fdr_coef[,17:21], MARGIN = c(1,2), sign))
# Annotate node ID
nodebias$nodeID <- glms2_fdr_coef$node_ID
nodebias <- melt(nodebias)
# Annotate stage 
nodebias$stage <- factor(x = substr(nodebias$variable, start = 0, stop = 5), levels = c("emb10","emb18","lar51","pupye","adult"))


# # Filter only significant contrasts
# nodebias <- nodebias[which(nodebias$value!=0),]
# nodebias <- nodebias[,-which(names(nodebias)=="variable")]
# # Order constitutives first in node + stage ~ non-zero coef table
# nodebias <- recast(data = nodebias, formula = nodeID + stage ~ .)
# names(nodebias)[3] <- "coef"

# Merge nodes with transcript IDs (accounts for CCREs)
nodebias <- merge(nodebias, transcriptdata[,c("nodeID","eigenexonID","geneID")], all.x = T)
# Merge information on parent gene transcription bias
nodebias <- merge(nodebias, transcriptdata[,c("geneID","sexbiasedgene")], all.x = T)
# If the node is a transcription node, replace sexbiasedgene with transcription
nodebias$class <- ifelse(grepl(pattern = "_con", x = nodebias$eigenexonID), "Transcription", as.character(nodebias$sexbiasedgene))
nodebias$class <- factor(x = nodebias$class, levels = c("Transcription","Unbiased","Sexbiased"), labels = c("Transcription","Splicing no Trans", "Splicing and Trans"))
# nodebias <- nodebias[which(nodebias$value!=0),]
# Plot
nodebias2 <- data.frame(table(nodebias$stage, nodebias$class, sign(nodebias$value)))
names(nodebias2) <- c("stage","class","sexbias","counts")
nodebias2$grouper <- paste(nodebias2$class, nodebias2$sexbias)
ggplot(data = nodebias2, mapping = aes(x = stage, y = counts, lty = class, col = sexbias, group = grouper)) + geom_line() + geom_point()+ scale_y_log10()
# Simplify: remove stages and directions
nodebias3 <- ddply(.data = nodebias, .variables = .(eigenexonID), .fun = summarize, 
                   sexbiasedgene = sexbiasedgene[1],
                   value = ifelse(any(value!=0), "Sexbiased","Unbiased")
                   )
nodebias3$class <- ifelse(grepl(pattern = "_con", x = nodebias3$eigenexonID), "Transcription", as.character(nodebias3$sexbiasedgene))
nodebias3 <- data.frame(table(nodebias3$class, nodebias3$value))
names(nodebias3) <- c("class","sexbias","counts")
ggplot(data = nodebias3[which(nodebias3$sexbias!=0),], mapping = aes(x = sexbias, y = counts, fill = class)) + geom_line() + geom_bar(stat = "identity", position = "dodge")
ggplot(data = nodebias2[which(nodebias2$sexbias!=0),], mapping = aes(x = class, y = counts, fill = class)) + geom_line() + geom_bar(stat = "identity", position = "dodge") + theme_bw()

# Collapse at the GENE level
nodebias_con <- ddply(.data = nodebias[grep(pattern = "_con", x = nodebias$eigenexonID),], .variables = .(geneID), .fun = summarize, 
                   sexbiasedgene = ifelse(any(value!=0), "Sexbiased","Unbiased")
)
nodebias_fac <- ddply(.data = nodebias[grep(pattern = "_fac", x = nodebias$eigenexonID),], .variables = .(geneID), .fun = summarize, 
                      sexbiasedspl = ifelse(any(value!=0), "Sexbiased","Unbiased")
)
nodebias4 <- merge(nodebias_con, nodebias_fac, by = "geneID", all = T)
table(nodebias4$sexbiasedgene, nodebias4$sexbiasedspl, useNA = "ifany", dnn = c("sexbiasedgene", "sexbiasedspl"))
nodebias4$sexbiasedspl <- ifelse(is.na(nodebias4$sexbiasedspl), "Unspliced", nodebias4$sexbiasedspl)
nodebiascounts <- data.frame(table(nodebias4$sexbiasedgene, nodebias4$sexbiasedspl, useNA = "ifany", dnn = c("sexbiasedgene", "sexbiasedspl")))
nodebiascounts$grouper <- paste(nodebiascounts$sexbiasedgene, nodebiascounts$sexbiasedspl)

ggplot(data = nodebiascounts, mapping = aes(x = sexbiasedgene, y = Freq, fill = sexbiasedspl)) + geom_bar(stat = "identity", position = "dodge")

nodebias5 <- nodebias4
nodebias5$sexbiasedspl <- as.factor(ifelse(nodebias5$sexbiasedspl=="Sexbiased", "Sexbiased\nSplicing", "Unbiased\n Splicing"))
nodebias5$sexbiasedgene <- as.factor(ifelse(nodebias5$sexbiasedgene=="Sexbiased", "Sexbiased\nTranscription", "Unbiased\nTranscription"))

summary(nodebias5)
table(nodebias5$sexbiasedgene, nodebias5$sexbiasedspl)
prop.table(table(nodebias5$sexbiasedgene, nodebias5$sexbiasedspl))
prop.table(table(nodebias5$sexbiasedgene, nodebias5$sexbiasedspl), margin = 1)
prop.table(table(nodebias5$sexbiasedgene, nodebias5$sexbiasedspl), margin = 2)

mosaicplot(x = table(nodebias5$sexbiasedgene, nodebias5$sexbiasedspl), shade = T, main = "Proportions of Genes with \nSex-Biased Transcription and Splicing")

## Proportion plot before manual edits

pdf(file = file.path(graphdir, "TransVsSpl.pdf"), width = textWidthHeightmm[1]*.039/2, height = textWidthHeightmm[1]*.039/2, useDingbats = F)
mosaicplot(x = table(nodebias5$sexbiasedspl, nodebias5$sexbiasedgene), main = NULL, color = c(cbPalette[6:8]))
dev.off()

# # Retain only gene ID duplicates with the smallest index
# nodebias <- nodebias[order(nodebias$geneID, grepl(pattern = "fac", x = nodebias$eigenexonID)),]
# unique_nodebias <- nodebias[which(duplicated(nodebias[,c("geneID","stage")])==F),]
# # Annotate spl vs trans
# unique_nodebias$type <- factor(ifelse(grepl(pattern = "fac", x = unique_nodebias$eigenexonID), "Splicing", "Transcription"))
# # Plot tabulating by coef sign
# table(unique_nodebias$type, sign(unique_nodebias$coef))
# prop.table(table(unique_nodebias$type, sign(unique_nodebias$coef)))
# prop.table(table(unique_nodebias$type, sign(unique_nodebias$coef)), margin = 2)
# prop.table(table(unique_nodebias$type, sign(unique_nodebias$coef)), margin = 1)
# ## Male biased nodes are more than female-biased ones??
# table(unique_nodebias$stage, unique_nodebias$type, sign(unique_nodebias$coef))
# # Cp with whole dataset
# table(transcriptdata$event_type, grepl("[f,m]", transcriptdata$devsexbias))
# 
# duplicated_nodebias <- nodebias[which(duplicated(nodebias[,c("geneID","stage")])==T),]
# table(sign(duplicated_nodebias$coef))
# 
# # Print non-DE splicing nodes
# nonDEspl <- unique_nodebias[grep(pattern = "fac", x = unique_nodebias$eigenexonID), "eigenexonID"]
# write(x = nonDEspl, file = file.path(newdir, "nonDEspl.txt"))
# write.csv(x = unique_nodebias[grep(pattern = "fac", x = unique_nodebias$eigenexonID), c("nodeID","eigenexonID")], file = file.path(newdir, "IndepSplicing.csv"))


# Plot number of unbiased/DE/DSDE/DS genes vs dev time vs sex