### Summary statistics on sexbiased nodes
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
textWidthHeightmm <- c(152.73, 216.00)
source(file = "./Scripts/multiplot.R")
# The colorblind friendly palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

graphdir<-file.path(getwd(), "Graphics/NodeStats")
dir.create(graphdir)

# ## Load DE data for nodes
# nodebias <- read.csv(file = "./Output/sex_biased_nodes_CCRE/glms2_fdr_coef.csv")[,-1]
# ## Reduce only to coefficients for nodes with significant lfdr (5E-2)
# # Create coef table
# nodebias_2 <- nodebias[,12:21]
# row.names(nodebias_2) <- nodebias$node_ID
# # Create significance matrix
# nodebias_sign <- nodebias[,2:11]
# nodebias_sign <- apply(X = as.matrix(nodebias_sign), MARGIN = c(1,2), function(x){ifelse(x<=5E-2, 1, 0)})
# # Multiply coefficients by significance (reduces NS to zero)
# nodebias <- nodebias_2*nodebias_sign
# 
# # Plot distribution of node bias vs dev stage
# nodebias_stage <- nodebias[,1:5]
# nodebias_stage$type <- "stage"
# names(nodebias_stage) <- str_extract(string = names(nodebias_stage), pattern = "^[^_]*")
# 
# nodebias_sex <- nodebias[,6:10]
# nodebias_sex$type <- "sex"
# names(nodebias_sex) <- str_extract(string = names(nodebias_sex), pattern = "^[^_]*")
# 
# nodebias_split <- rbind(nodebias_sex, nodebias_stage)
# nodebias_split$nodeID <- row.names(nodebias_split)
# nodebias_split <- recast(data = nodebias_split, formula = type+nodeID+variable~.)
# names(nodebias_split) <- c("type","nodeID","stage","E")
# nodebias_split$sign <- factor(sign(nodebias_split$E), labels = c("Female","Unbiased","Male"))
# 
# 
# pdf(file = file.path(graphdir, "Sexbiased_Nodecounts.pdf"))
# ggplot(data = nodebias_split[which(nodebias_split$E!=0&nodebias_split$type=="sex"),], mapping = aes(x = stage, fill = sign))+geom_bar(position = "dodge")+theme_bw()+ggtitle(label = "Sex-Biased Nodes per Stage")
# dev.off()

## Number of biased splicing nodes which do not also have sexbiased transcription at the same stage


## Same, with number of clusters
clusterGLMsummaries <- read.csv(file = "./Output/DEcluster_analysis_CCRE/clusterGLMsummaries.csv")[,-1]
clusterGLMsummaries$stage <- factor(x = clusterGLMsummaries$stage, levels = c("emb10","emb18","lar51","pupyel","adult"))
clusterGLMsummaries$sex <- sign(clusterGLMsummaries$coeff*(clusterGLMsummaries$fdr<0.05))
clusterGLMsummaries$sex <- factor(ifelse(clusterGLMsummaries$contrast=="Sex", clusterGLMsummaries$sex, NA), labels = c("Female","Unbiased","Male"))
clusterGLMsummaries$clusterID <- substring(text = clusterGLMsummaries$clusterID, first = 3)

pdf(file = file.path(graphdir,"clusterBiasvsStage.pdf"))
ggplot(data = clusterGLMsummaries[which(clusterGLMsummaries$contrast=="Sex"&clusterGLMsummaries$sex!="Unbiased"),], mapping = aes(x = stage, fill = sex))+geom_bar(position = "dodge")+theme_bw()+ggtitle(label = "Number of Sex-Biased Clusters across Development")+scale_fill_brewer(type = "qual", palette =  2, name = "Sex Bias")
dev.off()

## Add DI coefficients and plot in series of 2D graphs
# import and format data as cluster+stage~coeff+fdr
DI_data <- read.csv("./Output/sex_modularity_GLMs_CCRE/fdr_perm_dWithin.csv")[,-1]
DI_data <- DI_data[,c("coefWithinFactornames",".id","Estimate","lFDR_adjusted")]
DI_data$stage <- substring(text = str_extract(string = DI_data$coefWithinFactornames, pattern = "^[^:]*"), first = 6)
DI_data$stage <- factor(x = DI_data$stage, levels = c("emb10","emb18","lar51","pupyel","adult"))
DI_data <- DI_data[,c(".id","stage","Estimate","lFDR_adjusted")]
names(DI_data) <- c("clusterID","stage","coeff","fdr")

## Calculate significance of DIDC
clusterdata <- read.csv(
  file = "./Output/Results_compiler_CCRE/clusterdata_full.csv", 
  header = T, row.names = 1)
fisher.test(table(clusterdata$cluster_devsexbias!=".....", is.na(clusterdata$DiffIntegrated)==F, dnn = c("DE","DI")))
mosaicplot(table(clusterdata$cluster_devsexbias!=".....", is.na(clusterdata$DiffIntegrated)==F, dnn = c("DE","DI")), shade = T, main = "Convergence between\n Differential Expression and \nDifferential Correlation")
table(clusterdata$cluster_devsexbias, clusterdata$DiffIntegrated, useNA = "ifany")
# See if there is any sign matching
mosaicplot(table(grepl(pattern = "[m]", x = clusterdata$cluster_devsexbias), grepl(pattern = "male", x = clusterdata$DiffIntegrated)), shade = T)
mosaicplot(table(grepl(pattern = "[f]", x = clusterdata$cluster_devsexbias), grepl(pattern = "Female", x = clusterdata$DiffIntegrated)), shade = T)

# merge with reduced clusterDE data
DIDE_data <- merge(x = clusterGLMsummaries[which(clusterGLMsummaries$contrast=="Sex"),c("clusterID","stage","fdr","coeff")], 
                   y = DI_data, 
                   by = c("clusterID","stage"),
                   suffixes = c("_DE","_DI"))

DIDE_data$sex_DE <- factor(sign((DIDE_data$fdr_DE<0.05)*DIDE_data$coeff_DE), labels = c("Female","Unbiased","Male"))
DIDE_data$sex_DI <- factor(sign((DIDE_data$fdr_DI<0.05)*DIDE_data$coeff_DI), labels = c("Female","Unbiased","Male"))
DIDE_data$minfdr <- pmin(DIDE_data$fdr_DE, DIDE_data$fdr_DI)
DIDE_data$DIDE <- paste(DIDE_data$fdr_DI<0.05, DIDE_data$fdr_DE<0.05)
DIDE_data$DIDE <- factor(x = DIDE_data$DIDE, labels = c("No Bias", "DE Only", "DC Only", "DE and DC"))

DIDE_data <- merge(DIDE_data, clusterdata[,c("clusterID","nTranscripts")])

# ggplot(data = DIDE_data, mapping = aes(x = coeff_DE, y = coeff_DI, col = minfdr<.05, size = nTranscripts))+
#   geom_point()+facet_wrap(~stage)+theme_bw()+scale_color_brewer(type = "seq", palette = 1, name = "fdr")

library(ggrepel)
myClusters <- c("navajowhite3", "lavenderblush3", "palevioletred2", "darkseagreen2", "antiquewhite4", "plum", "plum3", "thistle3", "antiquewhite2")
DIDE_data_2 <- DIDE_data
levels(DIDE_data_2$stage) <- c('Embryo, Early', 'Embryo, Late', 'Larva', 'Pupa', 'Adult') 

# pdf(file = file.path(graphdir, "effectsizes_DIDE.pdf"))
ggplot(data = DIDE_data_2, mapping = aes(x = coeff_DE, y = coeff_DI, col = DIDE, alpha = DIDE=="No Bias", label = clusterID)) +
  geom_point() + 
  facet_wrap(~stage) + 
  geom_text_repel(data = subset(DIDE_data_2, clusterID%in%myClusters & DIDE!="No Bias" & coeff_DI > 0.08), 
    box.padding = 0.2, point.padding = 0.3, nudge_y = 0.15, force = 50, color='black') +
  geom_text_repel(data = subset(DIDE_data_2, clusterID%in%myClusters & DIDE!="No Bias" & coeff_DI < 0.08), 
    box.padding = 0.2, point.padding = 0.3, nudge_y = -0.12, nudge_x = 0.4, force = 2, color='black', min.segment.length = 0.02) +
  theme_bw() + 
  scale_color_manual(values = cbPalette, name = "Sex-Bias Type") + 
  scale_alpha_manual(values = c(.6,.1), guide = F) + 
  scale_size(name = "Cluster Size\nin Transcripts", breaks = c(20,50,100,250,500)) + 
  scale_x_continuous(name = "Differential Expression Coefficient") + 
  scale_y_continuous(name = "Differential Correlation Coefficient") + 
  theme(legend.position = c(.8,.2), legend.key.size = unit(3, "mm")) + 
  ggtitle(label = "Differential Expression and Differential Correlation Coefficients within Stages") + 
  guides(colour = guide_legend(nrow = 2, override.aes=list(size=5)), size = guide_legend(nrow = 2)) 
# dev.off()
ggsave(filename = file.path(graphdir, "effectsizes_DIDE.pdf"), device = "pdf", width = 200, height = 150, units = 'mm') 
ggsave(filename = file.path(graphdir, "effectsizes_DIDE.png"), device = "png", width = 200, height = 150, units = 'mm') 
  #width = textWidthHeightmm[1], height = textWidthHeightmm[2]/2, units = "mm")

## q-value 2d plot
ggplot(data = DIDE_data, mapping = aes(x = fdr_DE, y = fdr_DI, shape = sex_DI, col = sex_DE))+geom_point()+facet_wrap(~stage)+theme_bw()+scale_x_log10()+scale_y_log10()+geom_hline(yintercept=.05)+geom_vline(xintercept=.05)

## Histogram of cluster counts across stages for DIDE
ggplot(data = DIDE_data[which(DIDE_data$DIDE!="No Bias"),], mapping = aes(x = stage, fill = DIDE))+geom_bar(position = "dodge")+theme_bw()+ggtitle(label = "Number of Sex-Biased Clusters across Development")+scale_fill_manual(values = cbPalette[-1], name = "Sex Bias")


# Check proportion of transcription vs splicing and dupl vs single-copy nodes in sexDE clusters
transcriptdata <- read.csv(file = "./Output/Results_compiler_CCRE/transcriptdata_full.csv", row.names = 1)
clusterdata <- read.csv(file = "./Output/Results_compiler_CCRE/clusterdata_full.csv", row.names = 1)

sexDE_clusters <- clusterdata$clusterID[grep(pattern = "[m,f]", clusterdata$cluster_devsexbias)]
sexDI_clusters <- clusterdata$clusterID[which(x = is.na(clusterdata$DiffIntegrated)==F)]
transprop_cluster <- as.data.frame(table(transcriptdata$clusterID, transcriptdata$event_type))
names(transprop_cluster) <- c("ClusterID","EventType","Counts")
transprop_cluster$sexDE <- transprop_cluster$ClusterID%in%sexDE_clusters
transprop_cluster$sexDI <- transprop_cluster$ClusterID%in%sexDI_clusters
transprop_cluster$DEDI <- paste(transprop_cluster$ClusterID%in%sexDE_clusters, transprop_cluster$ClusterID%in%sexDI_clusters)
transprop_cluster <- recast(data = transprop_cluster, formula = ClusterID + sexDE + sexDI + DEDI ~ EventType)

ggplot(data = transprop_cluster, mapping = aes(x = sexDE, y = splicing/(splicing+transcription))) + geom_boxplot(notch = T, varwidth = T)
ggplot(data = transprop_cluster, mapping = aes(x = sexDI, y = splicing/(splicing+transcription))) + geom_boxplot(notch = T, varwidth = T)

# pdf(file = file.path(graphdir, "SplProp_cluster.pdf"))
transprop_boxplot <- ggplot(data = transprop_cluster, mapping = aes(x = DEDI, y = splicing/(splicing+transcription), group = DEDI, linetype = sexDI, col = sexDE)) + geom_boxplot(varwidth = T) + theme_bw() + ggtitle(label = "Proportions of Splicing Nodes in Clusters with\n Sex Biased Expression and Correlation") + ylab(label = "Proportion of Splicing Nodes") + xlab(label = "") + theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_color_brewer(type = "qual", name = "Differentially\nExpressed") + scale_linetype_discrete(name = "Differentially\nCorrelated")
# dev.off()  

## Control transcription/splicing proportions via GLM
transGLM_binom <- glm(formula = cbind(splicing, transcription) ~ sexDE*sexDI, family = binomial, data = transprop_cluster)
plot(transGLM_binom)
step(transGLM_binom)
summary(transGLM_binom)

duplprop_cluster <- as.data.frame(table(transcriptdata$clusterID, transcriptdata$quality7))
names(duplprop_cluster) <- c("ClusterID","Orthology","Counts")
duplprop_cluster$sexDE <- duplprop_cluster$ClusterID%in%sexDE_clusters
duplprop_cluster$sexDI <- duplprop_cluster$ClusterID%in%sexDI_clusters
duplprop_cluster$DEDI <- paste(duplprop_cluster$ClusterID%in%sexDE_clusters, duplprop_cluster$ClusterID%in%sexDI_clusters)
duplprop_cluster <- recast(data = duplprop_cluster, formula = ClusterID + sexDE + sexDI + DEDI ~ Orthology)

ggplot(data = duplprop_cluster, mapping = aes(x = sexDE, y = Paralog/(Ortholog+Paralog+None))) + geom_boxplot(notch = T, varwidth = T)
ggplot(data = duplprop_cluster, mapping = aes(x = sexDI, y = Paralog/(Ortholog+Paralog+None))) + geom_boxplot(notch = T, varwidth = T)

ggplot(data = duplprop_cluster, mapping = aes(x = sexDE, y = Ortholog/(Ortholog+Paralog+None))) + geom_boxplot(notch = T, varwidth = T)
ggplot(data = duplprop_cluster, mapping = aes(x = sexDI, y = Ortholog/(Ortholog+Paralog+None))) + geom_boxplot(notch = T, varwidth = T)

ggplot(data = duplprop_cluster, mapping = aes(x = sexDE, y = None/(Ortholog+Paralog+None))) + geom_boxplot(notch = T, varwidth = T)
ggplot(data = duplprop_cluster, mapping = aes(x = sexDI, y = None/(Ortholog+Paralog+None))) + geom_boxplot(notch = T, varwidth = T)

ggplot(data = duplprop_cluster, mapping = aes(x = DEDI, y = Paralog/(Ortholog+Paralog+None))) + geom_boxplot(notch = T, varwidth = T) + theme_bw()

## Check proportions of orthologs vs paralogs using OGS family counts

orthoprop_cluster <- as.data.frame(table(transcriptdata$clusterID, transcriptdata$Copynumber>1, useNA = "always"))
names(orthoprop_cluster) <- c("ClusterID","Orthologs","Counts")
orthoprop_cluster$sexDE <- orthoprop_cluster$ClusterID%in%sexDE_clusters
orthoprop_cluster$sexDI <- orthoprop_cluster$ClusterID%in%sexDI_clusters
orthoprop_cluster$DEDI <- paste(orthoprop_cluster$ClusterID%in%sexDE_clusters, orthoprop_cluster$ClusterID%in%sexDI_clusters)
orthoprop_cluster <- recast(data = orthoprop_cluster, formula = ClusterID + sexDE + sexDI + DEDI ~ Orthologs)
names(orthoprop_cluster)[5:7] <- c("Single","Multi","NoOG")
orthoprop_cluster <- droplevels(na.exclude(orthoprop_cluster))

pairs(log10(orthoprop_cluster[,c("Single","Multi","NoOG")]))

ggplot(data = orthoprop_cluster, mapping = aes(x = sexDE, y = Multi/(Single+Multi+NoOG))) + geom_boxplot(notch = T, varwidth = T)
ggplot(data = orthoprop_cluster, mapping = aes(x = sexDI, y = Multi/(Single+Multi+NoOG))) + geom_boxplot(notch = T, varwidth = T)

ggplot(data = orthoprop_cluster, mapping = aes(x = sexDE, y = NoOG/(Single+Multi+NoOG))) + geom_boxplot(notch = T, varwidth = T)
ggplot(data = orthoprop_cluster, mapping = aes(x = sexDI, y = NoOG/(Single+Multi+NoOG))) + geom_boxplot(notch = T, varwidth = T)

ggplot(data = orthoprop_cluster, mapping = aes(x = sexDE, y = Single/(Single+Multi+NoOG))) + geom_boxplot(notch = T, varwidth = T)
ggplot(data = orthoprop_cluster, mapping = aes(x = sexDI, y = Single/(Single+Multi+NoOG))) + geom_boxplot(notch = T, varwidth = T)

# Use GLM to disentangle individual factors
paraGLMbinom <- glm(formula = cbind(Multi,Single+NoOG) ~ sexDE*sexDI, family = binomial, data = orthoprop_cluster)
plot(paraGLMbinom)
step(paraGLMbinom)
summary(paraGLMbinom)

orthoGLMbinom <- glm(formula = cbind(Single,Multi+NoOG) ~ sexDE*sexDI, family = binomial, data = orthoprop_cluster)
plot(orthoGLMbinom)
step(orthoGLMbinom)
summary(orthoGLMbinom)

NoOGGLMbinom <- glm(formula = cbind(NoOG,Multi+Single) ~ sexDE*sexDI, family = binomial, data = orthoprop_cluster)
plot(NoOGGLMbinom)
step(NoOGGLMbinom)
summary(NoOGGLMbinom)

orthoprop_cluster[,5:7] <- orthoprop_cluster[,5:7]/(orthoprop_cluster[,5]+orthoprop_cluster[,6]+orthoprop_cluster[,7])

pairs(log10(orthoprop_cluster[,c("Single","Multi","NoOG")]))
pairs(orthoprop_cluster[,c("Single","Multi","NoOG")])

cor.test(orthoprop_cluster$Multi, orthoprop_cluster$NoOG)
cor.test(log10(orthoprop_cluster$Multi+1), log10(orthoprop_cluster$NoOG+1))
cor.test(orthoprop_cluster$Multi, orthoprop_cluster$Single)
cor.test(orthoprop_cluster$Single, orthoprop_cluster$NoOG)

orthoprop_cluster <- recast(data = orthoprop_cluster, formula = ClusterID + sexDE + sexDI + DEDI + variable ~ .)
names(orthoprop_cluster)[5:6] <- c("Orthologs","Counts")
orthoprop_cluster <- droplevels(na.exclude(orthoprop_cluster))

ggplot(data = orthoprop_cluster, mapping = aes(x = Orthologs, y = Counts, col = sexDE)) + geom_boxplot(notch = T) + theme_bw()
ggplot(data = orthoprop_cluster, mapping = aes(x = Orthologs, y = Counts, col = sexDI)) + geom_boxplot(notch = T) + theme_bw()
ggplot(data = orthoprop_cluster, mapping = aes(x = Orthologs, y = Counts, col = DEDI, linetype = sexDE)) + geom_boxplot(notch = T) + theme_bw()

# pdf(file = file.path(graphdir, "MultProp_cluster.pdf"))
multiprop_boxplot <- ggplot(data = orthoprop_cluster[which(orthoprop_cluster$Orthologs=="Multi"),], mapping = aes(x = DEDI, y = Counts, col = DEDI)) + geom_boxplot(varwidth = T) + theme_bw() + ggtitle(label = "Proportions of Nodes from Paralog Genes in \nClusters with Sex Biased Expression and Correlation") + ylab(label = "Proportion of Nodes from Paralog Genes") + xlab(label = "") + theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_color_manual(values = cbPalette, name = "Differentially\nExpressed")
# dev.off()  
multiplot(transprop_boxplot, multiprop_boxplot)

ggplot(data = orthoprop_cluster[which(orthoprop_cluster$Orthologs=="NoOG"),], mapping = aes(x = DEDI, y = Counts, col = sexDE, linetype = sexDI)) + geom_boxplot(varwidth = T) + theme_bw() + ggtitle(label = "Proportions of Nodes from Orphan Genes in \nClusters with Sex Biased Expression and Correlation") + ylab(label = "Proportion of Nodes from Orphan Genes") + xlab(label = "") + theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_color_brewer(type = "qual", name = "Differentially\nExpressed") + scale_linetype_discrete(name = "Differentially\nCorrelated")

ggplot(data = orthoprop_cluster[which(orthoprop_cluster$Orthologs=="Single"),], mapping = aes(x = DEDI, y = Counts, col = sexDE, linetype = sexDI)) + geom_boxplot(varwidth = T) + theme_bw() + ggtitle(label = "Proportions of Nodes from Single Copy Genes in \nClusters with Sex Biased Expression and Correlation") + ylab(label = "Proportion of Nodes from Single Copy Genes") + xlab(label = "") + theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_color_brewer(type = "qual", name = "Differentially\nExpressed") + scale_linetype_discrete(name = "Differentially\nCorrelated")

# Check and plot proportions of different strata
taxprop_cluster <- as.data.frame(table(transcriptdata$clusterID, transcriptdata$strata, useNA = "always"))
names(taxprop_cluster) <- c("ClusterID","Stratum","Counts")
taxprop_cluster$sexDE <- taxprop_cluster$ClusterID%in%sexDE_clusters
taxprop_cluster$sexDI <- taxprop_cluster$ClusterID%in%sexDI_clusters
taxprop_cluster$DEDI <- paste(taxprop_cluster$ClusterID%in%sexDE_clusters, taxprop_cluster$ClusterID%in%sexDI_clusters)
taxprop_cluster <- recast(data = taxprop_cluster, formula = ClusterID + sexDE + sexDI + DEDI ~ Stratum)
names(taxprop_cluster)[10] <- c("None")
taxprop_cluster <- droplevels(na.exclude(taxprop_cluster))
taxprop_cluster$total <- apply(X = taxprop_cluster[,5:10], MARGIN = 1, FUN = sum)
taxprop_cluster <- taxprop_cluster[,c("ClusterID", "sexDE", "sexDI", "DEDI", "total", "None", "Wasp", "Hymenoptera", "Insect", "Arthropod", "Metazoa")]

pairs(taxprop_cluster[6:11]/taxprop_cluster$total)
plot(hclust(as.dist(cor(taxprop_cluster[,6:11]))))

# Fit one glm per stratum, apply fdr and check if correlated with DEDI
WaspGLM <- glm(formula = cbind(Wasp,total-Wasp) ~ sexDE*sexDI, family = binomial, data = taxprop_cluster)
HymGLM <- glm(formula = cbind(Hymenoptera,total-Hymenoptera) ~ sexDE*sexDI, family = binomial, data = taxprop_cluster)
InsGLM <- glm(formula = cbind(Insect,total-Insect) ~ sexDE*sexDI, family = binomial, data = taxprop_cluster)
ArtGLM <- glm(formula = cbind(Arthropod,total-Arthropod) ~ sexDE*sexDI, family = binomial, data = taxprop_cluster)
MetGLM <- glm(formula = cbind(Metazoa,total-Metazoa) ~ sexDE*sexDI, family = binomial, data = taxprop_cluster)
NonGLM <- glm(formula = cbind(None,total-None) ~ sexDE*sexDI, family = binomial, data = taxprop_cluster)

taxprop_cluster_2 <- taxprop_cluster
taxprop_cluster_2[6:11] <- taxprop_cluster_2[,6:11]/taxprop_cluster_2$total
taxprop_cluster_2 <- taxprop_cluster_2[,-grep(pattern = "total", x = names(taxprop_cluster_2))]
taxprop_cluster_2 <- recast(data = taxprop_cluster_2, formula = ClusterID + sexDE + sexDI + DEDI + variable ~ .)
names(taxprop_cluster_2)[5:6] <- c("Stratum", "Proportion")

ggplot(data = taxprop_cluster_2, mapping = aes(x = DEDI, y = Proportion, col = sexDE, linetype = sexDI)) + geom_boxplot(notch = T) + scale_y_log10() + theme_bw() + facet_grid(. ~ Stratum)
ggplot(data = taxprop_cluster_2, mapping = aes(x = DEDI, y = Proportion, col = sexDE, linetype = sexDI)) + geom_boxplot(notch = T) + theme_bw() + facet_grid(. ~ Stratum)

# Normalize by overall observed proportions
taxmeans <- apply(taxprop_cluster[,5:11], 2, sum)
taxmeans <- taxmeans[2:7]/taxmeans[1]

taxprop_cluster_2 <- taxprop_cluster
taxprop_cluster_2[,6:11] <- taxprop_cluster_2[,6:11]/taxprop_cluster_2$total
taxprop_cluster_2 <- taxprop_cluster_2[,-grep(pattern = "total", x = names(taxprop_cluster_2))]
taxprop_cluster_2[,5:10] <- apply(X = taxprop_cluster_2[,5:10], MARGIN = 1, FUN = function(x){x/taxmeans})

taxprop_cluster_2 <- recast(data = taxprop_cluster_2, formula = ClusterID + sexDE + sexDI + DEDI + variable ~ .)
names(taxprop_cluster_2)[5:6] <- c("Stratum", "Proportion")

ggplot(data = taxprop_cluster_2, mapping = aes(x = Stratum, y = Proportion, col = sexDE, linetype = sexDI, group = DEDI)) + geom_smooth()
ggplot(data = taxprop_cluster_2, mapping = aes(x = DEDI, y = Proportion, col = sexDE, linetype = sexDI)) + geom_boxplot(notch = T) + theme_bw() + facet_grid(. ~ Stratum)

pdf(file = file.path(graphdir, "TaxDepth_Cluster.pdf"))
ggplot(data = taxprop_cluster_2, mapping = aes(x = DEDI, y = Proportion, col = sexDE, linetype = sexDI)) + geom_boxplot(varwidth = T) + theme_bw() + facet_grid(. ~ Stratum) + ggtitle(label = "Taxonomic Stratum Enrichment in Clusters with \nDifferential Expression and Differential Correlation") + ylab(label = "Excess Proportion of Nodes from Stratum") + xlab(label = "") + theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_color_brewer(type = "qual", name = "Differentially\nExpressed") + scale_linetype_discrete(name = "Differentially\nCorrelated")
dev.off()

# Strong significant enrichment of Wasp stratum genes for DEDI, even distribution of DE, strong depletion of DI in metazoan

## Which stage DE patterns are more frequent in male and femalebiased clusters?
# x axis stageDE, y axis counts, col = sexbias (m,f)
sexDE <- ifelse(grepl(pattern = "f", clusterdata$cluster_devsexbias), ifelse(grepl(pattern = "m", clusterdata$cluster_devsexbias), "Both", "Female"), ifelse(grepl(pattern = "m", clusterdata$cluster_devsexbias), "Male", "Unbiased"))
sexDE <- factor(x = sexDE, levels = c("Unbiased", "Female", "Male", "Both"))
sexDE <- data.frame(clusterID = clusterdata$clusterID, sexbias = sexDE)
sexDE <- sexDE[which(sexDE$sexbias!="Both"),]

devbiascounts <- clusterGLMsummaries[which(clusterGLMsummaries$contrast=="Stage"&clusterGLMsummaries$fdr<.05),]
devbiascounts <- merge(devbiascounts, sexDE)

ggplot(data = devbiascounts, mapping = aes(x = stage, fill = sexbias)) + geom_bar(stat = "count", position = "dodge") + scale_y_log10() + theme_bw()

pdf(file = file.path(graphdir, "SexBias_DevCounts.pdf"))
ggplot(data = devbiascounts, mapping = aes(x = stage, col = sexbias, group = sexbias)) + geom_smooth(stat = "count") + theme_bw() + ggtitle(label = "Developmental Differential Expression of Clusters\n with Male and Female Differential Expression") + scale_y_continuous(name = "Number of Differentially Expressed Clusters") + scale_x_discrete(name = "Stage") + scale_color_manual(values = cbPalette, name = "Sex Bias")
dev.off()
# Most female biased clusters are differentially expressed in pupal/adult stages
# Most male biased clusters are expressed in early embryos

## How many nodes are transcriptionally sexbiased if the main transcript isn't already sexbiased?
# Three proportions unbiased splbiased transbiased tra+spl biased
# Unit of study: transcript
# create 2 matrices splicing x stage and gene x stage (with duplicated rows)
# assing logical sexDE to splicing (1 is DE) and logical sexDE to transcription (0 is DE)
# Multiply (retains only spl with no DE trans)


ggplot(data = nodedata, mapping = aes(y = RelKME_wc, x = clusterID, col = clusterID=="grey"))+geom_boxplot()
table(clusterdata$cluster_devsexbias)
## Cluster navajowhite3 has female-bias in emb10!
# Cluster characterization: high density, small size, low heterogeneity