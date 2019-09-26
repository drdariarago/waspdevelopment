#### Summaries from compiled cluster data
# Version without dOut
## Compile results from the different analyses
# Clean workspace
rm(list=ls())
# Initialize script
library(reshape2)
library(plyr)
library(stringr)
library(vcd)
library(ggplot2)
library(lattice)
library(fdrtool)
# library(MuMIn)
# options(na.action="na.fail")
newdir<-file.path(getwd(),"Output/Results_compiler_summaries")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/Results_compiler_summaries")
dir.create(graphdir)
source(file = "./Scripts/multiplot.R")

# Load data
moduledata_short<-read.csv(file = "./Output/Results_compiler/moduledata_short.csv")


################## Begin summaries
# How many nodes are present in sex-biased clusters in total?
sum(moduledata_short[which(moduledata_short$signif==T), "Nodes"])
sum(moduledata_short[which(moduledata_short$signif==T), "Nodes"])/sum(moduledata_short$Nodes)
# How many nodes are present in sex-biased clusters in total?
sum(moduledata_short[which(moduledata_short$signif==T), "Genes_n"])
sum(moduledata_short[which(moduledata_short$signif==T), "Genes_n"])/sum(moduledata_short$Genes_n)

# Plot size distribution of clusters (in genes and nodes)
pdf(file = file.path(graphdir, "Node_distribution.pdf"), height = 3.153, width = 6.072)
ggplot(moduledata_short, aes(x=Nodes))+geom_dotplot()+scale_x_log10(breaks = c(1,5,10,20,50,100,250,500,1000,2000,5000,10000), limits=c(20,10000))+theme_bw()+ggtitle(label = "Size of Transcriptional Clusters in Nodes\n")+facet_grid(signif~.)
ggplot(moduledata_short, aes(x=Genes_n))+geom_dotplot()+scale_x_log10(breaks = c(1,5,10,20,50,100,250,500,1000,2000,5000,10000), limits=c(20,10000))+theme_bw()+ggtitle(label = "Size of Transcriptional Clusters in Genes\n")+facet_grid(signif~.)
dev.off()

# Joint distribution
ggplot(moduledata_short, aes(x=Nodes/Genes_n, fill=signif))+geom_dotplot()+theme_bw()+facet_grid(signif~.)+scale_x_log10(breaks=seq(from = 1, to = 1.9, by = 0.1))
ggplot(moduledata_short, aes(y=Nodes/Genes_n, x=signif))+geom_boxplot(notch = T, varwidth = T)+theme_bw()+scale_y_log10()+scale_y_log10(breaks=seq(from = 1, to = 1.9, by = 0.1))

# Compare size of sexb vs non
ggplot(moduledata_short, aes(x=signif, y=Nodes))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+scale_y_log10()
ggplot(moduledata_short, aes(x=signif, y=Genes_n))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+scale_y_log10()

# Plot number of clusters dI in each stage
plot_dWithin_clusters<-ggplot(moduledata_short, aes(x=stagedW, y=as.numeric(signif), fill=as.factor(sexdW)))+theme_bw()+stat_summary(fun.y = "sum", geom = "bar", position="dodge")+scale_y_continuous(name="Clusters", limits = c(0,10))+scale_x_continuous(name="Stage", breaks = 1:5, labels = c("embryo\n10", "embryo\n18", "larva\n51", "pupa", "adult"))+theme(legend.position="bottom", legend.title=element_blank())+scale_fill_discrete(labels=c("Female biased", "Male biased"))+ggtitle(label = "Number of clusters \n with Sex-Specific bias in Integration\n")

dWithin_nodecounts<-na.exclude(moduledata_short[,c("Nodes", "stagedW", "sexdW")])
plot_dWithintNodes<-ggplot(dWithin_nodecounts, aes(x=stagedW, y=Nodes, fill=as.factor(sexdW)))+theme_bw()+stat_summary(fun.y = "sum", geom = "bar", position="dodge")+ylim(c(0,4000))+ggtitle(label = "Number of nodes present in clusters \n with Sex-Specific bias in Integration\n")+theme(legend.position="bottom", legend.title=element_blank())+scale_fill_discrete(labels=c("Female biased", "Male biased"))+scale_x_continuous(name="Stage", breaks = c(1:5), labels = c("embryo\n10", "embryo\n18", "larva\n51", "pupa", "adult"))

# # Compare proportions of paralogs in sexb vs non (includes a NA cathegory, consider restricting only to genes with defined ID)
# ggplot(moduledata_short, aes(x=signif, y=Paralogs_OGS_p))+geom_boxplot(varwidth = T, notch = T)+theme_bw()
# ggplot(moduledata_short, aes(x=signif, y=Orthologs_OGS_p))+geom_boxplot(varwidth = T, notch = T)+theme_bw()

# Compare proportions of paralogs in sexb vs non (excluding non-classified genes)
ggplot(moduledata_short, aes(x=signif, y=Paralogs_OGS_p2))+geom_boxplot(varwidth = T, notch = T)+theme_bw()

# Compare proportions of isoforms in sexb vs non
ggplot(moduledata_short, aes(x=signif, y=Splicing_p))+geom_boxplot(varwidth = T, notch = T)+theme_bw()
ggplot(moduledata_short, aes(x=signif, y=Trans_p))+geom_boxplot(varwidth = T, notch = T)+theme_bw()

## Plot for paper
para<-ggplot(moduledata_short, aes(x=signif, y=Paralogs_OGS_p2))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+scale_x_discrete(labels=c("Unbiased \n Clusters", "Sex-Biased \n Clusters"))+theme(axis.title.x = element_blank())+scale_y_continuous(name = "Proportion of Nodes in each Cluster\n belonging to Duplicated Genes", limits = c(0,1))

iso<-ggplot(moduledata_short, aes(x=signif, y=Splicing_p))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+scale_x_discrete(labels=c("Unbiased \n Clusters", "Sex-Biased \n Clusters"))+theme(axis.title.x = element_blank())+scale_y_continuous(name = "Proportion of Splicing Nodes in each Cluster", limits = c(0,1))

pdf(file = file.path(graphdir, "para_iso_boxplots.pdf"), width = 5.5, height = 5.5, useDingbats = F)
multiplot(para, iso, cols = 2)
dev.off()

# Compare joint distribution of proportion of paralogs vs isoforms in signif vs non signif clusters
ggplot(moduledata_short, aes(x=Paralogs_OGS_p2, y=Splicing_p, col=signif))+geom_point()+theme_bw()+geom_density2d()+geom_rug()+facet_wrap(~signif)

pdf(file = file.path(graphdir, "Joint_para_iso_sexbias.pdf"), width = 5.5, height = 5.5, useDingbats = F)
ggplot(moduledata_short, aes(x=Paralogs_OGS_p2, y=Splicing_p, col=signif))+geom_point(alpha="0.5")+theme_bw()+geom_density2d()+ggtitle(label = "Proportions of Splicing Nodes and Nodes from Duplicated Genes\n in Clusters with and without Significant Sex-Bias\n")+scale_x_continuous(limits=c(0,1), name = "Proportion of Nodes in each Cluster belonging to Duplicated Genes")+scale_y_continuous(limits=c(0,1), name = "Proportion of Splicing Nodes in each Cluster")+theme(legend.title = element_blank(), legend.position = "bottom")+scale_color_brewer(type = "qual", palette = 7, labels = c("Non Sex-Biased","Sex-Biased"))
dev.off()


# Compare proportions of adult meth in sexb vs non
ggplot(moduledata_short, aes(x=signif, y=Meth_p2))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+geom_point(position="jitter")
ggplot(moduledata_short, aes(x=signif, y=Meth_p2))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+geom_point(position="jitter")+scale_y_log10()
ggplot(moduledata_short, aes(fill=signif, x=Meth_p2))+geom_histogram(position="dodge")+theme_bw()+scale_x_log10()
ggplot(moduledata_short, aes(col=signif, x=Meth_p2))+geom_density()+theme_bw()
#########
# Check evolutionary distances of clusters
ggplot(moduledata_short, aes(col=signif, x=medianNresiduals))+geom_density()+theme_bw()
ggplot(moduledata_short, aes(col=signif, x=medianRatio))+geom_density()+theme_bw()

ggplot(moduledata_short, aes(x=Paralogs_OGS_p2, y=medianNresiduals, col=signif))+geom_point(alpha="1")+theme_bw()+geom_rug(alpha="0.5")+ggtitle(label = "Proportions of Multicopy Genes and Evolutionary Rates \n in Clusters with and without significant Sex-Bias\n")+scale_x_continuous(name = "Proportion of Multicopy Gene Nodes")+scale_y_continuous(name = "Wasp Evolutionary rate")+theme(legend.title = element_blank(), legend.position = "bottom")+scale_color_brewer(type = "qual", palette = 7, labels = c("Non Sex-Biased","Sex-Biased"))+geom_smooth(method="lm")
