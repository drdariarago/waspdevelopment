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
library(pscl)
library(MuMIn)
options(na.action="na.fail")
newdir<-file.path(getwd(),"Output/Results_compiler")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/Results_compiler")
dir.create(graphdir)
source(file = "./Scripts/multiplot.R")

## Load OGS2 genome
NVIT_OGS2_goodannotcomplete <- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")

## Load eigenexon assignments
eigenexons_assignments <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexons_assignments.csv")

# ## Load eigenexon E values
eigenexon_evalues <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexon_evalues.csv")

## Load cluster assignments
clusterassignments_full <- read.csv("./Output/WGCNA_clustering_biweight/clusterassignments_full.csv")
names(clusterassignments_full)[1]<-"eigenexonID"

# ## Load Taxonomic profile (from Sackton, 2013)
tax_depth <- read.delim("./Input/DatasetS1.tsv")
tax_depth <- tax_depth[,c("ogs2", "strata")]
names(tax_depth)<-c("geneID","strata")
tax_depth$strata<-factor(tax_depth$strata, levels = c("Metazoa", "Arthropod", "Insect", "Hymenoptera", "Wasp"))

#### Merge datasets
splicing_cat<-unique(eigenexons_assignments[,c("geneID","splicing_cathegory","clusters")])
splicing_cat<-ddply(splicing_cat, .(geneID), summarize, splicing_cathegory=splicing_cathegory[1], clusters=max(clusters))

NVIT_OGS2<-merge(splicing_cat, NVIT_OGS2_goodannotcomplete, by="geneID", all.x=F)
# gene Nasvi2EG021272 isn't present in the genome annotation file, most likely removed because of low quality

NVIT_OGS2<-merge(NVIT_OGS2, tax_depth, by = "geneID")
  
# NVIT_OGS2<-merge(NVIT_OGS2, clusterassignments_full[,c("eigenexonID","clusterID")])
# NVIT_OGS2<-merge(NVIT_OGS2, OGS2_Taxascan, by.x="geneID", by.y="Gene_ID", all.x=T)

clusterassignments_short<-clusterassignments_full[,c("eigenexonID", "clusterID")]
clusterassignments_short$geneID<-str_extract(string = clusterassignments_short$eigenexonID, pattern = "^[^_]*")

clusterdata<-merge(clusterassignments_short, NVIT_OGS2[,c("geneID","ODB6_OG_ID","quality7","IntMatch","IntPred","OG_Copynumber","ratio","Nresiduals","splicing_cathegory","clusters","adult_female_meth_status","strata")], by="geneID")

# Save as .csv
write.csv(clusterdata, file=file.path(newdir, "clusterdata_compiled.csv"))

# Count number of clusters per gene (assesses diversification between isoforms)

cluster_isoform_count<-ddply(clusterdata, .(geneID), summarize, 
                             nClusters = length(unique(clusterID)), 
                             nIsoforms = length(unique(eigenexonID)),
                             .progress = "text")

densityplot(cluster_isoform_count[which(cluster_isoform_count$nIsoforms>1),]$nClusters/cluster_isoform_count[which(cluster_isoform_count$nIsoforms>1),]$nIsoforms)
table(cluster_isoform_count$nClusters/cluster_isoform_count$nIsoforms<1)
prop.table(table(cluster_isoform_count$nClusters/cluster_isoform_count$nIsoforms<1))

ggplot(cluster_isoform_count, aes(x=nClusters/nIsoforms))+geom_density()+theme_bw()+scale_x_log10()

## Summarize by module
moduledata<-ddply(clusterdata, .(clusterID), summarize, .progress = "text",
                  Nodes=length(unique(eigenexonID)),
                  Genes=length(unique(geneID)),
                  Splicing=sum(grepl("con", eigenexonID)),
                  Trans=sum(grepl("fac", eigenexonID)),
                  Paralogs_OGS=sum(grepl("Paralog",quality7)),
                  Orthologs_OGS=sum(grepl("Ortholog",quality7)),
                  Paralogs_ODB=sum(is.na(OG_Copynumber)==T),
                  Orthologs_ODB=sum(is.na(OG_Copynumber)==F),
                  Meth=sum(adult_female_meth_status=="Methylated", na.rm = T),
                  Unmeth=sum(adult_female_meth_status=="Unmethylated", na.rm = T),
                  medianRatio=median(ratio, na.rm = T),
                  medianNresiduals=median(Nresiduals, na.rm = T),
                  taxaNasonia=table(strata)[5],
                  taxaHymenoptera=table(strata)[4],
                  taxaInsecta=table(strata)[3],
                  taxaArthropoda=table(strata)[2],
                  taxaMetazoa=table(strata)[1]
)

# Normalize by number of nodes in cluster
prop_moduledata<-cbind(clusterID = moduledata[,1],moduledata[,-c(1,2,12,13)]/moduledata[,2])
moduledata<-merge(moduledata, prop_moduledata, by="clusterID", suffixes = c("_n","_p"))

# Normalize OGS orthologs/paralogs by number of nodes with annotation in cluster
moduledata$Paralogs_OGS_p2<-moduledata$Paralogs_OGS_n/(moduledata$Paralogs_OGS_n+moduledata$Orthologs_OGS_n)

# Normalize methylation status by number of nodes with annotation in cluster
moduledata$Meth_p2<-moduledata$Meth_n/(moduledata$Meth_n+moduledata$Unmeth_n)

## Load stage-specific sexbias coefficients
fdr_perm_dWithin <- read.csv("./Output/sex_modularity_GLMs/fdr_perm_dWithin.csv")
names(fdr_perm_dWithin)[2:3]<-c("Factornames","clusterID")
fdr_perm_dOut <- read.csv("./Output/sex_modularity_GLMs/fdr_perm_dOut.csv")
names(fdr_perm_dOut)[2:3]<-c("Factornames","clusterID")
fdr_perm_data<-merge(fdr_perm_dOut[,c(2,3,4,9,10)], fdr_perm_dWithin[,c(2,3,4,9,10)], by=c("Factornames","clusterID"), suffixes = c("_dOut","_dWithin"))

# Merge with annotation data
moduledata<-merge(moduledata, fdr_perm_data, by = "clusterID", all=T)

# Save compiled dataset
write.csv(moduledata, file=file.path(newdir, "moduledata.csv"))

# How many clusters have sex-biased changes in either parameter?
length(unique(fdr_perm_data[fdr_perm_data$lFDR_adjusted_dOut<0.1|fdr_perm_data$lFDR_adjusted_dWithin<0.05,"clusterID"]))
length(unique(fdr_perm_data[fdr_perm_data$lFDR_adjusted_dOut<0.1|fdr_perm_data$lFDR_adjusted_dWithin<0.05,"clusterID"]))/length(levels(fdr_perm_data$clusterID))

## Create data.frame with stages with significant contrasts for each cluster
# significant_contrasts<-ddply(fdr_perm_data, .(clusterID), summarize, stage = Factornames[which(lFDR_adjusted_dOut<0.05|lFDR_adjusted_dWithin<0.05)], .progress = "text")

significant_contrasts<-ddply(fdr_perm_data, .(clusterID), summarize, 
                             stagedO = ifelse(any(lFDR_adjusted_dOut<0.1),Factornames[which(lFDR_adjusted_dOut<0.1)],NA),
                             sexdO = ifelse(any(lFDR_adjusted_dOut<0.1),sign(Estimate_dOut[which(lFDR_adjusted_dOut<0.1)]),NA), 
                             stagedW = ifelse(any(lFDR_adjusted_dWithin<0.05),Factornames[which(lFDR_adjusted_dWithin<0.05)],NA),
                             sexdW = ifelse(any(lFDR_adjusted_dWithin<0.05),sign(Estimate_dWithin[which(lFDR_adjusted_dWithin<0.05)]),NA), 
                             .progress = "text")

# Select only one row for each cluster, dividing between clusters with at least one significant stage:sex interaction and no sexbias
moduledata_short<-ddply(moduledata, .(clusterID), summarize, signif = any(lFDR_adjusted_dOut<0.1|lFDR_adjusted_dWithin<0.05), .progress = "text")
moduledata_short<-unique(merge(moduledata_short, moduledata[,c("clusterID","Nodes","Genes_n","Splicing_p","Trans_p","Paralogs_OGS_p2","Paralogs_OGS_p","Orthologs_OGS_p","Meth_p","Unmeth_p","Meth_p2","medianRatio","medianNresiduals","taxaNasonia_p","taxaHymenoptera_p","taxaInsecta_p","taxaArthropoda_p","taxaMetazoa_p","taxaNasonia_n","taxaHymenoptera_n","taxaInsecta_n","taxaArthropoda_n","taxaMetazoa_n")], by = "clusterID", all.x=F, all.y=F))

# Add stage and direction of significance
moduledata_short<-merge(moduledata_short, significant_contrasts, by="clusterID")

# How many nodes are present in sex-biased clusters in total?
sum(moduledata_short[which(moduledata_short$signif==T), "Nodes"])
# How many nodes are present in sex-biased clusters in total?
sum(moduledata_short[which(moduledata_short$signif==T), "Genes_n"])


# Plot size distribution of clusters (in genes and nodes)
pdf(file = file.path(graphdir, "Node_distribution.pdf"), height = 3.153, width = 6.072)
ggplot(moduledata_short, aes(x=Nodes))+geom_histogram()+scale_x_log10(breaks = c(1,5,10,20,50,100,250,500,1000,2000,5000,10000), limits=c(20,10000))+theme_bw()+ggtitle(label = "Size of Transcriptional Clusters in Nodes\n")
ggplot(moduledata_short, aes(x=Nodes))+geom_dotplot()+scale_x_log10(breaks = c(1,5,10,20,50,100,250,500,1000,2000,5000,10000), limits=c(20,10000))+theme_bw()+ggtitle(label = "Size of Transcriptional Clusters in Nodes\n")+facet_grid(.~signif)
dev.off()

ggplot(moduledata_short, aes(x=Genes_n))+geom_dotplot()+scale_x_log10()+theme_bw()+facet_grid(.~signif)

ggplot(moduledata_short, aes(x=Nodes, y=Genes_n))+geom_density2d()+geom_point()+scale_x_log10()+scale_y_log10()+theme_bw()+geom_rug()
# ggplot(moduledata_short, aes(x=Nodes, y=Splicing_p))+geom_density2d()+geom_point()+scale_x_log10()+scale_y_log10()+theme_bw()+geom_rug()+geom_smooth()

# Compare size of sexb vs non
ggplot(moduledata_short, aes(x=signif, y=Nodes))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+scale_y_log10()
ggplot(moduledata_short, aes(x=signif, y=Genes_n))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+scale_y_log10()

# Plot number of clusters dI/dC in each stage
# ggplot(moduledata_short, aes(x=stagedO, fill=as.factor(sexdO)))+geom_histogram(position="dodge")+theme_bw()
plot_dOut_clusters<-ggplot(moduledata_short, aes(x=stagedO, y=as.numeric(signif), fill=as.factor(sexdO)))+theme_bw()+stat_summary(fun.y = "sum", geom = "bar", position="dodge")+scale_y_continuous(name="Clusters", limits = c(0,30))+scale_x_continuous(name="Stage", breaks = 1:5, labels = c("embryo\n10", "embryo\n18", "larva\n51", "pupa", "adult"))+theme(legend.position="bottom", legend.title=element_blank())+scale_fill_discrete(labels=c("Female biased", "Male biased"))+ggtitle(label = "Number of clusters \n with Sex-Specific bias in Constraint\n")

# ggplot(moduledata_short, aes(x=stagedW, fill=as.factor(sexdW)))+geom_histogram(position="dodge")+theme_bw()
plot_dWithin_clusters<-ggplot(moduledata_short, aes(x=stagedW, y=as.numeric(signif), fill=as.factor(sexdW)))+theme_bw()+stat_summary(fun.y = "sum", geom = "bar", position="dodge")+scale_y_continuous(name="Clusters", limits = c(0,30))+scale_x_continuous(name="Stage", breaks = 1:5, labels = c("embryo\n10", "embryo\n18", "larva\n51", "pupa", "adult"))+theme(legend.position="bottom", legend.title=element_blank())+scale_fill_discrete(labels=c("Female biased", "Male biased"))+ggtitle(label = "Number of clusters \n with Sex-Specific bias in Integration\n")

# Plot number of nodes dI/dC in each stage
source('./Scripts/multiplot.R', echo=F)

dOut_nodecounts<-na.exclude(moduledata_short[,c("Nodes", "stagedO", "sexdO")])

plot_dOutNodes<-ggplot(dOut_nodecounts, aes(x=stagedO, y=Nodes, fill=as.factor(sexdO)))+theme_bw()+stat_summary(fun.y = "sum", geom = "bar", position="dodge")+ylim(c(0,4000))+ggtitle(label = "Number of nodes present in clusters \n with Sex-Specific bias in Constraint\n")+theme(legend.position="bottom", legend.title=element_blank())+scale_fill_discrete(labels=c("Female biased", "Male biased"))+scale_x_continuous(name="Stage", breaks = 1:5, labels = c("embryo\n10", "embryo\n18", "larva\n51", "pupa", "adult"))

dWithin_nodecounts<-na.exclude(moduledata_short[,c("Nodes", "stagedW", "sexdW")])
plot_dWithintNodes<-ggplot(dWithin_nodecounts, aes(x=stagedW, y=Nodes, fill=as.factor(sexdW)))+theme_bw()+stat_summary(fun.y = "sum", geom = "bar", position="dodge")+ylim(c(0,4000))+ggtitle(label = "Number of nodes present in clusters \n with Sex-Specific bias in Integration\n")+theme(legend.position="bottom", legend.title=element_blank())+scale_fill_discrete(labels=c("Female biased", "Male biased"))+scale_x_continuous(name="Stage", breaks = c(1:5), labels = c("embryo\n10", "embryo\n18", "larva\n51", "pupa", "adult"))

pdf(file = file.path(graphdir, "cluster_dynamics.pdf"), width = 8.605, height = 5.378)
multiplot(plot_dOutNodes, plot_dOut_clusters, plot_dWithintNodes, plot_dWithin_clusters, cols =2)
dev.off()

# Compare proportions of paralogs in sexb vs non (includes a NA cathegory, consider restricting only to genes with defined ID)
ggplot(moduledata_short, aes(x=signif, y=Paralogs_OGS_p))+geom_boxplot(varwidth = T, notch = T)+theme_bw()
ggplot(moduledata_short, aes(x=signif, y=Orthologs_OGS_p))+geom_boxplot(varwidth = T, notch = T)+theme_bw()

# Compare proportions of paralogs in sexb vs non (excluding non-classified genes)
ggplot(moduledata_short, aes(x=signif, y=Paralogs_OGS_p2))+geom_boxplot(varwidth = T, notch = T)+theme_bw()

# Compare proportions of isoforms in sexb vs non
ggplot(moduledata_short, aes(x=signif, y=Splicing_p))+geom_boxplot(varwidth = T, notch = T)+theme_bw()
ggplot(moduledata_short, aes(x=signif, y=Trans_p))+geom_boxplot(varwidth = T, notch = T)+theme_bw()

## Plot for SMBE poster
para<-ggplot(moduledata_short, aes(x=signif, y=Paralogs_OGS_p))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+scale_x_discrete(labels=c("Unbiased \n Clusters", "Sex-Biased \n Clusters"))+theme(axis.title.x = element_blank())+scale_y_continuous(name = "Proportion of Nodes in each Cluster\n belonging to Duplicated Genes", limits = c(0,1))

iso<-ggplot(moduledata_short, aes(x=signif, y=Splicing_p))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+scale_x_discrete(labels=c("Unbiased \n Clusters", "Sex-Biased \n Clusters"))+theme(axis.title.x = element_blank())+scale_y_continuous(name = "Proportion of Splicing Nodes in each Cluster", limits = c(0,1))

pdf(file = file.path(graphdir, "para_iso_boxplots.pdf"), width = 5.5, height = 5.5, useDingbats = F)
multiplot(para, iso, cols = 2)
dev.off()

# Compare joint distribution of proportion of paralogs vs isoforms in signif vs non signif clusters
ggplot(moduledata_short, aes(x=Paralogs_OGS_p2, y=Splicing_p, col=signif))+geom_point()+theme_bw()+geom_density2d()+geom_rug()+facet_wrap(~signif)

pdf(file = file.path(graphdir, "Joint_para_iso_sexbias.pdf"), width = 5.5, height = 5.5, useDingbats = F)
ggplot(moduledata_short, aes(x=Paralogs_OGS_p2, y=Splicing_p, col=signif))+geom_point(alpha="0.5")+theme_bw()+geom_density2d()+ggtitle(label = "Proportions of Splicing Nodes and Nodes from Duplicated Genes\n in Clusters with and without Significant Sex-Bias\n")+scale_x_continuous(name = "Proportion of Nodes in each Cluster belonging to Duplicated Genes")+scale_y_continuous(name = "Proportion of Splicing Nodes in each Cluster")+theme(legend.title = element_blank(), legend.position = "bottom")+scale_color_brewer(type = "qual", palette = 7, labels = c("Non Sex-Biased","Sex-Biased"))
dev.off()

# ggplot(moduledata_short, aes(x=Paralogs_OGS_p2+1, y=Splicing_p+1, col=signif))+geom_point()+theme_bw()+geom_density2d()+geom_rug()+facet_wrap(~signif)+scale_x_log10()+scale_y_log10()

# Compare proportions of adult meth in sexb vs non
ggplot(moduledata_short, aes(x=signif, y=Meth_p2))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+geom_point(position="jitter")
ggplot(moduledata_short, aes(x=signif, y=Meth_p2))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+geom_point(position="jitter")+scale_y_log10()
ggplot(moduledata_short, aes(fill=signif, x=Meth_p2))+geom_histogram(position="dodge")+theme_bw()+scale_x_log10()
ggplot(moduledata_short, aes(col=signif, x=Meth_p2))+geom_density()+theme_bw()

# ggplot(moduledata_short, aes(x=signif, y=Meth_p))+geom_boxplot(varwidth = T, notch = T)+theme_bw()

# Check evolutionary distances of clusters
ggplot(moduledata_short, aes(col=signif, x=medianNresiduals))+geom_density()+theme_bw()
ggplot(moduledata_short, aes(col=signif, x=medianRatio))+geom_density()+theme_bw()

ggplot(moduledata_short, aes(x=Paralogs_OGS_p2, y=medianNresiduals, col=signif))+geom_point(alpha="0.5")+theme_bw()+geom_density2d()+geom_rug(alpha="0.5")+ggtitle(label = "Proportions of Multicopy Genes and Splicing Nodes \n in Clusters with and without significant Sex-Bias\n")+scale_x_continuous(name = "Proportion of Multicopy Gene Nodes")+scale_y_continuous(name = "Proportion of Splicing Nodes")+theme(legend.title = element_blank(), legend.position = "bottom")+scale_color_brewer(type = "qual", palette = 7, labels = c("Non Sex-Biased","Sex-Biased"))


# Compare taxonomic depth of clusters

# Compare proprtion of nodes for each stratum, split by sexbias
taxdepth_boxplot<-recast(moduledata_short[,c("clusterID","signif","taxaNasonia_p","taxaHymenoptera_p","taxaInsecta_p","taxaArthropoda_p","taxaMetazoa_p")], formula = clusterID+signif+variable~.)
names(taxdepth_boxplot)<-c("clusterID","signif","stratum","proportion")
ggplot(taxdepth_boxplot, aes(x=signif, y=proportion))+geom_boxplot(notch = T, varwidth = T)+facet_grid(.~stratum)+theme_bw()

# Compare proportion of sexbiased vs non-sexbiased nodes in each stratum
taxdepth_barplot<-recast(moduledata_short[,c("clusterID","signif","taxaNasonia_n","taxaHymenoptera_n","taxaInsecta_n","taxaArthropoda_n","taxaMetazoa_n")], formula = signif+variable~., fun.aggregate=sum)
names(taxdepth_barplot)<-c("signif","stratum","number")
# taxdepth_barplot<-ddply(taxdepth_barplot, .variables = .(signif, stratum), .fun = summarize, value=sum(value))
ggplot(taxdepth_barplot, aes(x=stratum, col=signif, y=number))+geom_point()+theme_bw()
taxdepth_barplot<-recast(taxdepth_barplot, formula = stratum~signif)
taxdepth_barplot<-rbind(taxdepth_barplot, c(NA,sum(taxdepth_barplot$"FALSE"),sum(taxdepth_barplot$"TRUE")))
taxdepth_barplot$stratum<-as.factor(c(as.character(taxdepth_barplot$stratum)[-6],"allnodes_n"))
taxdepth_barplot$prop_sexb<-taxdepth_barplot$"TRUE"/(taxdepth_barplot$"FALSE"+taxdepth_barplot$"TRUE")
taxdepth_barplot$prop_unb<-taxdepth_barplot$"FALSE"/(taxdepth_barplot$"FALSE"+taxdepth_barplot$"TRUE")
taxdepth_barplot<-melt(taxdepth_barplot[,c("stratum","prop_sexb","prop_unb")], variable.name = "sexbias")
ggplot(taxdepth_barplot, aes(x=stratum, col=variable, y=value))+geom_point()+geom_line(aes(group=variable))+theme_bw()
ggplot(taxdepth_barplot, aes(x=stratum, fill=variable, y=value))+geom_bar(stat="identity")+theme_bw()


taxdepth_mosaicplot<-recast(moduledata_short[,c("clusterID","signif","taxaNasonia_n","taxaHymenoptera_n","taxaInsecta_n","taxaArthropoda_n","taxaMetazoa_n")], 
                            formula = variable~signif, fun.aggregate=sum)
row.names(taxdepth_mosaicplot)<-taxdepth_mosaicplot$variable
taxdepth_mosaicplot<-taxdepth_mosaicplot[,-1]
mosaicplot(x = taxdepth_mosaicplot, shade = T)
######################## Outdated code snippets


# Recode as cluster~factor table
# dWithin_avg<-recast(data = dWithin_avg, formula = clusterID~factorName)
# dOut_avg<-recast(data = dOut_avg, formula = clusterID~factorName)

# 
# ## Merge datasets
# eigenexon_data<-merge(clusterassignments_full[,c("eigenexonID", "clusterID"),],eigenexons_assignments[c("eigenexonID","geneID","constitutive")],all.x=T, all.y=T)
# eigenexon_data<-merge(eigenexon_data, dOut_avg, by="clusterID", all.x=T, all.y=F)
# eigenexon_data<-merge(eigenexon_data, dWithin_avg, by="clusterID", all.x=T, all.y=F, suffixes = c(".dOut",".dWithin"))
# eigenexon_data<-merge(eigenexon_data, OGS2_Taxascan, by.x="geneID", by.y="Gene_ID", all.x=T, all.y=F)
# eigenexon_data$Taxascan_max<-ifelse(is.na(eigenexon_data$Taxascan_max),0,eigenexon_data$Taxascan_max)