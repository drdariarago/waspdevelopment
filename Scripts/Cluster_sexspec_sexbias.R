### Calculate overrepresentation in sex-biased and sex-specific genes
# Clean workspace
rm(list=setdiff(ls(),"ODB8_EukOGs_genes"))
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
newdir<-file.path(getwd(),"Output/Cluster_sexspec_sexbias")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/Cluster_sexspec_sexbias")
dir.create(graphdir)
source(file = "./Scripts/multiplot.R")

# Load dataset

moduledata_short<-read.csv(file = "./Output/Results_compiler/moduledata_for_enrichment.csv")[,-1]

## Compare proportions of sex-specific nodes in clusters with and without dI
# ggplot(moduledata_short, aes(y=sexspec_p, x=as.factor(signif)))+geom_boxplot(notch = T, varwidth = T)+geom_point(position="jitter")+theme_bw()
# ggplot(moduledata_short[which(moduledata_short$sexspec_p>0),], aes(y=sexspec_p, x=as.factor(signif)))+geom_boxplot(notch = T, varwidth = T)+geom_point(position="jitter")+theme_bw()
# 
# ## Compare proportions of sexbiased nodes in clusters with and without dI
# ggplot(moduledata_short, aes(y=sexbias_p, x=as.factor(signif)))+geom_boxplot(notch = T, varwidth = T)+geom_point(position="jitter")+theme_bw()+scale_y_log10()
# ggplot(moduledata_short[which(moduledata_short$sexbias_p>0),], aes(y=sexbias_p, x=as.factor(signif)))+geom_boxplot(notch = T, varwidth = T)+geom_point(position="jitter")+theme_bw()+scale_y_log10()
# 
# # Compare joint distribution split by cluster dI
# ggplot(moduledata_short, aes(x=sexbias_p, y=sexspec_p, col=signif))+geom_point()+theme_bw()+scale_x_log10()+scale_y_log10()+geom_rug()

# Perform enrichment test for sxb nodes on clusters
d_sxbenrichment<-moduledata_short[,c("sexbias_n", "Nodes")]
module<-1:nrow(d_sxbenrichment)
d_sxbenrichment<-rbind(d_sxbenrichment, apply(d_sxbenrichment, 2, sum))

sxbenrichment_pval<-sapply(X = module, FUN = function(x){
  fisher.test(
    rbind(d_sxbenrichment[x,],d_sxbenrichment[nrow(d_sxbenrichment),])
  )$p.value
})
sxbenrichment_pval<-sxbenrichment_pval*max(module) # Bonferroni correction

d_sxbenrichment$prop<-d_sxbenrichment$sexbias_n/d_sxbenrichment$Nodes
d_sxbenrichment$enriched<-sign(d_sxbenrichment$prop-d_sxbenrichment$prop[nrow(d_sxbenrichment)])
d_sxbenrichment$enriched<-ifelse(sign(d_sxbenrichment$enriched)>0, "Enriched", "Depleted")

d_sxbenrichment<-d_sxbenrichment[-nrow(d_sxbenrichment),]
d_sxbenrichment$pval<-sxbenrichment_pval

d_sxbenrichment$enriched<-ifelse(d_sxbenrichment$pval<0.01, d_sxbenrichment$enriched, "Not Significant")
d_sxbenrichment<-data.frame(clusterID = moduledata_short$clusterID, d_sxbenrichment)

# Order by category and save as csv
write.csv(d_sxbenrichment[order(d_sxbenrichment$enriched, d_sxbenrichment$pval),], file = "./Output/sxb_enriched_clusters_e-1.csv")


# Perform enrichment test for sxsp nodes on clusters
d_sxspenrichment<-moduledata_short[,c("sexspec_n", "Nodes")]
module<-1:nrow(d_sxspenrichment)
d_sxspenrichment<-rbind(d_sxspenrichment, apply(d_sxspenrichment, 2, sum))

sxspenrichment_pval<-sapply(X = module, FUN = function(x){
  fisher.test(
    rbind(d_sxspenrichment[x,],d_sxspenrichment[nrow(d_sxspenrichment),])
  )$p.value
})
sxspenrichment_pval<-sxspenrichment_pval*max(module) # Bonferroni correction

d_sxspenrichment$prop<-d_sxspenrichment$sexspec_n/d_sxspenrichment$Nodes
d_sxspenrichment$enriched<-sign(d_sxspenrichment$prop-d_sxspenrichment$prop[nrow(d_sxspenrichment)])
d_sxspenrichment$enriched<-ifelse(sign(d_sxspenrichment$enriched)>0, "Enriched", "Depleted")

d_sxspenrichment<-d_sxspenrichment[-nrow(d_sxspenrichment),]
d_sxspenrichment$pval<-sxspenrichment_pval

d_sxspenrichment$enriched<-ifelse(d_sxspenrichment$pval<0.01, d_sxspenrichment$enriched, "Not Significant")
d_sxspenrichment<-data.frame(clusterID=moduledata_short$clusterID, d_sxspenrichment)

# Order by category and save as csv
write.csv(d_sxspenrichment[order(d_sxspenrichment$enriched, d_sxspenrichment$pval),], file = "./Output/sxsp_enriched_clusters_e-1.csv")

# ## Merge and compare
# 
# d_sxsp_sxb<-merge(d_sxbenrichment, d_sxspenrichment, by=c("clusterID","Nodes"), suffixes = c("_sxb","_sxsp"))
# table(d_sxsp_sxb$enriched_sxb, d_sxsp_sxb$enriched_sxsp)
# with(droplevels(d_sxsp_sxb[which(d_sxsp_sxb$enriched_sxb!="Not Significant"|d_sxsp_sxb$enriched_sxsp!="Not Significant"),]),
#      mosaic(table(enriched_sxb, enriched_sxsp), shade = T))
# with(droplevels(d_sxsp_sxb[which(d_sxsp_sxb$enriched_sxb!="Not Significant"&d_sxsp_sxb$enriched_sxsp!="Not Significant"),]),
#      mosaic(table(enriched_sxb, enriched_sxsp), shade = T))
# 
# # Include into module data, and check for distribution of nodes in each
# moduledata_short<-merge(moduledata_short, d_sxsp_sxb[,c("clusterID","enriched_sxb","enriched_sxsp")], by="clusterID")
# 
# table(moduledata_short$signif, moduledata_short$enriched_sxb)
# mosaic(table(moduledata_short$signif, moduledata_short$enriched_sxb), shade = T)
# mosaic(table(moduledata_short[which(moduledata_short$enriched_sxb!="Not Significant"),]$signif, moduledata_short[which(moduledata_short$enriched_sxb!="Not Significant"),]$enriched_sxb), shade = T)
# fisher.test(table(moduledata_short[which(moduledata_short$enriched_sxb!="Not Significant"),]$signif, moduledata_short[which(moduledata_short$enriched_sxb!="Not Significant"),]$enriched_sxb))
# 
# table(moduledata_short$signif, moduledata_short$enriched_sxsp)
# mosaic(table(moduledata_short[which(moduledata_short$enriched_sxsp!="Not Significant"),]$signif, moduledata_short[which(moduledata_short$enriched_sxsp!="Not Significant"),]$enriched_sxsp), shade = T)
# fisher.test(table(moduledata_short[which(moduledata_short$enriched_sxsp!="Not Significant"),]$signif, moduledata_short[which(moduledata_short$enriched_sxsp!="Not Significant"),]$enriched_sxsp))

### Cluster enrichment checks
threshold<-0.01

## Check for clusters enriched in male-specific genes

d_malespecenrichment<-moduledata_short[,c("malespec_n", "Nodes")]
module<-1:nrow(d_malespecenrichment)
d_malespecenrichment<-rbind(d_malespecenrichment, apply(d_malespecenrichment, 2, sum))

malespecenrichment_pval<-sapply(X = module, FUN = function(x){
  fisher.test(
    rbind(d_malespecenrichment[x,],d_malespecenrichment[nrow(d_malespecenrichment),])
  )$p.value
})
malespecenrichment_pval<-malespecenrichment_pval*max(module) # Bonferroni correction

d_malespecenrichment$prop<-d_malespecenrichment$malespec_n/d_malespecenrichment$Nodes
d_malespecenrichment$enriched<-sign(d_malespecenrichment$prop-d_malespecenrichment$prop[nrow(d_malespecenrichment)])
d_malespecenrichment$enriched<-ifelse(sign(d_malespecenrichment$enriched)>0, "Enriched", "Depleted")

d_malespecenrichment<-d_malespecenrichment[-nrow(d_malespecenrichment),]
d_malespecenrichment$pval<-malespecenrichment_pval

d_malespecenrichment$enriched<-ifelse(d_malespecenrichment$pval<threshold, d_malespecenrichment$enriched, "Not Significant")
d_malespecenrichment<-data.frame(clusterID=moduledata_short$clusterID, d_malespecenrichment)

# Order by category and save as csv
write.csv(d_malespecenrichment[order(d_malespecenrichment$enriched, d_malespecenrichment$pval),], file = "./Output/malesp_enriched_clusters_e-1.csv")

## Check for clusters enriched in female-specific genes

d_femalespecenrichment<-moduledata_short[,c("femalespec_n", "Nodes")]
module<-1:nrow(d_femalespecenrichment)
d_femalespecenrichment<-rbind(d_femalespecenrichment, apply(d_femalespecenrichment, 2, sum))

femalespecenrichment_pval<-sapply(X = module, FUN = function(x){
  fisher.test(
    rbind(d_femalespecenrichment[x,],d_femalespecenrichment[nrow(d_femalespecenrichment),])
  )$p.value
})
femalespecenrichment_pval<-femalespecenrichment_pval*max(module) # Bonferroni correction

d_femalespecenrichment$prop<-d_femalespecenrichment$femalespec_n/d_femalespecenrichment$Nodes
d_femalespecenrichment$enriched<-sign(d_femalespecenrichment$prop-d_femalespecenrichment$prop[nrow(d_femalespecenrichment)])
d_femalespecenrichment$enriched<-ifelse(sign(d_femalespecenrichment$enriched)>0, "Enriched", "Depleted")

d_femalespecenrichment<-d_femalespecenrichment[-nrow(d_femalespecenrichment),]
d_femalespecenrichment$pval<-femalespecenrichment_pval

d_femalespecenrichment$enriched<-ifelse(d_femalespecenrichment$pval<threshold, d_femalespecenrichment$enriched, "Not Significant")
d_femalespecenrichment<-data.frame(clusterID=moduledata_short$clusterID, d_femalespecenrichment)

# Order by category and save as csv
write.csv(d_femalespecenrichment[order(d_femalespecenrichment$enriched, d_femalespecenrichment$pval),], file = "./Output/femalesp_enriched_clusters_e-1.csv")

## Check for clusters enriched in male-biased genes

d_malebiasenrichment<-moduledata_short[,c("malebias_n", "Nodes")]
module<-1:nrow(d_malebiasenrichment)
d_malebiasenrichment<-rbind(d_malebiasenrichment, apply(d_malebiasenrichment, 2, sum))

malebiasenrichment_pval<-sapply(X = module, FUN = function(x){
  fisher.test(
    rbind(d_malebiasenrichment[x,],d_malebiasenrichment[nrow(d_malebiasenrichment),])
  )$p.value
})
malebiasenrichment_pval<-malebiasenrichment_pval*max(module) # Bonferroni correction

d_malebiasenrichment$prop<-d_malebiasenrichment$malebias_n/d_malebiasenrichment$Nodes
d_malebiasenrichment$enriched<-sign(d_malebiasenrichment$prop-d_malebiasenrichment$prop[nrow(d_malebiasenrichment)])
d_malebiasenrichment$enriched<-ifelse(sign(d_malebiasenrichment$enriched)>0, "Enriched", "Depleted")

d_malebiasenrichment<-d_malebiasenrichment[-nrow(d_malebiasenrichment),]
d_malebiasenrichment$pval<-malebiasenrichment_pval

d_malebiasenrichment$enriched<-ifelse(d_malebiasenrichment$pval<threshold, d_malebiasenrichment$enriched, "Not Significant")
d_malebiasenrichment<-data.frame(clusterID=moduledata_short$clusterID, d_malebiasenrichment)

# Order by category and save as csv
write.csv(d_malebiasenrichment[order(d_malebiasenrichment$enriched, d_malebiasenrichment$pval),], file = file.path(newdir,"malebias_enriched_clusters_e-1.csv"))

## Check for clusters enriched in female-biased genes

d_femalebiasenrichment<-moduledata_short[,c("femalebias_n", "Nodes")]
module<-1:nrow(d_femalebiasenrichment)
d_femalebiasenrichment<-rbind(d_femalebiasenrichment, apply(d_femalebiasenrichment, 2, sum))

femalebiasenrichment_pval<-sapply(X = module, FUN = function(x){
  fisher.test(
    rbind(d_femalebiasenrichment[x,],d_femalebiasenrichment[nrow(d_femalebiasenrichment),])
  )$p.value
})
femalebiasenrichment_pval<-femalebiasenrichment_pval*max(module) # Bonferroni correction

d_femalebiasenrichment$prop<-d_femalebiasenrichment$femalebias_n/d_femalebiasenrichment$Nodes
d_femalebiasenrichment$enriched<-sign(d_femalebiasenrichment$prop-d_femalebiasenrichment$prop[nrow(d_femalebiasenrichment)])
d_femalebiasenrichment$enriched<-ifelse(sign(d_femalebiasenrichment$enriched)>0, "Enriched", "Depleted")

d_femalebiasenrichment<-d_femalebiasenrichment[-nrow(d_femalebiasenrichment),]
d_femalebiasenrichment$pval<-femalebiasenrichment_pval

d_femalebiasenrichment$enriched<-ifelse(d_femalebiasenrichment$pval<threshold, d_femalebiasenrichment$enriched, "Not Significant")
d_femalebiasenrichment<-data.frame(clusterID=moduledata_short$clusterID, d_femalebiasenrichment)

# Order by category and save as csv
write.csv(d_femalebiasenrichment[order(d_femalebiasenrichment$enriched, d_femalebiasenrichment$pval),], file = "./Output/femalebias_enriched_clusters_e-1.csv")

# Merge the four tables

specenrichment <- merge(x = d_femalespecenrichment, y = d_malespecenrichment, by=c("clusterID", "Nodes"), suffixes = c("_fem", "_mal"))
biasenrichment <- merge(x = d_femalebiasenrichment, y = d_malebiasenrichment, by=c("clusterID", "Nodes"), suffixes = c("_fem", "_mal"))
enrichment <- merge (x = specenrichment, y = biasenrichment, by=c("clusterID", "Nodes"), suffixes = c("_spec", "_bias"))

sign_enrichment <- enrichment[,c("clusterID", "Nodes", "enriched_fem_spec", "enriched_mal_spec", "enriched_fem_bias", "enriched_mal_bias")]
sign_enrichment <- sign_enrichment[which(apply(sign_enrichment[,-(1:2)], 1, function(x){any(x!="Not Significant")})),]
sign_enrichment <- sign_enrichment[order(sign_enrichment$enriched_fem_spec, sign_enrichment$enriched_mal_spec, sign_enrichment$enriched_fem_bias, sign_enrichment$enriched_mal_bias),]

sign_enrichment[order(apply(sign_enrichment[,-(1:2)], 1, function(x){sum(as.numeric(as.factor(x)))})),]
table(enrichment$enriched_fem_spec, enrichment$enriched_mal_spec)
table(enrichment$enriched_fem_bias, enrichment$enriched_mal_bias)

## Save table

write.csv(enrichment[,c(1,grep("enriched",names(enrichment)))], file = file.path(newdir, "sex_enriched_clusters.csv"))

# # Check transcription in probelmatic clusters
# eigenexon_evalues <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexon_evalues.csv")
# clusterdata<- read.csv(file = "./Output/WGCNA_clustering_biweight/clusterassignments_simple.csv")
# sexbiased_nodes_short<-read.csv(file = "./Output/sex_biased_nodes/sexbiased_nodes_short_1E-5.csv")[,-1]
# names(sexbiased_nodes_short)<-c("eigenexonID","DE_pattern")
# 
# cluster<-"peru"
# 
# sexbiased_nodes_short[which(sexbiased_nodes_short$eigenexonID%in%clusterdata$transcriptID[which(clusterdata$clusterID==cluster)]),]
# 
# a<-melt(eigenexon_evalues[which(eigenexon_evalues$eigenexonID%in%clusterdata$transcriptID[which(clusterdata$clusterID==cluster)]),])
# a$female<-grepl("female", a$variable)
# a$stage<-as.factor(str_extract(string = a$variable, pattern = "^[^_]*"))
# a$stage<- factor(x = a$stage, levels=c("emb10","emb18","lar51","pupyel","adult"))
# a$replicate<- str_extract(string = a$variable, pattern = "[^.]$")
# a$grouper<-paste(a$eigenexonID, a$replicate, sep="_")
# 
# # ggplot(data = a, aes(x=stage, y=value+0.0000001, col=female, group=grouper))+geom_point(alpha=0.5)+geom_line()+scale_y_log10()
# ggplot(data = a, aes(x=stage, y=value+1e-16, col=female))+geom_boxplot(notch = T)+scale_y_log10()+geom_point(position = "jitter", alpha=0.5)
# ggplot(data = a, aes(x=stage, y=value, col=female))+geom_boxplot(notch = T)+scale_y_log10()+geom_point(position = "jitter", alpha=0.5)

# ## Cluster with IDs navajowhite1 and darkgoldenrod4 have genes with male bias in larvae and female bias in adults, hence their enrichment in both types of biased genes
# pdf(file=file.path(graphdir, "navajowhite1_darkgoldenrod4_dev_expr.pdf"), paper = "a4r")
# multiplot(nw1, dgr4, cols = 2)
# dev.off()