## Perform enrichment test for splicing nodes
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
newdir<-file.path(getwd(),"Output/Cluster_splicing_enrichment")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/Cluster_splicing_enrichment")
dir.create(graphdir)
source(file = "./Scripts/multiplot.R")

## Load dataset

moduledata_short<-read.csv(file = "./Output/Results_compiler/moduledata_for_enrichment.csv")[,-1]

# Perform enrichment test for sxb nodes on clusters
d_splenrichment<-moduledata_short[,c("Splicing_n", "Nodes")]
module<-1:nrow(d_splenrichment)
d_splenrichment<-rbind(d_splenrichment, apply(d_splenrichment, 2, sum))

sxbenrichment_pval<-sapply(X = module, FUN = function(x){
  fisher.test(
    rbind(d_splenrichment[x,],d_splenrichment[nrow(d_splenrichment),])
  )$p.value
})
sxbenrichment_pval<-sxbenrichment_pval*max(module) # Bonferroni correction

d_splenrichment$prop<-d_splenrichment$Splicing_n/d_splenrichment$Nodes
d_splenrichment$enriched<-sign(d_splenrichment$prop-d_splenrichment$prop[nrow(d_splenrichment)])
d_splenrichment$enriched<-ifelse(sign(d_splenrichment$enriched)>0, "Enriched", "Depleted")

d_splenrichment<-d_splenrichment[-nrow(d_splenrichment),]
d_splenrichment$pval<-sxbenrichment_pval

d_splenrichment$enriched<-ifelse(d_splenrichment$pval<0.01, d_splenrichment$enriched, "Not Significant")
d_splenrichment<-data.frame(clusterID = moduledata_short$clusterID, d_splenrichment)
names(d_splenrichment)[5]<-"enriched_spl"
d_splenrichment<-d_splenrichment[order(d_splenrichment$clusterID, d_splenrichment$enriched_spl),c("clusterID","enriched_spl")]
# Order by category and save as csv
write.csv(x = d_splenrichment, file = file.path(newdir, "spl_enrichment.csv"))
