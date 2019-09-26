# Version without dOut
## Compile results from the different analyses
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

## Load sex-specific node assignments
sexspec <- read.csv("./Output/sex_specific_nodes/sexspec_allstages.csv")[,-1]
sexspec <- sexspec[,c("eigenexonID","spec")]
# sexspec$spec[which(sexspec$spec=="Aspecific")]<-NA

## Load sexbiased node assignments
sexbiased_nodes <- read.csv("./Output/sex_biased_nodes/sexbiased_nodes.csv")
sexbiased_nodes <- sexbiased_nodes[,c("eigenexonID","sexbias")]
sexbiased_nodes <- na.exclude(sexbiased_nodes)
sexbiased_nodes <- ddply(.data = sexbiased_nodes, .variables = .(eigenexonID), .fun = summarize, sexbias= paste(levels(droplevels(sexbias)), collapse=", "))
sexbiased_nodes$sexbias <- as.factor(sexbiased_nodes$sexbias)

## Load developmental sexbias expression profiles
devsexbiased_nodes <- read.csv(file = "./Output/sex_biased_nodes/sexbiased_nodes_short_5E-2.csv")[,-1]
names(devsexbiased_nodes)<-c("eigenexonID","devsexbias")

## Load stage expression profile
stagexpr_nodes <- read.csv(file = "./Output/sex_specific_nodes/stagexpr_nodes_short.csv")[,-1]

## Load network concepts
NetConcepts <- read.csv(file = "./Output/Cluster_Coefs/emp_modularconnectivities.csv")[,-1]

# ## Load OrthoDBv8 OG assignments (large file, requires time)
# ODB8_EukOGs_genes <- read.delim("./Input/ODB8_EukOGs_genes_ALL_levels.txt")
# 
# ODB8_Nas_EukOGs_genes<-ODB8_EukOGs_genes[which(ODB8_EukOGs_genes$organism=="Nasonia vitripennis"),c("odb8_og_id","protein_id","odb8_level","uniprot_description")]# select only Nasonia entries, with geneID, ODB ID and uniprot description
# ODB8_Nas_EukOGs_genes<-ODB8_Nas_EukOGs_genes[grep("Arthropoda", ODB8_Nas_EukOGs_genes$odb8_level),c("odb8_og_id","protein_id","uniprot_description")]# select only Arthropoda level
# ODB8_Nas_EukOGs_genes<-ODB8_Nas_EukOGs_genes[order(ODB8_Nas_EukOGs_genes$protein_id),]
# 
# ODB8_Nas_EukOGs_genes$geneID<-str_extract(string = ODB8_Nas_EukOGs_genes$protein_id, pattern = "^[^t]*")
# ODB8_Nas_EukOGs_genes<-droplevels(ODB8_Nas_EukOGs_genes[,c("geneID","odb8_og_id","uniprot_description")])
# ODB8_Copynumber<-data.frame(table(ODB8_Nas_EukOGs_genes$odb8_og_id)) # annotate copy number per OG
# names(ODB8_Copynumber)<-c("odb8_og_id","Copynumber")
# ODB8_Nas_EukOGs_genes<-merge(ODB8_Nas_EukOGs_genes, ODB8_Copynumber, by="odb8_og_id", all.x=F) # add onto OG information
# ## Save as csv 
# write.csv(ODB8_Nas_EukOGs_genes, file = file.path(newdir, "ODB8_Nas_EukOGs_genes.csv"))

# Load precompiled ODB8 counts from previous commented code
ODB8_Nas_EukOGs_genes <- read.csv(file = "./Output/Results_compiler/ODB8_Nas_EukOGs_genes.csv")[,-1]

# Load information on stage-specific changes in integration
dWithin_signClusters<-read.csv(file="./Output/sex_modularity_GLMsummaries/dWithin_signClusters.csv")

# Load information on stage-specific cluster expression bias
sexbiased_clusters_short <- read.csv(file = "./Output/DEcluster_analysis/sexbiased_clusters_short_5E-2.csv")[,-1]

# # load gonad Bias
# gonadannot<-read.csv(file = "./Output/splicing_gonadvalidation/gonadannot_interaction.csv")[,-1]

#### Merge datasets
splicing_cat<-unique(eigenexons_assignments[,c("geneID","splicing_cathegory","clusters")])
splicing_cat<-ddply(splicing_cat, .(geneID), summarize, splicing_cathegory=splicing_cathegory[1], n_isoforms=max(clusters))

NVIT_OGS2<-merge(splicing_cat, NVIT_OGS2_goodannotcomplete, by="geneID", all.x=F)
# gene Nasvi2EG021272 isn't present in the genome annotation file, most likely removed because of low quality

NVIT_OGS2<-merge(NVIT_OGS2, tax_depth, by = "geneID")
NVIT_OGS2<-merge(NVIT_OGS2, ODB8_Nas_EukOGs_genes, by = "geneID", all.x=T)

clusterassignments_short<-clusterassignments_full[,c("eigenexonID", "clusterID")]
clusterassignments_short$geneID<-str_extract(string = clusterassignments_short$eigenexonID, pattern = "^[^_]*")

clusterdata<-merge(clusterassignments_short, NVIT_OGS2[,c("geneID","Name","odb8_og_id","quality7","IntMatch","IntPred","Copynumber","splicing_cathegory","n_isoforms","adult_female_meth_status","strata")], by="geneID")

clusterdata<-merge(clusterdata, sexspec, all.x=T)

clusterdata<-merge(clusterdata, sexbiased_nodes, all.x=T)

clusterdata<-droplevels(clusterdata)

# Merge developmental sexbias profile
clusterdata<-merge(clusterdata, devsexbiased_nodes, by = "eigenexonID", all.x = T)
# Merge developmental expression profile
clusterdata<-merge(clusterdata, stagexpr_nodes, by="eigenexonID", all.x=T)

# Add explicit transcription/splicing label
clusterdata$node_type<-as.factor(ifelse(grepl("_fac", clusterdata$eigenexonID), "splicing", "transcription"))

# Add gonad vs adult biased
clusterdata <- merge(clusterdata, gonadannot, all.x=T, all.y=F)

# Duplications caused by merging contrasts with genes, therefore one row per gene plus one row per sign sex contrast
# solution: convert to one row per node, add status M/F/NA/Both

# # Add developmental changes in integration
# dWithin_signClusters_2 <- dWithin_signClusters[,c("clusterID","Stage","Sex","Estimate")]
# dWithin_signClusters_2$dev_integration_change<-paste(dWithin_signClusters_2$Stage, dWithin_signClusters_2$Sex, sep="_")
# dWithin_signClusters_2 <- dWithin_signClusters_2[,c("clusterID","dev_integration_change","Estimate")]
# names(dWithin_signClusters_2) <- c("clusterID","dev_integration_change","integration_change_coef")
# clusterdata<-merge(clusterdata, dWithin_signClusters_2, by="clusterID", all.x=T, all.y=F)

# Collapse all nodes within a gene that have the same devsexbias and sexbias pattern
duplicatednodes <- clusterdata[which(duplicated(clusterdata[,c("geneID","dev_expr","devsexbias","spec","node_type","gonadbias","tissuebias")])==T),]
clusterdata <- clusterdata[which(duplicated(clusterdata[,c("geneID","dev_expr","devsexbias","spec","node_type","gonadbias","tissuebias")])==F),]

# Remove splicing nodes from within sex-specific genes (can't differentiate between sexes, already expressed across relevant stages)
specgenes <- unique(clusterdata[which(clusterdata$spec%in%c("Male","Female")&clusterdata$node_type=="transcription"),"geneID"])
clusterdata <- clusterdata[-which(clusterdata$geneID%in%specgenes&clusterdata$node_type=="splicing"),]

# Save as .csv
write.csv(clusterdata, file=file.path(newdir, "clusterdata_compiled.csv"))

# cluster_isoform_count<-ddply(clusterdata, .(geneID), summarize, 
#                              nClusters = length(unique(clusterID)), 
#                              nIsoforms = length(unique(eigenexonID)),
#                              .progress = "text")
# 
# densityplot(cluster_isoform_count[which(cluster_isoform_count$nIsoforms>1),]$nClusters/cluster_isoform_count[which(cluster_isoform_count$nIsoforms>1),]$nIsoforms)
# table(cluster_isoform_count$nClusters/cluster_isoform_count$nIsoforms<1)
# prop.table(table(cluster_isoform_count$nClusters/cluster_isoform_count$nIsoforms<1))
# 
# ggplot(cluster_isoform_count, aes(x=nClusters/nIsoforms))+geom_density()+theme_bw()+scale_x_log10()

## Summarize by module
moduledata<-ddply(clusterdata, .(clusterID), summarize, .progress = "text",
                  Nodes=length(unique(eigenexonID)),
                  Genes=length(unique(geneID)),
                  Splicing=sum(grepl("con", eigenexonID)),
                  Trans=sum(grepl("fac", eigenexonID)),
                  Paralogs_OGS=sum(grepl("Paralog",quality7)),
                  Orthologs_OGS=sum(grepl("Ortholog",quality7)),
                  Paralogs_ODB=sum(is.na(Copynumber)==T),
                  Orthologs_ODB=sum(is.na(Copynumber)==F),
                  Meth=sum(adult_female_meth_status=="Methylated", na.rm = T),
                  Unmeth=sum(adult_female_meth_status=="Unmethylated", na.rm = T),
#                   medianRatio=median(ratio, na.rm = T),
#                   medianNresiduals=median(Nresiduals, na.rm = T),
                  taxaNasonia=table(strata)[5],
                  taxaHymenoptera=table(strata)[4],
                  taxaInsecta=table(strata)[3],
                  taxaArthropoda=table(strata)[2],
                  taxaMetazoa=table(strata)[1],
                  sexspec=sum(spec%in%c("Male","Female")),
                  malespec=summary(spec)[2],
                  femalespec=summary(spec)[1],
                  sexbias=sum(sexbias%in%c("Male","Female"))
                  )
moduledata

sxb_moduledata<-ddply(clusterdata, .(clusterID), summarize, 
      malebias=summary(sexbias)[3]+summary(sexbias)[2],
      femalebias=summary(sexbias)[1]+summary(sexbias)[2],
      conflictbias=summary(sexbias)[2]
      )

moduledata<-merge(moduledata, sxb_moduledata, by="clusterID")

# Normalize by number of nodes in cluster
prop_moduledata<-cbind(clusterID = moduledata[,1],moduledata[,-c(1,2,12,13)]/moduledata[,2])
moduledata<-merge(moduledata, prop_moduledata, by="clusterID", suffixes = c("_n","_p"))

# Normalize OGS orthologs/paralogs by number of nodes with annotation in cluster
moduledata$Paralogs_OGS_p2<-moduledata$Paralogs_OGS_n/(moduledata$Paralogs_OGS_n+moduledata$Orthologs_OGS_n)

# Normalize methylation status by number of nodes with annotation in cluster
moduledata$Meth_p2<-moduledata$Meth_n/(moduledata$Meth_n+moduledata$Unmeth_n)

# Save as csv
write.csv(moduledata, file=file.path(newdir, "moduledata_for_enrichment.csv"))

############## Check for enrichment in DE and SS genes as well as spl nodes (run script Cluster_sexspec_sexbias.R and Cluster_splicing_enrichment.R)

# Load and merge results from sex enrichment
enrichment<-read.csv(file = "./Output/Cluster_sexspec_sexbias/sex_enriched_clusters.csv")[,-1]
moduledata<-merge(moduledata, enrichment, by="clusterID", all.x = T, all.y = F)
clusterdata<-merge(clusterdata, enrichment, by = "clusterID", all.x = T, all.y = F)
# Load and merge results from spl enrichment
enrichment<- read.csv(file = "./Output/Cluster_splicing_enrichment/spl_enrichment.csv")[,-1]
moduledata<-merge(moduledata, enrichment[,c("clusterID","enriched_spl")], by="clusterID", all.x = T, all.y = F)
clusterdata<-merge(clusterdata, enrichment[,c("clusterID","enriched_spl")], by = "clusterID", all.x = T, all.y = F)

# ## Load stage-specific sexbias coefficients
# fdr_perm_dWithin <- read.csv("./Output/sex_modularity_GLMs/fdr_perm_dWithin.csv")
# names(fdr_perm_dWithin)[2:3]<-c("Factornames","clusterID")
# 
# # Merge with annotation data
# moduledata<-merge(moduledata, fdr_perm_dWithin, by = "clusterID", all=T)
# clusterdata <- merge(clusterdata, fdr_perm_dWithin, by = "clusterID", all=T)

# Save compiled dataset
write.csv(moduledata, file=file.path(newdir, "moduledata.csv"))
write.csv(clusterdata, file = file.path(newdir, "clusterdata.csv"))

############# Merge integration data
dWithin_signClusters$dWithinBias<-as.factor(paste(dWithin_signClusters$Stage, dWithin_signClusters$Sex, sep = "_"))
dWithin_signClusters<-dWithin_signClusters[,c("clusterID","dWithinBias")]
clusterdata<-merge(clusterdata, dWithin_signClusters, by = "clusterID", all.x = T)
moduledata<-merge(moduledata, dWithin_signClusters, by = "clusterID", all.x = T)

####### Merge cluster bias data
clusterdata <- merge(clusterdata, sexbiased_clusters_short, by = "clusterID", all.x = T)
moduledata<-merge(moduledata, sexbiased_clusters_short, by = "clusterID", all.x = T)

## Merge Network Concepts
clusterdata2 <- merge(clusterdata, NetConcepts, by="eigenexonID")

# # Remove gonadbiased genes from femalebiased ones
# clusterdata$devsexbias <- as.character(clusterdata$devsexbias)
# clusterdata$devsexbias[which(grepl(pattern = "[.]{4}f$", x = clusterdata$devsexbias)&grepl(pattern = "Testes", x = clusterdata$gonadbias))] <- "....g"
# clusterdata$devsexbias <- as.factor(clusterdata$devsexbias)

# Save full compiled dataset
write.csv(clusterdata, file=file.path(newdir, "clusterdata_full.csv"))
write.csv(moduledata, file = file.path(newdir, "moduledata_full.csv"))

# Save dataset restricted to nodes in clusters with sex-biased changes in dWithin
write.csv(clusterdata[which(is.na(clusterdata$dWithinBias)==F),], file=file.path(newdir, "dWithin_signNodes.csv"))

# Save gene sets of clusters with significant changes in dWithin
write.csv(unique(clusterdata[which(is.na(clusterdata$dWithinBias)==F),c("clusterID","geneID","dWithinBias")]), file=file.path(newdir, "dWithin_signGenes.csv"))

# Save gene ID lists
# dWithinSignClusters2<-unique(clusterdata[which(is.na(clusterdata$dWithinBias)==F),c("clusterID","geneID","dWithinBias")])
# dWithinSignClusters2$ID<-paste(dWithinSignClusters2$dWithinBias, dWithinSignClusters2$clusterID, sep = "_")
# dWithinSignClusters2<-dWithinSignClusters2[,c("ID","geneID")]
# # sapply(unique(dWithinSignClusters2$ID), function(x){write(dWithinSignClusters2[which(dWithinSignClusters2$ID==x),"geneID"], file = file.path(newdir, paste(x,"txt",sep=".")))})
# 
# # Write cluster sizes in separate file
# write.csv(as.data.frame(table(dWithinSignClusters2$ID)), file.path(newdir, "dWithin_ClusterSizes.csv"))

# Save concatenated gene ID lists in one fasta file (each annotated by >clusterID, bias and size)
dWithin_signClusters2<-merge(dWithinSignClusters2,as.data.frame(table(dWithinSignClusters2$ID)), by.x="ID", by.y="Var1", all.x=T, all.y=F)
dWithin_signClusters2$ID<-as.factor(paste(dWithin_signClusters2$ID, dWithin_signClusters2$Freq, sep="_"))
headers<-unique(dWithin_signClusters2$ID)


dWithin_signClusters2<- dWithin_signClusters2[,-3]
dWithin_signClusters2<-dlply(dWithin_signClusters2, .variables = .(ID), .fun = summarize, geneID=geneID)
# Note, file is edited after printing to remove the extra output and conform to the list of gene IDs separated by >listName
writeLines(file=file.path(newdir, "test.txt"))
dWithin_signClusters2
sink()

# Print all clusters
allcluster_genelist<-dlply(.data = unique(clusterdata[,c("clusterID","geneID")]), .variables = .(clusterID), .fun = summarize, geneID=geneID)
sink(file = file.path(newdir, "all_modules_enrichment.txt"))
print(allcluster_genelist)
sink()

# Print only clusters with significant dI, merged within each stage:sex contrast
dIclusters_genelist<-dlply(.data = clusterdata[which(is.na(clusterdata$dWithinBias)==F),], .variables = .(dWithinBias), .fun = summarize, geneID=unique(geneID))
sink(file = file.path(newdir, "dIclusters_genelist.txt"))
print(dIclusters_genelist)
sink()

# How many clusters have sex-biased changes in dWithin?
lFDR_threshold<-0.05

length(unique(fdr_perm_dWithin[fdr_perm_dWithin$lFDR_adjusted<lFDR_threshold,"clusterID"]))
length(unique(fdr_perm_dWithin[fdr_perm_dWithin$lFDR_adjusted<lFDR_threshold,"clusterID"]))/length(levels(fdr_perm_dWithin$clusterID))

# Plot dI clusters by stage
moduledata$Factornames<-factor(moduledata$Factornames, levels = c("stageemb10:sexbias","stageemb18:sexbias","stagelar51:sexbias","stagepupyel:sexbias","stageadult:sexbias"))

ggplot(moduledata[which(moduledata$lFDR_adjusted<lFDR_threshold),], aes(x=Factornames, fill=as.factor(sign(Estimate)), group=sign(Estimate)))+geom_histogram(position="dodge")+theme_bw()+scale_x_discrete(labels=c("emb10","emb18","larva","pupa","adult"))


## Create data.frame with stages with significant contrasts for each cluster
# significant_contrasts<-ddply(fdr_perm_data, .(clusterID), summarize, stage = Factornames[which(lFDR_adjusted_dOut<lFDR_threshold|lFDR_adjusted_dWithin<lFDR_threshold)], .progress = "text")

significant_contrasts<-ddply(fdr_perm_dWithin, .(clusterID), summarize, 
                             stagedW = ifelse(any(lFDR_adjusted<lFDR_threshold),Factornames[which(lFDR_adjusted<lFDR_threshold)],NA),
                             sexdW = ifelse(any(lFDR_adjusted<lFDR_threshold),sign(Estimate[which(lFDR_adjusted<lFDR_threshold)]),NA), 
                             .progress = "text")

# Select only one row for each cluster, dividing between clusters with at least one significant stage:sex interaction and no sexbias
moduledata_short<-ddply(moduledata, .(clusterID), summarize, signif = any(lFDR_adjusted<lFDR_threshold), .progress = "text")
moduledata_short<-unique(merge(moduledata_short, moduledata[,c("clusterID","Nodes","Genes_n","Splicing_p","Trans_p","Paralogs_OGS_p2","Paralogs_OGS_p","Orthologs_OGS_p","Meth_p","Unmeth_p","Meth_p2","medianRatio","medianNresiduals","taxaNasonia_p","taxaHymenoptera_p","taxaInsecta_p","taxaArthropoda_p","taxaMetazoa_p","taxaNasonia_n","taxaHymenoptera_n","taxaInsecta_n","taxaArthropoda_n","taxaMetazoa_n","sexspec_n","malespec_n","femalespec_n","sexbias_n","malebias_n","femalebias_n","sexspec_p","malespec_p","femalespec_p","sexbias_p","malebias_p","femalebias_p")], by = "clusterID", all.x=F, all.y=F))

# Add stage and direction of significance
moduledata_short<-merge(moduledata_short, significant_contrasts, by="clusterID")

# Save as csv file
write.csv(moduledata_short, file=file.path(newdir, "moduledata_short.csv"))

# Save list of nodes present in network (to use as universe set for enrichment tests)
write(x = unique(clusterdata$geneID), file = file.path(newdir, "geneIDs_fullnetwork.txt"))

# Save random cluster assignments

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