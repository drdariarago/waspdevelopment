###### Compile data from CCRE network

# Clean workspace
rm(list=setdiff(ls(),"ODB8_EukOGs_genes"))
# Load Packages
library(reshape2)
library(plyr)
library(magrittr)
library(stringr)
library(vcd)
library(ggplot2)
library(lattice)
library(fdrtool)
# Initialize output path
newdir<-file.path(getwd(),"Output/Results_compiler_CCRE")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/Results_compiler_CCRE")
dir.create(graphdir)
source(file = "./Scripts/multiplot.R")

# Function to convert pasted output of comma separated entries into single character string with no whitespaces
# input must be QUOTED
trimmer <- function(x){
  trimws(strsplit(x = x, split = ",")[[1]])
}

#### Load datasets
## Load OGS2 genome
NVIT_OGS2_goodannotcomplete <- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")
## Load eigenexon assignments
eigenexons_assignments <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexons_assignments.csv")
## Load CCRE assignments
CCREs <- read.csv(file = "./Output/CCREs/CCRE_assignments.csv")[,-1]
names(CCREs) <- c("eigenexonID", "CCRE_ID", "Representative_Node")
## Load CCRE metadata
CCRE_metadata <- read.csv(file = "./Output/CCREs/CCRE_metadata.csv")[,-1]
## Load cluster assignments
clusterassignments_full <- read.csv("./Output/WGCNA_CCRE_clustering_biweight/CCRE_clusterassignments_full.csv")
names(clusterassignments_full)[1]<-"eigenexonID"
## Load Taxonomic profile (from Sackton, 2013)
tax_depth <- read.delim("./Input/DatasetS1.tsv")
tax_depth <- tax_depth[,c("ogs2", "strata")]
names(tax_depth)<-c("geneID","strata")
tax_depth$strata<-factor(tax_depth$strata, levels = c("Metazoa", "Arthropod", "Insect", "Hymenoptera", "Wasp"))
## Load within-Hymenoptera strata
# swap file to NasStrata_Pteromalid to include unpubl data
NasStrata <- read.csv(file = "./Output/Wasp_Phylostratigraphy/NasStrata.csv", header = T, row.names = 1)

# ## Load sex-specific node assignments
# sexspec <- read.csv("./Output/sex_specific_nodes/sexspec_allstages.csv")[,-1]
# sexspec <- sexspec[,c("eigenexonID","spec")]
# sexspec$spec[which(sexspec$spec=="Aspecific")]<-NA # optional
## Load developmental sexbias expression profiles
devsexbiased_nodes <- read.csv(file = "./Output/sex_biased_nodes_CCRE/sexbiased_nodes_short_5E-2.csv")[,-1]
names(devsexbiased_nodes)<-c("nodeID","devsexbias")
## Load developmental stagebias expression profiles
devstagebiased_nodes <- read.csv(file = "./Output/sex_biased_nodes_CCRE/stagebiased_nodes_short_5E-2.csv")[,-1]
names(devstagebiased_nodes)<-c("nodeID","devstagebias")
# ## Load stage specific expression profile
# stagexpr_nodes <- read.csv(file = "./Output/sex_specific_nodes/stagexpr_nodes_short.csv")[,-1]
# Load all within cluster network concepts
WithinClusterNetworkConcepts_nodes <- read.csv(file = "./Output/KME_correction/WithinCluster_NetworkParams_rescaled.csv")[,-1]

## Load whole-network network concepts
GlobalNetworkConcepts_nodes <- read.csv(file = "./Output/NetworkConcepts_CCRE/CCRE_nodeparams_global.csv")
names(GlobalNetworkConcepts_nodes)[1] <- "nodeID"

## Load precompiled ODB8 counts from commented code at bottom script
ODB8_Nas_EukOGs_genes <- read.csv(file = "./Output/Results_compiler_CCRE/ODB8_Nas_EukOGs_genes.csv")[,-1]

## Load simplified gonad expression annotation
gonad_annotation <- read.csv(file = "./Output/gonads_DE/gonads_DE_summary.csv", row.names = 1)



#### Merge datasets (for each node or CCRE) (multiple nodes per CCREs are representative, check code)
## Merge gene annotation from different sources
NVIT_OGS2<- NVIT_OGS2_goodannotcomplete[,c("geneID","Name","quality7","IntPred","LinkageCluster","Scaffold","recombrate","adult_female_meth_status")]
# NVIT_OGS2<- merge(NVIT_OGS2, tax_depth, by = "geneID")
NVIT_OGS2<- merge(NVIT_OGS2, ODB8_Nas_EukOGs_genes, by = "geneID", all.x=T)


## Merge node annotation from different sources
nodedata<-clusterassignments_full[,c("eigenexonID", "clusterID")]
names(nodedata)[1] <- "CCRE_ID"
nodedata <- merge(CCREs, nodedata, all.y = T)
names(nodedata)[1] <- "nodeID"
nodedata$eigenexonID <- as.factor(ifelse(is.na(nodedata$eigenexonID),
                               as.character(nodedata$nodeID),
                               as.character(nodedata$eigenexonID)))
nodedata$geneID <- str_extract(string = nodedata$eigenexonID, pattern = "^[^_]*")
nodedata <- merge(nodedata, NVIT_OGS2, by="geneID", all.x=F)
# gene Nasvi2EG021272 isn't present in the genome annotation file, removed because of low annotation quality

## Recombination rates are stored as cM/Kb, this line converts cM/Mb
nodedata$recombrate <- nodedata$recombrate*1000
# ## This one converts to Mb/cM
# nodedata$recombrate <- 1/nodedata$recombrate

# Label grey nodes
nodedata$Grey <- as.factor(ifelse(nodedata$clusterID=="grey","grey","cluster"))
# nodedata<-merge(nodedata, sexspec, all.x=T)

nodedata<-droplevels(nodedata)
# Merge developmental sexbias profile
nodedata<-merge(nodedata, devsexbiased_nodes, by = "nodeID", all.x = T, all.y = F)
# Merge developmental expression profile
nodedata<-merge(nodedata, devstagebiased_nodes, by="nodeID", all.x=T, all.y = F)
# Merge taxonomic depth
nodedata <- merge(nodedata, tax_depth, by = "geneID", all.x = T, all.y = F)
# Merge within-hymenoptera stratum
nodedata <- merge(nodedata, NasStrata, by = "geneID", all.x = T, all.y = F)

# Add explicit node type label
nodedata$node_type<- as.factor(
  sapply(X = nodedata$nodeID, FUN = function(x){
    if(grepl(pattern = "CCRE", x = x)) {"CCRE"} else 
      if (grepl(pattern = "_fac", x = x)) {"splicing"} else
      {"transcription"}
  }))
# Add explicit event type label
nodedata$event_type <- as.factor(ifelse(grepl("_fac", nodedata$eigenexonID), "splicing", "transcription"))
## Annotate genes with any form of transcriptional sexbias
sexbgenes <- nodedata[which(nodedata$event_type=="transcription"), c("geneID","devsexbias")]
sexbgenes <- sexbgenes[grep(pattern = "[f,m]", x = sexbgenes$devsexbias),"geneID"]
nodedata$sexbiasedgene <- ifelse(test = nodedata$geneID%in%sexbgenes, "Sexbiased","Unbiased")
## Annotate gonad vs soma and testes vs ovaries expression
gonad_annotation = data.frame(
  eigenexonID = row.names(gonad_annotation),
  gonad_bias = gonad_annotation$gonad_bias
)
nodedata <- merge(nodedata, gonad_annotation, 
  by.x = "nodeID", by.y = "eigenexonID",
  all.x = T, all.y = F)



# Save dataset with all transcripts
transcriptdata <- nodedata
write.csv(x = transcriptdata, file = file.path(newdir, "transcriptdata.csv"))

# # Save dataset with only network nodes
# nodedata <- nodedata[which(nodedata$Representative_Node==T|is.na(nodedata$Representative_Node)),-which(names(nodedata)=="Representative_Node")]
# write.csv(x = nodedata, file = file.path(newdir, "nodedata.csv"))

# Add within network concepts
names(WithinClusterNetworkConcepts_nodes)[-c(1:2)] <- paste(names(WithinClusterNetworkConcepts_nodes)[-c(1:2)], "wc", sep = "_")
nodedata <- merge(nodedata, WithinClusterNetworkConcepts_nodes, by.x = c("nodeID","clusterID"), by.y = c("nodeID","clusterID"), all.x = T, all.y =F)

# Add global network concepts
names(GlobalNetworkConcepts_nodes)[-c(1)] <- paste(names(GlobalNetworkConcepts_nodes)[-c(1)], "Global", sep = "_")
nodedata <- merge(nodedata, GlobalNetworkConcepts_nodes, by.x = c("nodeID"), by.y = c("nodeID"), all.x = T, all.y = F)

# Return devconflict nodes
transcriptdata[grep(pattern = "f.*m|m.*f", x = transcriptdata$devsexbias),]

# Add CCRE metadata
nodedata <-   merge(x = nodedata, y = CCRE_metadata, by.x = "nodeID", by.y = "CCREclust", all.x = T, all.y = F)

## Save as csv
write.csv(x = nodedata, file = file.path(newdir, "transcriptdata_full.csv"))
write.csv(x = nodedata[which(nodedata$Representative_Node==T|is.na(nodedata$Representative_Node)),-which(names(nodedata)=="Representative_Node")], file = file.path(newdir, "nodedata_full.csv"))

### Create module dataset (fix cluster sizes in nodes)
## Load within module network parameters
clusterparams <- read.csv(file = "./Output/NetworkConcepts_CCRE/CCRE_clusterparams.csv")[,-1]
## Load module DE
clusterDE <-  read.csv(file = "./Output/DEcluster_analysis_CCRE/biased_clusters.csv")[,-1]
## Calculate cluster size
clustersize_transcripts <- as.data.frame(table(transcriptdata$clusterID))
names(clustersize_transcripts) <- c("clusterID","nTranscripts")
clustersize_genes <- ddply(.data = transcriptdata, .variables = .(clusterID), summarize, 
                           nGenes = length(unique(geneID)))
clustersize <- merge(clustersize_genes, clustersize_transcripts)
## Calculate cluster size in nodes
clustersize_nodes <- ddply(.data = transcriptdata[which(transcriptdata$Representative_Node==T|is.na(transcriptdata$Representative_Node)),], 
                           .variables = .(clusterID), summarize, 
                           nNodes = length(nodeID))
clustersize <- merge(clustersize, clustersize_nodes)
## and proportion of transcription nodes (normalized to global proportion of splicing nodes)
cluster_spl <- prop.table(table(nodedata$clusterID, nodedata$event_type), margin = 1)[,1]
cluster_spl <- data.frame(clusterID = names(cluster_spl), excessSpl = cluster_spl/prop.table(table(nodedata$event_type))[1])
## add proportion of nodes from duplicated genes
cluster_dup <- prop.table(table(nodedata$clusterID, nodedata$quality7), margin = 1)[,"Paralog"]
cluster_dup <- data.frame(clusterID = names(cluster_dup), excessDup = cluster_dup/prop.table(table(nodedata$quality7))["Paralog"])
## Add proportion of nodes from methylated genes
cluster_met <- prop.table(table(nodedata$clusterID, nodedata$adult_female_meth_status), margin = 1)[,"Methylated"]
cluster_met <- data.frame(clusterID = names(cluster_met), excessMet = cluster_met/prop.table(table(nodedata$adult_female_meth_status))["Methylated"])
## Add enrichment for testes and ovaries genes
# Test if the proportion of testes-biased genes is greater in cluster than in genome

gonad_enrichment <- function(clusterIDs, test_type) {
  nodedata %$%
    table(
      clusterID == clusterIDs,  
      gonad_bias == test_type, 
      dnn = c("In_Cluster", "Biased")) %>%
      {if (extract(., 4) < 5)  {
        NA
      } else {
        fisher.test(., alternative = "greater") %$% 
          p.value
      }}
}


testes = list(levels(nodedata$clusterID)) 
for (c in levels(nodedata$clusterID)) {
  testes[c] = gonad_enrichment(clusterID = c, test_type = "Testes")
}
ovaries = list(levels(nodedata$clusterID)) 
for (c in levels(nodedata$clusterID)) {
  ovaries[c] = gonad_enrichment(clusterID = c, test_type = "Ovaries")
}
gonads = list(levels(nodedata$clusterID)) 
for (c in levels(nodedata$clusterID)) {
  gonads[c] = gonad_enrichment(clusterID = c, test_type = "Gonads")
}

# Compile intod data.frame and apply FDR correction

gonad_cluster_enrichment = data.frame(
  clusterID = names(testes[-1]),
  testes_enrichment = testes[-1] %>% unlist() %>% p.adjust(method = "fdr"),
  ovaries_enrichment = ovaries[-1] %>% unlist() %>% p.adjust(method = "fdr"),
  gonads_enrichment = gonads[-1] %>% unlist() %>% p.adjust(method = "fdr")
)


# ## Add metanetwork parameters
# metanetwork <- read.csv(file = "./Output/Metanetwork_CCRE/AllMetanetworks.csv")[,-1]
# names(metanetwork)[1] <- "clusterID"
# names(metanetwork)[2:ncol(metanetwork)] <- paste("meta_", names(metanetwork)[2:ncol(metanetwork)], sep = "")

## Load information on stage-specific changes in integration (need to code GLM analyses for CCRE network)
DIDC <- read.csv(file = "./Output/sex_modularity_GLMs_CCRE/DIDC_signclusters.csv")[,-1]
## Create average cluster KME, Cluster Coef and Diameter
medianKME <- ddply(.data = nodedata[which(nodedata$Representative_Node==T|is.na(nodedata$Representative_Node)),], .variables = .(clusterID), summarize, 
                   medianKME = median(KME_wc, na.rm = T),
                   medianClusterCoef = median(ClusterCoef_wc, na.rm = T),
                   diameter = max(Betweenness_wc))

## Merge all module data
clusterdata <- merge(clusterparams, clustersize, by = "clusterID", all = T)
clusterdata <- merge(clusterdata, cluster_spl, by = "clusterID")
clusterdata <- merge(clusterdata, cluster_dup, by = "clusterID")
clusterdata <- merge(clusterdata, cluster_met, by = "clusterID")
clusterdata <- merge(clusterdata, medianKME, by = "clusterID", all = T)
# clusterdata <- merge(clusterdata, metanetwork, by = "clusterID", all = T)
clusterdata <- merge(clusterdata, clusterDE, by = "clusterID", all = T)
clusterdata <- merge(clusterdata, DIDC, by = "clusterID", all = T)
clusterdata <- merge(clusterdata, gonad_cluster_enrichment, by = "clusterID", all = T)
head(clusterdata)

## Save as csv
write.csv(x = clusterdata, file = file.path(newdir, "clusterdata_full.csv"))

### Merge module and node data
networkdata <- merge(clusterdata, nodedata, by = "clusterID", all = T)
## Save as csv
write.csv(x = networkdata, file = file.path(newdir, "networkdata.csv"))

# Save lists of geneIDs in each cluster
genelists <- dlply(.data = transcriptdata, .variables = .(clusterID), .fun = summarize, 
      geneIDs = geneID)

dir.create(path = file.path(newdir, "genelists"))

sapply(X = names(genelists), FUN = function(x){
  write(x = as.character(genelists[[x]][,1]), 
        file = file.path(newdir, "genelists", paste(x,".txt", sep = ""))
  )
})
# Universe of geneIDs in final network
write(x = as.character(transcriptdata$geneID), file = file.path(newdir, "genelists", "universe.txt"))


## Restrict only to constitutive nodes
genelists_con <- dlply(.data = transcriptdata[grep(pattern = "_con", x = transcriptdata$eigenexonID),], .variables = .(clusterID), .fun = summarize, 
                   geneIDs = geneID)

dir.create(path = file.path(newdir, "genelists_con"))

sapply(X = names(genelists), FUN = function(x){
  write(x = as.character(genelists_con[[x]][,1]), 
        file = file.path(newdir, "genelists_con", paste(x,".txt", sep = ""))
  )
})
# Universe of geneIDs in final network
write(x = as.character(transcriptdata[grep(pattern = "_con", x = transcriptdata$eigenexonID),]$geneID), file = file.path(newdir, "genelists_con", "universe.txt"))

# Gene counts of sexbiased nodes (at least one per gene)
genesummary <- ddply(
  .data = transcriptdata, 
  .variables = .(geneID), 
  .fun = summarize, 
  malebiased = any(grepl(pattern = "m", x = devsexbias)),
  femalebiased = any(grepl(pattern = "f", x = devsexbias)),
  testesbiased = any(gonad_bias == "Testes"),
  ovariesbiased = any(gonad_bias == "Ovaries"),
  gonadsbiased = any(gonad_bias == "Gonads")
)
genesummary$bothbiased <- with(genesummary, malebiased&femalebiased)
summary(genesummary)
write.csv(x = genesummary, file = file.path(newdir, "genesummary.csv"))


## Store list with top 10 conected and/or betweenness genes in each cluster
# Connectivity
top10K <- ddply(.data = nodedata, .variables = .(clusterID), .fun = summarize, 
                nodeID = nodeID,
                K_rank = length(AbsConnectivity_wc)-rank(AbsConnectivity_wc, ties.method = "max"))
top10K <- top10K[order(top10K$clusterID, top10K$K_rank),]
top10K <- ddply(.data = top10K, .variables = .(clusterID), .fun = summarize, 
                nodeID = nodeID[which(K_rank<10)]
)
# Betweenness
top10B <- ddply(.data = nodedata, .variables = .(clusterID), .fun = summarize, 
                nodeID = nodeID,
                B_rank = length(AbsBetweenness_wc)-rank(AbsBetweenness_wc, ties.method = "max"))
top10B <- top10B[order(top10B$clusterID, top10B$B_rank),]
top10B <- ddply(.data = top10B, .variables = .(clusterID), .fun = summarize, 
                nodeID = nodeID[which(B_rank<10)]
                )
# Merge sets, filter nodedata and order final dataset
top10 <- union(top10B$nodeID, top10K$nodeID)
top10 <- nodedata[which(nodedata$nodeID%in%top10),]
top10 <- top10[order(top10$clusterID, top10$AbsConnectivity_wc, top10$AbsBetweenness_wc),]
top10 <- top10[,c("nodeID","clusterID","geneID","eigenexonID","Name","quality7","devsexbias","devstagebias","AbsConnectivity_wc","AbsBetweenness_wc","KME_wc","ClusterCoef_wc","MAR_wc","CCREsize","adult_female_meth_status")]



clusterdata[grep(pattern = "emb", x = clusterdata),]
summary(clusterdata)

top10[which(top10$clusterID=="brown2"),]
top10[which(top10$clusterID=="grey60"),]
top10[which(top10$clusterID=="honeydew1"),]
top10[which(top10$clusterID=="magenta"),]
top10[which(top10$clusterID=="mediumpurple3"),]
top10[which(top10$clusterID=="orange"),]
top10[which(top10$clusterID=="skyblue"),]

# Save full phylostratigraphic annotation, with cluster assignments (representative nodes only)
StratAnnot <- transcriptdata[, c("geneID", "Representative_Node", "clusterID", "strata", "Stratum")]
StratAnnot <- StratAnnot[which(StratAnnot$Representative_Node==T|is.na(StratAnnot$Representative_Node)),]
StratAnnot$FullStrata <- ifelse(StratAnnot$strata%in%c("Metazoa", "Arthropod", "Insect"), 
  as.character(StratAnnot$strata),
  as.character(StratAnnot$Stratum)
)
StratAnnot$FullStrata <- factor(
  x = StratAnnot$FullStrata, 
  levels = c("Metazoa", "Arthropod", "Insect", "Hymenoptera", "Apocrita", "Chalcid", "Pteromalid", "Nasonia")
)
StratAnnot$FullStrata[which(is.na(StratAnnot$FullStrata))] <- "Nasonia"

StratAnnot <- StratAnnot
table(StratAnnot$FullStrata, useNA = 'ifany')
table(StratAnnot$FullStrata, StratAnnot$clusterID=="grey")

write.csv(
  x = unique(StratAnnot[,c("geneID", "FullStrata")]),
  file = file.path(newdir, "StratumAnnotation.csv")
)
## Improved transcriptdata version, 
# merge all gene annotation metadata
# merge all node annotation metadata
# merge CCRE and nodes 
# merge previous set to gene set

# # General module overview
# ggplot(clusterdata, aes(x = clustersize, y = Density, alpha = Centralization, size = Heterogeneity)) + geom_point() + scale_x_log10() + scale_y_log10() + geom_smooth()
# ggplot(clusterdata, aes(x = clustersize, y = Heterogeneity)) + geom_point() + scale_x_log10() + geom_smooth()
# ggplot(clusterdata, aes(x = clustersize, y = meta_ScaledBetweenness)) + geom_point() + scale_x_log10() + geom_smooth()
# ggplot(clusterdata, aes(x = meta_ScaledConnectivity, y = meta_Betweenness, alpha = meta_MAR, size = clustersize)) + geom_point() + scale_x_log10() + geom_smooth()

# ############# Generate abridged ODB8 version
## Load OrthoDBv8 OG assignments (large file, requires time)
# ODB8_EukOGs_genes <- read.delim("./Input/ODB8_EukOGs_genes_ALL_levels.txt")

# ODB8_Nas_EukOGs_genes<-ODB8_EukOGs_genes[which(ODB8_EukOGs_genes$organism=="Nasonia vitripennis"),c("odb8_og_id","protein_id","odb8_level","uniprot_description")]# select only Nasonia entries, with geneID, ODB ID and uniprot description
# ODB8_Nas_EukOGs_genes<-ODB8_Nas_EukOGs_genes[grep("Arthropoda", ODB8_Nas_EukOGs_genes$odb8_level),c("odb8_og_id","protein_id","uniprot_description")]# select only Arthropoda level
# ODB8_Nas_EukOGs_genes<-ODB8_Nas_EukOGs_genes[order(ODB8_Nas_EukOGs_genes$protein_id),]
# ODB8_Nas_EukOGs_genes$geneID<-str_extract(string = ODB8_Nas_EukOGs_genes$protein_id, pattern = "^[^t]*")
# ODB8_Nas_EukOGs_genes<-droplevels(ODB8_Nas_EukOGs_genes[,c("geneID","odb8_og_id","uniprot_description")])
# ODB8_Copynumber<-data.frame(table(ODB8_Nas_EukOGs_genes$odb8_og_id)) # annotate copy number per OG
# names(ODB8_Copynumber)<-c("odb8_og_id","Copynumber")
# ODB8_Nas_EukOGs_genes<-merge(ODB8_Nas_EukOGs_genes, ODB8_Copynumber, by="odb8_og_id", all.x=F) # add onto OG information
# ## Save as csv 
# write.csv(ODB8_Nas_EukOGs_genes, file = file.path(newdir, "ODB8_Nas_EukOGs_genes.csv"))

######### Load sexbiased node assignments
# sexbiased_nodes <- read.csv("./Output/sex_biased_nodes_CCRE/sexbiased_nodes.csv")
# sexbiased_nodes <- sexbiased_nodes[,c("eigenexonID","sexbias")]
# sexbiased_nodes <- na.exclude(sexbiased_nodes)
# sexbiased_nodes <- ddply(.data = sexbiased_nodes, .variables = .(eigenexonID), .fun = summarize, sexbias= paste(levels(droplevels(sexbias)), collapse=", "))
# sexbiased_nodes$sexbias <- as.factor(sexbiased_nodes$sexbias)

####### Load eigenexon E values
# eigenexon_evalues <- read.csv("./Output/CCREs/CCRE_nasdevgeneexon.csv")