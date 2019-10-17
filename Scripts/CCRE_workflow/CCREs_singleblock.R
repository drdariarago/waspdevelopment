### Same as CCREs but without splitting dataset in blocks for CCRE calculations
## Collapse all nodes with >95% correlation before DE and network construction analyses
## This generates a list of CCREs (Constitutively Coexpressed Reguatory Events)

## Load libraries and clean space
date()
rm(list=ls())
library(plyr)
library(stringr)
library(WGCNA)
library(ggplot2)
allowWGCNAThreads()
library(reshape2)
library(flashClust)
source(file = "./Scripts/multiplot.R")
# logit<-function(x){
#   log10(x/(1-x))
# }
z<-function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
}


# initialize output path
newdir<-file.path(getwd(), "Output/CCREs")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/CCREs")
dir.create(graphdir)

# import dataset
nasdevgeneexon<-read.csv(file="./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexon_evalues.csv")
row.names(nasdevgeneexon)<-nasdevgeneexon$eigenexonID
# nasdevgeneexon<-nasdevgeneexon[1:1000,] # reduce size for testing
nasdevgeneexon<-nasdevgeneexon[,which(colnames(nasdevgeneexon)%in%setdiff(names(nasdevgeneexon),"eigenexonID"))]
t_nasdevgeneexon<-transposeBigData(nasdevgeneexon, blocksize = 50000)

### Apply noise thresholding and normalize
## Separate splicing and transcription for normalization procedures
m_nasdevgeneexon<-data.frame(melt(t_nasdevgeneexon),
                             sampleID=row.names(t_nasdevgeneexon))
m_nasdevgeneexon$con<-grepl("_con", m_nasdevgeneexon$variable)

ddply(.data = m_nasdevgeneexon, .variables = .(con), summarize, 
      mean=mean(value, na.rm = T),
      sd=sd(value, na.rm = T),
      zero=(0-mean(value, na.rm = T))/sd(value, na.rm = T),
      min = min(value, na.rm = T),
      max = max(value, na.rm = T)
)

# Rescale constitutive nodes in 0-1 space relative to most expressed node
m_nasdevgeneexon[which(m_nasdevgeneexon$con==T),"value"] <- 
  m_nasdevgeneexon[which(m_nasdevgeneexon$con==T),"value"]/max(m_nasdevgeneexon[which(m_nasdevgeneexon$con==T),"value"], na.rm = T)

# Recast into sample x node format
t_nasdevgeneexon<-recast(m_nasdevgeneexon[,1:3], sampleID~variable)
row.names(t_nasdevgeneexon)<-t_nasdevgeneexon$sampleID
t_nasdevgeneexon<-t_nasdevgeneexon[,-1]

# QC filtering of genes with low variance or low expression
# Check for zero variance nodes and nodes with too many missing values
goodCols<-goodSamplesGenes(datExpr = t_nasdevgeneexon, minNSamples = 3, minFraction = 1/10)
lapply(goodCols, function(x) {table(x)})
lapply(goodCols, function(x) {prop.table(table(x))})
# Save discarded nodes 
discarded <- t_nasdevgeneexon[,which(goodCols$goodGenes==F)]
write.csv(x = discarded, file = file.path(newdir, "discarded_nodes.csv"))
# Remove discarded nodes
t_nasdevgeneexon<-t_nasdevgeneexon[,which(goodCols$goodGenes==T)]
# How many gene and how many transcript values were removed?
table(grepl("con", row.names(nasdevgeneexon)))-table(grepl("con", colnames(t_nasdevgeneexon)))

# # Print distributions
# pdf(file = file.path(graphdir, "distr_scaling.pdf"), paper = "a4r", useDingbats = F)
# distr_post<-ggplot(data=m_nasdevgeneexon_z, aes(x=value, col=con))+geom_density()+theme_bw()
# multiplot(distr_pre, distr_post, cols = 1)
# dev.off()

# Save transposed and normalized expression dataset
save(t_nasdevgeneexon, file = file.path(newdir, "t_nasdevgeneexon.RData"))
write.csv(nasdevgeneexon, file = file.path(newdir, "nasdevgeneexon.csv"))

### CCRE computation
### Create correlation matrix with bicor function 
# cormat <- bicor(t_nasdevgeneexon,  quick = 0, use = "pairwise.complete.obs")
cormat <- bicor(t_nasdevgeneexon,  quick = 0, use = "all")

## Save correlation matrix
save(cormat, file=file.path(newdir,"cormat.RData"))
# load(file = "./Output/CCREs/cormat.RData")

## Create distance matrix for tree
# Create binary distance matrix with only distances above threshold ranking as 1 and convert to tree
CCREclust <- hclust(d = as.dist(cormat < .95), method = "complete")
# Create list of node to cluster 
CCREclust <- cutree(tree = CCREclust, h = 0.05)
# Subselect only CCREs with >1 member
CCREs <- names(table(CCREclust))[(table(CCREclust)>1)]
CCREs <- CCREclust[which(CCREclust%in%CCREs)]
# Name CCREs
CCREclust <- paste("CCRE_", as.numeric(as.factor(CCREs)), sep = "")
names(CCREclust) <- names(CCREs)
# Reduce dataset to CCREs
CCRE_nasdevgeneexon <- t(t_nasdevgeneexon[,which(colnames(t_nasdevgeneexon)%in%names(CCREs))])
# Add groupings
CCRE_nasdevgeneexon <- merge(as.data.frame(CCREclust), CCRE_nasdevgeneexon, by="row.names")
row.names(CCRE_nasdevgeneexon) <- CCRE_nasdevgeneexon$Row.names
CCRE_nasdevgeneexon <- CCRE_nasdevgeneexon[,-grep("Row.names", names(CCRE_nasdevgeneexon))]
# Calculate reduced set
CCRE_list <- collapseRows(datET = CCRE_nasdevgeneexon[,-1], 
                          rowGroup = CCRE_nasdevgeneexon$CCREclust, 
                          rowID = row.names(CCRE_nasdevgeneexon), 
                          connectivityBasedCollapsing = TRUE, 
                          method = "MaxMean")
save(CCRE_list, file = file.path(newdir, "CCRE_list.RData"))
# Create dataset of non-CCRE nodes
nasdevgeneexon <- t(t_nasdevgeneexon)
nasdevgeneexon <- nasdevgeneexon[which((row.names(nasdevgeneexon)%in%names(CCREclust))==F),]
# Add CCREs to eigenexons and save as compiled dataset
CCRE_expr <- CCRE_list$datETcollapsed
CCRE_nasdevgeneexon <- rbind(nasdevgeneexon, CCRE_expr)
write.csv(CCRE_nasdevgeneexon, file=file.path(newdir, "CCRE_nasdevgeneexon.csv"))
# Save list of gene used as representative one for the CCRE
CCREclust <- merge(as.data.frame(CCREclust),as.data.frame(CCRE_list$selectedRow), by="row.names")
CCREclust <- CCREclust[order(CCREclust$CCREclust,CCREclust$Row.names),]
write.csv(x = CCREclust, file = file.path(newdir, "CCRE_assignments.csv"))

save.image(file = file.path(newdir, "CCREs_workspace.RData"))

# Annotate CCRE size, proportion of transcripts and number of genes
# CCREclust <- read.csv(file = "./Output/CCREs/CCRE_assignments.csv")[,-1]
CCRE_metadata <- ddply(.data = CCREclust, .variables = .(CCREclust), .fun = summarize, 
                       CCREclust = CCREclust[1],
                       CCREsize = length(Row.names),
                       # CCREtrans = ifelse(any(grepl(pattern = "fac",x = Row.names)),sum(grepl(pattern = "fac", Row.names))/length(Row.names),0),
                       CCREtrans = sum(grepl(pattern = "fac", Row.names)),
                       CCREgenes = length(unique(str_extract(string = Row.names, pattern = "^[^_]*")))
)
write.csv(x = CCRE_metadata, file = file.path(newdir, "CCRE_metadata.csv"))

## Annotate CCREs with gene IDs and gene names
CCREannot <- CCREclust
names(CCREannot) <- c("nodeID","CCREID","RepresentativeNode")
CCREannot$geneID <- str_extract(string = CCREannot$nodeID, pattern = "^[^_]*")
# Load genome annotation
NVIT_OGS2_goodannotcomplete <- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")[,-1]
NVIT_annot <- NVIT_OGS2_goodannotcomplete[,c("geneID","Name")]
# Merge
CCREannot <- merge(CCREannot, NVIT_annot, by = "geneID", all.x = F, all.y = F)
# Reorder rows by CCRE/representative/gene
CCREannot <- CCREannot[order(CCREannot$CCREID, 1-CCREannot$RepresentativeNode, CCREannot$geneID, decreasing = F),]
# Reorder columns by CCRE/node/gene/repr
CCREannot <- CCREannot[,c("CCREID","geneID","nodeID","RepresentativeNode","Name")]
# Save as csv
write.csv(x = CCREannot, file = file.path(newdir, "CCRE_annotation.csv"))

# #### Post-hoc checks for CCREs with <1 selected row
# nodetable <- merge(CCREclust, CCRE_list$selectedRow, by.x="Row.names", by.y = "row.names")
# nodetable <- merge(nodetable, CCRE_list$group2row, by.x="x", by.y="group")
# nodetable <- nodetable[order(nodetable$x, nodetable$y),]
# write.csv(nodetable, file = file.path(newdir, "nodetable.csv"))
# table(nodetable$x, nodetable$y)
# CCRE_nasdevgeneexon <- CCRE_nasdevgeneexon[match(row.names(CCRE_nasdevgeneexon), names(CCREclust)),] # order according to list of CCREs

# ## Convert to binary distance items
# cormat <- lapply(x = cormat, function(x){
#   as.dist( x < .95 )
# })  