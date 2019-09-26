### early genes
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(limma)
library(fdrtool)
library(MuMIn)
rm(list = ls())
# initialize output path
newdir<-file.path(getwd(), "Output/sexbias_sexspec_early")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/sexbias_sexspec_early")
dir.create(graphdir)

# Load data
clusterdata <- read.csv(file = "./Output/Results_compiler/clusterdata_full.csv")[,-1]
genedata <- read.csv(file = "./Output/sex_biased_genes_summary/genedata.csv")[,-1]

# Subset early nodes
nodes_early <- grep(pattern = "[f|m]", x = strtrim(clusterdata$devsexbias, width = 3))
nodes_early <- union(nodes_early, 
                     which(grepl(pattern = "1", x = strtrim(clusterdata$dev_expr, width = 3))&is.na(clusterdata$spec)==F)
)
nodes_early <- clusterdata[nodes_early,]
# Subset early genes
genes_early <- genedata[which(genedata$geneID%in%nodes_early$geneID),]

# Plot regulatory modes as venn diagram
summary(genes_early)
sexregulated_genes <- data.frame(
  spectrans = is.na(genes_early$sexspec_transcription)==F, 
  biastrans = is.na(genes_early$sexbias_transcription)==F, 
  biasspl = genes_early$sexbias_splicing%in%c("Confbias","Femalebias","Malebias"), 
  specspl = genes_early$sexspec_splicing%in%c("Confspec","Femalespec","Malespec"))

vennDiagram(vennCounts(sexregulated_genes), main = "\n Genes with Early Sex-Linked Regulation",cex = 1.2,
            names = c("Specific\nTranscription", "Biased\nTranscription", "Biased\nSplicing", "Specific\nSplicing")
)

# Most frequent bias&expr patterns among sexbiased

# Most frequent expr patterns among specnodes
specnodes_early <- droplevels(nodes_early[which(is.na(clusterdata$spec)==F),])
table(specnodes_early$dev_expr, specnodes_early$spec, specnodes_early$node_type)

table(grepl("^1", specnodes_early$dev_expr),grepl("^.1",specnodes_early$dev_expr), useNA = "ifany", dnn = c("emb10","emb18"))
mosaic(table(grepl("^1", specnodes_early$dev_expr),grepl("^.1",specnodes_early$dev_expr), useNA = "ifany", dnn = c("emb10","emb18")), shade = T)
