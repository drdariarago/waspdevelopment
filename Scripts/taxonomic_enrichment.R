#### Calculate taxonomic level of module enrichment
# Ideal method Metazoa/Arthropoda/Insecta/Holometabola/Hymenoptera/Nasonia
# Use series of GLMs with gamma distribution (frequency) and resample labels for bootstrapping

## Initialize script
rm(list=ls())
library(stringr)
library(plyr)

## Load eigenexon assignments

## Load phyletic profiles
OrthoDB7_LEVELS <- read.delim("./Input/OrthoDB7_LEVELS_tabtext")
OrthoDB7_levels <- str_extract_all(string = OrthoDB7_LEVELS$Species, pattern = "[[:alpha:]]{5}")[1:295]

Taxa<-list(
  Hymenoptera=OrthoDB7_levels[which(grepl("NVITR",OrthoDB7_levels)&grepl("AMELL",OrthoDB7_levels))],
#   Hymenoptera=Hymenoptera[[which.min(sapply(Hymenoptera, length))]],
  Holometabola=OrthoDB7_levels[which(grepl("NVITR",OrthoDB7_levels)&grepl("DMELA",OrthoDB7_levels))],
#   Holometabola=Holometabola[[which.min(sapply(Holometabola, length))]],
  Insecta=OrthoDB7_levels[which(grepl("NVITR",OrthoDB7_levels)&grepl("ZNEVA",OrthoDB7_levels))],
#   Insecta=Insecta[[which.min(sapply(Insecta, length))]],
  Arthropoda=OrthoDB7_levels[which(grepl("NVITR",OrthoDB7_levels)&grepl("ISCAP",OrthoDB7_levels))],
  Metazoa=OrthoDB7_levels[[which(grepl("NVITR",OrthoDB7_levels)&grepl("HSAP",OrthoDB7_levels))]]
)

Taxa[-5]<-lapply(Taxa[-5], function(x){x<-x[which.min(sapply(x,length))]})
Taxa[-5]<-lapply(Taxa[-5], function(x){x<-x[[1]]})

Taxa2<-list(
  Hymenoptera=setdiff(Taxa$Hymenoptera,"NVITR"),
  Holometabola=setdiff(Taxa$Holometabola,Taxa$Hymenoptera),
  Insecta=setdiff(Taxa$Insecta,Taxa$Holometabola),
  Arthropoda=setdiff(Taxa$Arthropoda, Taxa$Insecta),
  Metazoa=setdiff(Taxa$Metazoa, Taxa$Arthropoda)
)

Taxa<-Taxa2

# # Remove Nasonia from the match table
# Taxa<-lapply(X = Taxa, FUN = function(x){x[-grep("NVITR",x)]})

# Import metazoan OGs {testing version, only Arthropoda}
OrthoDB7_Arthropoda <- read.delim("./Input/OrthoDB7_Arthropoda_tabtext", dec=",")

# Select only [OGs containing Nasonia genes,OG_ID, Species_ID, Gene_ID]
OrthoDB7_OGs <- levels(droplevels(OrthoDB7_Arthropoda[grep("Nasonia vitripennis",OrthoDB7_Arthropoda$Organism),"ODB7_OG_ID"]))
OrthoDB7_Arthropoda <- OrthoDB7_Arthropoda[which(OrthoDB7_Arthropoda$ODB7_OG_ID%in%OrthoDB7_OGs),c("ODB7_OG_ID","Gene_ID","Organism","UniProt_Species")]

# Attach ODB code
OrthoDB7_SPECIES <- read.delim("./Input/OrthoDB7_SPECIES_tabtext")
OrthoDB7_SPECIES <- OrthoDB7_SPECIES[,c("ODB_code","Organism")]
OrthoDB7_Arthropoda <- merge(OrthoDB7_Arthropoda, OrthoDB7_SPECIES, by.x="Organism", by.y="Organism", all.x=T, all.y=F)

# Create list of species per OG
OG_species_list<-dlply(OrthoDB7_Arthropoda, .(ODB7_OG_ID), .progress = "text", .fun = summarize, ODB_code=levels(droplevels(ODB_code)))
OG_species_list<-lapply(OG_species_list, FUN = function(x){x<-x[,1]})

# Scan OGs for presence of species from a taxonomic level

Taxascan<-data.frame(
  ODB7_OG_ID=names(OG_species_list),
  Hymenoptera=laply(OG_species_list,function(x){any(x%in%Taxa$Hymenoptera)}),
  Holometabola=laply(OG_species_list,function(x){any(x%in%Taxa$Holometabola)}),
  Insecta=laply(OG_species_list,function(x){any(x%in%Taxa$Insecta)}),
  Arthropoda=laply(OG_species_list,function(x){any(x%in%Taxa$Arthropoda)}),
  Metazoa=laply(OG_species_list,function(x){any(x%in%Taxa$Metazoa)})
)

# Select furthest match

Taxascan_max<-apply(Taxascan[,-1], 1, function(x){max(which(x==T))})
Taxascan<-cbind(Taxascan, Taxascan_max)

# Create Nasonia geneID/Taxonomic depth table
OGS2_Taxascan<-OrthoDB7_Arthropoda[which(OrthoDB7_Arthropoda$ODB_code=="NVITR"),c("Gene_ID","ODB7_OG_ID")]
OGS2_Taxascan<-merge(OGS2_Taxascan, Taxascan, all.x=T, all.y=F)


########Outdated code

# ## Same but for Insecta
# # Import Insecta OGs {testing version, only Arthropoda}
# OrthoDB7_Insecta <- read.delim("./Input/OrthoDB7_Insecta_tabtext", dec=",")
# 
# # Select only [OGs containing Nasonia genes,OG_ID, Species_ID, Gene_ID]
# OrthoDB7_OGs <- levels(droplevels(OrthoDB7_Insecta[grep("Nasonia vitripennis",OrthoDB7_Insecta$Organism),"ODB7_OG_ID"]))
# OrthoDB7_Insecta <- OrthoDB7_Insecta[which(OrthoDB7_Insecta$ODB7_OG_ID%in%OrthoDB7_OGs),c("ODB7_OG_ID","Gene_ID","Organism","UniProt_Species")]
# 
# # Attach ODB code
# OrthoDB7_Insecta <- merge(OrthoDB7_Insecta, OrthoDB7_SPECIES, by.x="Organism", by.y="Organism", all.x=T, all.y=F)
# 
# # Create list of species per OG
# OG_species_list_Insecta<-dlply(OrthoDB7_Insecta, .(ODB7_OG_ID), .progress = "text", .fun = summarize, ODB_code=levels(droplevels(ODB_code)))
# OG_species_list_Insecta<-lapply(OG_species_list_Insecta, FUN = function(x){x<-x[,1]})
# 
# # Scan OGs for presence of species from a taxonomic level
# 
# Taxascan_Insecta<-data.frame(
#   ODB7_OG_ID=names(OG_species_list_Insecta),
#   Hymenoptera=laply(OG_species_list_Insecta,function(x){any(x%in%Taxa$Hymenoptera)}),
#   Holometabola=laply(OG_species_list_Insecta,function(x){any(x%in%Taxa$Holometabola)}),
#   Insecta=laply(OG_species_list_Insecta,function(x){any(x%in%Taxa$Insecta)}),
#   Arthropoda=laply(OG_species_list_Insecta,function(x){any(x%in%Taxa$Arthropoda)})
# )
# 

### Outdated code
# ODB8_EukOGs_genes_Arthropoda <- read.delim("./Input/ODB8_EukOGs_genes_Arthropoda-6656.txt")
# ODB8_EukOGs_genes_Arthropoda<-ODB8_EukOGs_genes_Arthropoda[which(ODB8_EukOGs_genes_Arthropoda$organism=="Nasonia vitripennis"),]
