# ODB7 Taxonomic profile, highest and lowest level

## Initialize script
rm(list=ls())
library(stringr)
library(plyr)

# initialize output path
newdir<-file.path(getwd(), "Output/BEAR_ODB7_Taxonomic_enrichment")
dir.create(newdir)

# Import metazoan OGs
OrthoDB7_Metazoa <- read.delim("./Input/OrthoDB7_ALL_METAZOA_tabtext", dec=",")

# Select Nasonia genes
OrthoDB7_Nasonia<-OrthoDB7_Metazoa[which(OrthoDB7_Metazoa$Organism=="Nasonia vitripennis"),]

# Convert taxonomic level to number
OrthoDB7_Nasonia$ODB7_Level<-as.numeric(as.character(OrthoDB7_Nasonia$ODB7_Level))

# For every gene, pick lowest and highest taxonomic level
OrthoDB7_Summary<-ddply(OrthoDB7_Nasonia, .(Gene_ID), summarize, maxLevel=max(ODB7_Level), minLevel=min(ODB7_Level))

# And save as csv
write.csv(OrthoDB7_Summary, file=file.path(newdir, "OrthoDB7_Taxonsummary.csv"))
write.csv(OrthoDB7_Nasonia, file=file.path(newdir, "OrthoDB7_Nasonia.csv"))