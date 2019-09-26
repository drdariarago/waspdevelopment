# Taxonomic enrichment OrthoDB8
# Initialize script
rm(list=ls())
library(stringr)
library(plyr)
library(taxize)
# initialize output path
newdir<-file.path(getwd(), "Output/Taxonomic_enrichment_OrthoDB8")
dir.create(newdir)

# Load data
ODB8_Arthropoda <- read.delim("./Input/ODB8_EukOGs_genes_Arthropoda-6656.txt")

