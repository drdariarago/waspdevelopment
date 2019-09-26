#### Check oldest matches of Nasonia genes, using the Orthologous Groups from Lindsey et al. 2018 (Trichogramma genome)
rm(list=ls())

# Initialize output path
newdir<-file.path(getwd(),"Output/Wasp_Phylostratigraphy")
dir.create(newdir)

library(stringr)
library(plyr)

# Set taxonomic strata of interest
# Nasonia Only
# Pteromalid: Muscidifurax, Trichomalopsis?
# Chalcid: Trichogramma
# Apocrita: Apis, Micropolis
# Hymenoptera: Orossus, Athalia

Nasonia = "NAVI"
Pteromalid = c("MURA", "TRSA")
Chalcid = c("CESO", "COFL", "TPRE")
Apocrita = c("APME", "MIDE")
Hymenoptera = c("AROS", "ORAB")

# Import Taxonomic information for wasp genes
OGs = read.csv(file = "./Input/OrthoMCL/Amelia/hymenoptera_orthogroups_Lindsey2018BMCBiol.txt", header = F, sep = ":" ,stringsAsFactors = F)
names(OGs) = c("OG_ID", "Data")

# Create list of orthologous assignment per each nasonia gene
NasGenes = str_extract_all(string = OGs$Data, pattern = "Nasvi[A-z,0-9]*")
names(NasGenes) = OGs$OG_ID
NasGenes = stack(NasGenes)
names(NasGenes) = c("geneID", "OG_ID")

#### Integrate further matches from Amelia's trichogramma paper

# Import Trichogramma genes with other Hymenopteran matches
Trichogramma <- scan(
  file = "./Input/OrthoMCL/Amelia/rewasporthology/no_longer_Tpre_lineage_specific.txt",
  what = "character"
)
class(Trichogramma)

# Find which groups have extra Trichogramma matches
BLAST <- read.table(file = "./Input/OrthoMCL/Amelia/rewasporthology/tpre_lineage_vs_HYMN.6.txt")
head(BLAST)
BLAST <- BLAST [,c(1,2)]
names(BLAST) <- c("TPRE_ID", "gene_ID")
BLAST <- unique(BLAST)

# Annotate relatives of Trichogramma divergent genes
Divergent <- BLAST[which(BLAST$TPRE_ID%in%Trichogramma),] # Some TPRE genes match to more than 1 ortholog
Divergent$gene_ID <- ifelse(
  grepl(pattern = "XP_", x = Divergent$gene_ID),
  str_extract(string = Divergent$gene_ID,
    pattern = "XP_[A-z, 0-9, .]*"),
  as.character(Divergent$gene_ID)
)
Divergent$gene_ID <- ifelse(
  grepl(pattern = "NP_", x = Divergent$gene_ID),
  str_extract(
    string = Divergent$gene_ID,
    pattern = "NP_[A-z, 0-9, .]*"),
  as.character(Divergent$gene_ID)
)

# Search for the relatives of Divergent genes in OG members
# OG_Trichogramma = rep(NA, length(Divergent))
# names(OG_Trichogramma) = Divergent
Divergent$OG_ID  = NA
for (i in Divergent$gene_ID) {
  pos <- grep(pattern = i, x = OGs$Data)
  # print(paste0(c("Searching for gene ", i, "Found at ", pos)))
  ifelse(
    pos!=0,
    Divergent[which(Divergent$gene_ID==i), "OG_ID"] <- as.character(OGs$OG_ID[pos]),
    next)
}
head(Divergent)
DivergentOGs <- as.character(na.exclude(unique(Divergent$OG_ID)))

# Concatenate matches to respective entries of OrthoMCL lists
for (i in which(OGs$OG_ID%in%DivergentOGs)) { # Actual working loop region
# for (i in which(OGs$OG_ID%in%OGs$OG_ID)) { # Testing line, includes TPRE match in all OGs
  OGs$Data[i] <- paste(as.character(OGs$Data[i]), "TPRE", sep = "|")
}

####

# Create list of species mapped in each orthologous group
Taxa = str_extract_all(string = OGs$Data, pattern = "[A-Z]{4}")
names(Taxa) = OGs$OG_ID
Taxa = lapply(X = Taxa, FUN = unique)
unique(unlist(Taxa))

# Subset only OGs with Nasonia matches
Taxa = Taxa[which(sapply(Taxa, function(x) {any(Nasonia%in%x)}))]

# Annotate taxonomic match
Stratum_Chalcid = sapply(Taxa, function(x){any(Chalcid%in%x)})
Stratum_Apocrita = sapply(Taxa, function(x){any(Apocrita%in%x)})
Stratum_Hymenoptera = sapply(Taxa, function(x){any(Hymenoptera%in%x)})

summary(cbind(Stratum_Hymenoptera, Stratum_Apocrita, Stratum_Chalcid))

# Annotate furthest taxonomic match for each OG
Strata = rep("Nasonia", length(Taxa))
  
for (i in 1:length(Strata)){
  if (Stratum_Hymenoptera[i]) {
    Strata[i] = "Hymenoptera"
  } else if (Stratum_Apocrita[i]) {
    Strata[i] = "Apocrita"
  } else if (Stratum_Chalcid[i]) {
    Strata[i] = "Chalcid"
  } else {
    next
  }
}

Strata = data.frame(
  OG_ID = names(Taxa),
  Stratum = Strata
)

# Match each Nasonia gene with respective OG stratum
NasStrata = merge(NasGenes, Strata, by = "OG_ID", all.x = T, all.y = F)
summary(NasStrata)
write.csv(x = NasStrata, file = file.path(newdir, "NasStrata.csv"))

#### Integrate Muscidifurax data (Pteromalid level)

Mus_OGs = read.csv(file = "Input/OrthoMCL/Yogeshwar/mclOutput_3May2016.groups.Nasvionly.txt", header = F, sep = ":")
names(Mus_OGs) = c("OG_ID", "Data")

# Create list of orthologous assignment per each nasonia gene
NasGenes = str_extract_all(string = Mus_OGs$Data, pattern = "Nasvi[A-z,0-9]*")
names(NasGenes) = Mus_OGs$OG_ID
NasGenes = stack(NasGenes)
names(NasGenes) = c("geneID", "OG_ID")

# Create list of species mapped in each orthologous group
Taxa = str_extract_all(string = Mus_OGs$Data, pattern = "[A-Z]{4}(?=[|])")
names(Taxa) = Mus_OGs$OG_ID
Taxa = lapply(X = Taxa, FUN = unique)
unique(unlist(Taxa))

# Subset only OGs with Nasonia matches
Nasonia = "NVIT"
Taxa = Taxa[which(sapply(Taxa, function(x) {Nasonia%in%x}))]

# Annotate OGs with Pteromalid matches
Pteromalid = c("MURA", "MUNI", "MELL")
Stratum_Pteromalid = sapply(Taxa, function(x){any(Pteromalid%in%x)})

NasPteromalid = stack(Stratum_Pteromalid)
names(NasPteromalid) = c("Pteromalid", "OG_ID")

NasPteromalid = merge(NasPteromalid, NasGenes, by = "OG_ID", all.x = T, all.y = F)
NasPteromalid = NasPteromalid[,c("geneID", "Pteromalid")]

# Merge with older strata
NasStrata = NasStrata[,c("geneID", "Stratum")]
NasStrata = merge(NasStrata, NasPteromalid)
NasStrata$Stratum = as.character(NasStrata$Stratum)
NasStrata$Stratum = ifelse(NasStrata$Stratum=="Nasonia"&NasStrata$Pteromalid==T, "Pteromalid", NasStrata$Stratum)

NasStrata = NasStrata[,c("geneID", "Stratum")]
NasStrata$Stratum = factor(x = NasStrata$Stratum, levels = c("Hymenoptera", "Apocrita", "Chalcid", "Pteromalid", "Nasonia"))
summary(NasStrata)
write.csv(x = NasStrata, file = file.path(newdir, "NasStrata_pteromalid.csv"))


################### Outdated code snippets
#######All genes in this file are ortholog at the hymenoptera level
# Trichogramma <- scan(
#   file = "./Input/OrthoMCL/Amelia/rewasporthology/navi_genes_with_signif_asex_hit.txt", 
#   what = "character"
# )
# 
# table(NasStrata[which(NasStrata$geneID%in%Trichogramma),"Stratum"])
# # NasStrata$Stratum <- 
#  ifelse(NasStrata$geneID%in%Trichogramma&NasStrata$Stratum%in%c("Nasonia", "Pteromalid"), "SSSS", as.character(NasStrata$Stratum))


# Annotate which genes have older matches in Saxton
tax_depth <- read.delim("./Input/DatasetS1.tsv")
tax_depth <- tax_depth[,c("ogs2", "strata")]
names(tax_depth)<-c("geneID","strata")
tax_depth$strata<-factor(tax_depth$strata, levels = c("Metazoa", "Arthropod", "Insect", "Hymenoptera", "Wasp"))

tax_depth_2 <- merge(tax_depth, NasStrata, by = "geneID") 
head(tax_depth_2)
tax_depth_2$Full_strata <- ifelse(tax_depth_2$strata%in%c("Metazoa", "Arthropod", "Insect"), 
  as.character(tax_depth_2$strata),
  as.character(tax_depth_2$Stratum)
)
tax_depth_2$Full_strata <- factor(
  x = tax_depth_2$Full_strata, 
  levels = c("Metazoa", "Arthropod", "Insect", "Hymenoptera", "Apocrita", "Chalcid", "Pteromalid", "Nasonia"))
table(tax_depth_2$Full_strata, useNA = 'ifany')
