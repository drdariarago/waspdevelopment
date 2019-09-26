## Compile results from the different analyses
library(reshape2)
library(plyr)
library(stringr)
rm(list=ls())
newdir<-file.path(getwd(),"Output/Results_compiler")
dir.create(newdir)

## Load OGS2 genome
NVIT_OGS2_goodannotcomplete <- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")

## Load eigenexon assignments
eigenexons_assignments <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexons_assignments.csv")
  
# ## Load eigenexon E values
eigenexon_evalues <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexon_evalues.csv")

## Load cluster assignments
clusterassignments_full <- read.csv("./Output/WGCNA_clustering_biweight/clusterassignments_full.csv")
names(clusterassignments_full)[1]<-"eigenexonID"

# ## Load Taxonomic profile
# OGS2_Taxascan <- read.table("./Output/Taxonomic_enrichment_BEAR/OGS2_Taxascan.csv", header=TRUE, quote="\"")

## Merge datasets
splicing_cat<-unique(eigenexons_assignments[,c("geneID","splicing_cathegory","clusters")])
splicing_cat<-ddply(splicing_cat, .(geneID), summarize, splicing_cathegory=splicing_cathegory[1], clusters=max(clusters))

NVIT_OGS2<-merge(NVIT_OGS2_goodannotcomplete, splicing_cat, by="geneID", all.y=T)
NVIT_OGS2<-merge(NVIT_OGS2, clusterassignments_full[,c("eigenexonID","clusterID")])
# NVIT_OGS2<-merge(NVIT_OGS2, OGS2_Taxascan, by.x="geneID", by.y="Gene_ID", all.x=T)

clusterassignments_short<-clusterassignments_full[,c("eigenexonID", "clusterID")]
clusterassignments_short$geneID<-str_extract(string = clusterassignments_short$eigenexonID, pattern = "^[^_]*")

clusterdata<-merge(clusterassignments_short, NVIT_OGS2[,c("geneID","ODB6_OG_ID","quality7","IntMatch","IntPred","OG_Copynumber","ratio","Nresiduals","splicing_cathegory","clusters","adult_female_meth_status")], by="geneID", all.y=T)

## Summarize by module
moduledata<-ddply(clusterdata, .(clusterID), summarize, .progress = "text",
                  nNodes=length(unique(eigenexonID)),
                  nGenes=length(unique(geneID)),
                  nSplicing=sum(grepl("con", eigenexonID)),
                  nTrans=sum(grepl("fac", eigenexonID)),
                  nParalogs_OGS=sum(grepl("Paralog",quality7)),
                  nOrthologs_OGS=sum(grepl("Ortholog",quality7)),
                  nParalogs_ODB=sum(is.na(OG_Copynumber)==T),
                  nOrthologs_ODB=sum(is.na(OG_Copynumber)==F),
                  nMeth=sum(adult_female_meth_status=="Methylated", na.rm = T),
                  nUnmeth=sum(adult_female_meth_status=="Unmethylated", na.rm = T)
)

prop_moduledata<-moduledata[,-c(1,2)]/moduledata[,2]

moduledata<-merge(moduledata, prop_moduledata, by=clusterID)
## Load stage-specific sexbias coefficients
fdr_perm_dWithin <- read.csv("./Output/sex_modularity_GLMs/fdr_perm_dWithin.csv")
fdr_perm_dOut <- read.csv("./Output/sex_modularity_GLMs/fdr_perm_dOut.csv")

# Select only relevant factors 


# Recode as cluster~factor table
# dWithin_avg<-recast(data = dWithin_avg, formula = clusterID~factorName)
# dOut_avg<-recast(data = dOut_avg, formula = clusterID~factorName)


## Merge datasets
eigenexon_data<-merge(clusterassignments_full[,c("eigenexonID", "clusterID"),],eigenexons_assignments[c("eigenexonID","geneID","constitutive")],all.x=T, all.y=T)
eigenexon_data<-merge(eigenexon_data, dOut_avg, by="clusterID", all.x=T, all.y=F)
eigenexon_data<-merge(eigenexon_data, dWithin_avg, by="clusterID", all.x=T, all.y=F, suffixes = c(".dOut",".dWithin"))
eigenexon_data<-merge(eigenexon_data, OGS2_Taxascan, by.x="geneID", by.y="Gene_ID", all.x=T, all.y=F)
eigenexon_data$Taxascan_max<-ifelse(is.na(eigenexon_data$Taxascan_max),0,eigenexon_data$Taxascan_max)

# Save compiled dataset
write.csv(eigenexon_data, file=file.path(newdir, "eigenexon_data.csv"))