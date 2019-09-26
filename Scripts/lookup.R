## Lookup script to identify different clusters 
library(plyr)
library(stringr)
rm(list=ls())

# Load cluster assignments
clusterassignments_simple <- read.csv("~/My Box Files/Data/Nasonia_development/Output/WGCNA_clustering_biweight/clusterassignments_simple.csv")
# Load genome information
NVIT_OGS2_goodannotcomplete <- read.csv("~/My Box Files/Data/Nasonia_development/Input/NVIT_OGS2_goodannotcomplete.csv")
# select columns of interest
# For cluster assignments
clusterassignments_simple<-clusterassignments_simple[,c("clusterID","transcriptID")]
clusterassignments_simple$geneID<-str_extract(clusterassignments_simple$transcriptID, pattern = "[[:alnum:]]*")
clusterassignments_simple$spliced<-str_extract(clusterassignments_simple$transcriptID, pattern = "[[:alnum:]]*$")
# Calculate cluster size
clusterassignments_simple<-merge(clusterassignments_simple, as.data.frame(table(clusterassignments_simple$clusterID)), by.x = "clusterID", by.y = "Var1", all.x = T)
names(clusterassignments_simple)<-c("clusterID","transcriptID","geneID","spliced","clustersize")

# For genome data
OGS2<-NVIT_OGS2_goodannotcomplete[,c("geneID","Name","groupname1","ODB6_OG_ID","location","Transposon","multicopy","OG_Copynumber","Nresiduals","adult_female_meth_status")]

# Merge all data
metatable<-merge(clusterassignments_simple, OGS2, by="geneID", all.x=T, all.y=F)

## lookup by gene name
metatable[grep("histone", metatable$Name),]

# Add matrix of ternary integration and parcellation



## Test whether genes within the same cluster share similar protein evolutionary rates, compared with genes within the same OG and interaction (diversification)
library(lme4)
library(MuMIn)
evorates<-na.exclude(metatable[,c("Nresiduals","clusterID","ODB6_OG_ID")])
evorates_lmfull<-lmer(data=evorates, formula = Nresiduals~1+(1|clusterID)*(1|ODB6_OG_ID))

plot(evorates_lmfull)
summary(evorates_lmfull)

evorates_lm0<-lm(data=evorates, formula = Nresiduals~1)
evorates_lmOG<-lmer(data=evorates, formula = Nresiduals~1+(1|ODB6_OG_ID))
evorates_lmClus<-lmer(data=evorates, formula = Nresiduals~1+(1|clusterID))

model.sel(list(evorates_lm0, evorates_lmClus, evorates_lmOG, evorates_lmfull))

######### outdated code
# # For eigenexons
# eigenexons_assignments2<-eigenexons_assignments[,c("geneID","eigenexonID","splicing_cathegory","constitutive")]
# eigenexons_assignments2$splicing_cathegory<-as.factor(ifelse(eigenexons_assignments2$splicing_cathegory=="spliced"&eigenexons_assignments2$constitutive=="constitutive","constitutive",as.character(eigenexons_assignments2$splicing_cathegory)))
# eigenexons_assignments2<-eigenexons_assignments2[,c("geneID","eigenexonID","splicing_cathegory")]

# # How many exons of each type are assigned to a cluster?
# prop.table(table(eigenexons_assignments2$eigenexonID%in%clusterassignments_simple$eigenexonID, eigenexons_assignments2$splicing_cathegory), margin=2)
