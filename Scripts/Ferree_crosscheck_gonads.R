## Add data on testes/ovary specific expression
# Initialize script
library(plyr)
library(stringr)
newdir <- file.path(getwd(), "./Output/Ferree_gonads")
dir.create(newdir)

# Load main dataset
transcriptdata <- read.csv(file = "./Output/Results_compiler_CCRE/transcriptdata_full.csv")
genedata <- transcriptdata[which(transcriptdata$event_type=="transcription"),]
# Load OGS2 to NCBI mapping table
conversion <- read.csv("./Input/12864_2016_2886_MOESM6_ESM.csv")
# Load extra data from paper
germline <- read.csv(file = "./Input/Ferree2015_SM/Ferree_TableS10_germline.csv")
testes <- read.csv(file = "./Input/Ferree2015_SM/Ferree_TableS8_testes.csv")
ovaries <- read.csv(file = "./Input/Ferree2015_SM/Ferree_TableS9_ovaries.csv")
meiosis <- read.csv(file = "./Input/Ferree2015_SM/Ferree_TableS11_meiosis.csv")
centrosomal <- read.csv(file = "./Input/Ferree2015_SM/Ferree_TableS12_centrosomal.csv")

# Subset gene info
genedata <- genedata[,c("clusterID","geneID","nodeID","devsexbias")]
# Add testes info
head(testes[1:4])
genedata$testes <- genedata$geneID%in%testes$Nvit2.Gene.ID
# Add ovaries info
head(ovaries[,1:4])
ovaries <- merge(ovaries[,1:4], conversion[,c(1,6)], by.x = "Transcript.ID", by.y = "OGS1geneID", all.x = T, all.y = F)
ovaries$OGS2GeneID <- str_extract(string = ovaries$OGS2transcriptID, pattern = "^[^t]*")
genedata$ovaries <- genedata$geneID%in%ovaries$OGS2GeneID
# Add germline info
head(germline[,1:4])
genedata$germline <- genedata$geneID%in%germline$Nvit2.Gene.ID
# Add meiosis info
head(meiosis[,1:4])
genedata$meiosis <- genedata$geneID%in%meiosis$Nvit2.Gene.ID
# Add centrosomal info
head(centrosomal[,1:4])
genedata$centrosomal <- genedata$geneID%in%centrosomal$Nasonia.Ortholog

head(genedata)

## Save as CSV
write.csv(genedata, file = file.path(newdir, "genedata_gonads.csv"))


## Summarize by module
# library(reshape2)
# recast(data = genedata[,c("clusterID","testes","ovaries","germline","meiosis","centrosomal")], id.var = "clusterID",
#        formula = clusterID ~ testes + ovaries + germline + meiosis + centrosomal)

clusterdata <- ddply(.data = genedata[,c("clusterID","testes","ovaries","germline","meiosis","centrosomal")],
      .variables = .(clusterID), .fun = colwise(sum))
## Add module data
clusterdata2 <- read.csv("./Output/Results_compiler_CCRE/clusterdata_full.csv")[,-1]

clusterdata <- merge(clusterdata, 
                     clusterdata2[,c("clusterID", "nGenes", "excessDup","excessMet","cluster_devsexbias","cluster_stagebiased","DiffIntegrated")])

write.csv(clusterdata, file = file.path(newdir, "clusterdata_gonads.csv"))
