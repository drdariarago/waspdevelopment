## Check where the venom proteins are
library(stringr)
library(plyr)

newdir<-file.path(getwd(), "Output/venoms")
dir.create(newdir)


jackvenoms <- read.csv(file = "./Input/nVit_venoms-jack.csv")
jackvenoms<-jackvenoms[grep("Nasvi", jackvenoms[,1]),]

CCRE_assignments <- read.csv(file = "./Output/CCREs/CCRE_assignments.csv")[,-1]
CCRE_nasdevgeneexon <- read.csv(file = "./Output/CCREs/CCRE_nasdevgeneexon.csv")
row.names(CCRE_nasdevgeneexon) <- CCRE_nasdevgeneexon$X
CCRE_nasdevgeneexon <- CCRE_nasdevgeneexon[,-1]
eigenexon_assignments <- read.csv(file = "./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexons_assignments.csv")
nasdevexon_prefiltering <- read.csv(file = "./Output/exon_to_transcript_clustering_rankingonly_adjusted/nasdevexon_prefiltering.csv")[,-1]
row.names(nasdevexon_prefiltering) <- nasdevexon_prefiltering$ExonID
nasdevexon_postfiltering <- read.csv(file = "./Output/exon_to_transcript_clustering_rankingonly_adjusted/nasdevexon_postfiltering.csv")[,-1]
row.names(nasdevexon_postfiltering) <- nasdevexon_postfiltering$ExonID
nodedata <- read.csv(file = "./Output/Results_compiler_CCRE/nodedata_full.csv")[,-1]

# nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/2015_05_19_Justin/nasonia_devtesova.exon.r99.score2")
# row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON
# nasdevexon<-nasonia_devtesova.exon.r99[,grep("[.][123]",names(nasonia_devtesova.exon.r99))]
# nasdevexon<-nasdevexon[,-grep("(ovaries|testes)", colnames(nasdevexon))]
# nasdevexon$geneID<-str_extract(row.names(nasdevexon), "^[^t]*") # including the t postscript to disambiguate *b and non b gene names
# nasdevexon$exonID<-str_extract(row.names(nasdevexon), "[[:digit:]]*$")
# nasdevexon$ExonID<-paste(nasdevexon$geneID, nasdevexon$exonID, sep = "e")
# nasdevexon<-nasdevexon[which(duplicated(nasdevexon$ExonID)==F),]

## Check which genes are present in CCREs or as nodes
summary(jackvenoms[,1]%in%nasdevexon_prefiltering$geneID)
# 101 are present in raw dataset, 1 is missing (Nasvi2EG036525)
summary(jackvenoms[,1]%in%nasdevexon_postfiltering$geneID)
# 101 are present after sample expression filtering
summary(jackvenoms[,1]%in%str_extract(CCRE_assignments$Row.names, "^[^_]*"))
# 93 are present as CCREs
summary(jackvenoms[,1]%in%str_extract(row.names(CCRE_nasdevgeneexon), "^[^_]*"))
# 29 are present as individual nodes
summary(jackvenoms[,1]%in%nodedata$geneID)

# Annotate missing venoms
novenom <- which(jackvenoms[,1]%in%str_extract(CCRE_assignments$Row.names, "^[^_]*")==F&
                   jackvenoms[,1]%in%str_extract(row.names(CCRE_nasdevgeneexon), "^[^_]*")==F)
novenomIDs <- jackvenoms[novenom,1]
novenomIDs

# The only gene that is really missing is absent from the exon probes as well
novenom <- which(jackvenoms[,1]%in%str_extract(CCRE_assignments$Row.names, "^[^_]*")==F&
                   jackvenoms[,1]%in%nodedata$geneID==F)
novenomIDs <- jackvenoms[novenom,1]
novenomIDs

## Save list of CCREs or nodes containing at least one venom gene
venoms <- unique(CCRE_assignments[which(jackvenoms[,1]%in%str_extract(CCRE_assignments$Row.names, "^[^_]*")),"CCREclust"])
venoms <- nodedata[which(nodedata$geneID%in%jackvenoms[,1]|nodedata$geneID%in%venoms),]
## Save list of nodes containing at least one venom gene (or CCRE that contains a venom gene)

write.csv(x = venoms, "./Output/venoms/venoms.csv")

# # Add list of their raw exon expression values
# novenomlist <- lapply(X = novenomIDs, FUN = function(x){
#   nasdevexon[grep(x, nasdevexon$geneID),]
# })
# names(novenomlist)<-novenomIDs
# 
# novenomIDs[which(novenomIDs%in%eigenexon_assignments$geneID==F)]
# ## Nasvi2EG036525 is missing after exon collapsing, the remainder is still in
# nasdevexon[grep("Nasvi2EG036525", nasdevexon$geneID),] #  absent also from array file...
# 
# # Restrict expression set to venom exons present in raw data but not in processed data
# expr_novenomIDs<-novenomIDs[which(novenomIDs%in%eigenexon_assignments$geneID==T)]
# nasdevexon <- nasdevexon[which(nasdevexon$geneID%in%expr_novenomIDs),]
