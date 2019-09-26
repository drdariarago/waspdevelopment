## genome activation summary

## when are exons activated for the first time
nasoniadevexon>0

# Load gene data
nasonia_devtesova.gene.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.gene.r99.score")
row.names(nasonia_devtesova.gene.r99)<-nasonia_devtesova.gene.r99$GENE

# filter only rows with samples (exclude means)
nasoniadevgene<-nasonia_devtesova.gene.r99[,grep("[.][123]",names(nasonia_devtesova.gene.r99))]
# Exclude gonads 
nasoniadevgene<-nasoniadevgene[,-grep("(ovaries|testes)", names(nasoniadevgene))]



# Annotate samples with male, stage and replicate
geneData<-nasoniadevgene
testData$exonID<-row.names(testData)
# Gene names are defined as the portion of the exon name that goes from the beginning to the last non "t" character
testData$geneID<-unlist(regmatches(testData$exonID, gregexpr("^[^t]*",testData$exonID)))
library(reshape)
testData<-melt(testData)
testData$Male<-grepl("_male",testData$variable)
testData$Stage<-factor(as.character(regmatches(testData$variable, gregexpr("^[[:alnum:]]+",testData$variable))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
testData$replicate<-factor(as.character(regmatches(testData$variable, gregexpr("[[:digit:]]$",testData$variable))))
testData<-subset(testData, select=-c(variable))
