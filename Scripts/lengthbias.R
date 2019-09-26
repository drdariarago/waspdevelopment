### Test for signal bias in long vs short genes 

library(stringr)
library(reshape2)
library(lattice)
library(ggplot2)
library(plyr)

# Load exon expression data
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
# filter only rows with samples (exclude means)
row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON
nasdevexon<-nasonia_devtesova.exon.r99[,grep("[.][123]",names(nasonia_devtesova.exon.r99))]

# Load gene annotation data
OGS2 <- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")

# Load eigenexon assignments
eigenexons_assignments <- read.csv("./Output/exon_to_transcript_clustering_rankingonly/eigenexons_assignments.csv")

# Load exon annotation data
OGS2_GFF <- read.table("./Input/nvit2_evigenes_pub11u.good.gff", header=FALSE, comment.char = "#", sep = "\t")
OGS2_GFF <- droplevels(OGS2_GFF[which(OGS2_GFF$V3=="exon"),])
OGS2_GFF$transcriptID<-str_extract(string = OGS2_GFF$V9, pattern = "Nasvi[[:alnum:]]*")
# # Retain only one exon per gene (removes exons included mutiple times because of multiple transcripts)
# OGS2_GFF<-OGS2_GFF[which(duplicated(OGS2_GFF[,c("V4","V5")])==F),]

# Add unique exon ID (scaffold, start, stop and strand)
OGS2_GFF$ExonID<-paste(sep = "",
  "S",str_extract(string = OGS2_GFF$V1, pattern = "[[:digit:]]*$"),
  "b",OGS2_GFF$V4,
  "e",OGS2_GFF$V5,
  ifelse(OGS2_GFF$V7=="+","f","r")
  )
# Select only one exon per genomic span and orientation
OGS2_GFF<-OGS2_GFF[which(duplicated(OGS2_GFF$ExonID)==F),]

# Select only first time each exon is annotated
nasdevexon$geneID<-str_extract(row.names(nasdevexon), pattern = "Nasvi2EG[[:digit:]]*[^t]")
nasdevexon$exonID<-str_extract(row.names(nasdevexon), pattern = "[[:digit:]]*$")
nasdevexon$ExonID<-paste(nasdevexon$geneID, nasdevexon$exonID, sep = ".")

nasdevexon<-nasdevexon[which(duplicated(nasdevexon$ExonID)==F),]

## Match exons in GFF file and expression file, via probe assignments
# Import probe mapping table
nasonia.probe.exon <- read.delim("./Nasonia Score files/nasonia.probe.exon", header=FALSE)
# Filter only exons present in expression dataset
nasonia.probe.exon<-nasonia.probe.exon[which(nasonia.probe.exon$V2%in%row.names(nasdevexon)),]
# Annotate scaffold and start position
names(nasonia.probe.exon)[1]<-"exonID"
nasonia.probe.exon$scaffold<-str_extract(string = nasonia.probe.exon$V1, pattern = "SCAFFOLD[[:digit:]]*")
nasonia.probe.exon$start<-str_extract(string = nasonia.probe.exon$V1, pattern = "[[:digit:]]*$")
nasonia.probe.exon$start<-as.numeric(str_replace(string = head(nasonia.probe.exon$start), pattern = "^0*", replacement=""))

# Annotate expression data with respective GFF exon entries
# Filter only genes present in nasonia probes
which(str_extract(OGS2_GFF$transcriptID, "^[^t]*")%in%str_extract(nasonia.probe.exon$V2, "^[^t]*"))

# match scaffold, then look for exons comprised in gene span
which(OGS2_GFF[which(OGS2_GFF$V1%in%nasonia.probe.exon[100,]$scaffold),"V4"]>=nasonia.probe.exon[100,"start"]&&
        OGS2_GFF[which(OGS2_GFF$V1%in%nasonia.probe.exon[100,]$scaffold),"V5"]<=nasonia.probe.exon[100,"start"])

# Annotate relative position to the first exon (reversed for reverse genes)

# Annotate number of exons

nExons<-ddply(nasdevexon, .(geneID), summarize, nExons=length(exonID))

# Filter only genes with more than 1 exon
nasdevexon<-merge(nasdevexon, nExons, by="geneID", all.x=T, all.y=T)
nasdevexon<-nasdevexon[which(nasdevexon$nExons>1),]

# Measure expression compared to 1st exon
nasdevexon<-nasdevexon[order(nasdevexon$geneID, nasdevexon$exonID),]
diffirst<-function(x){x-x[1]}
nasdiffexon<-ddply(nasdevexon, .(geneID), numcolwise(diffirst))
nasdiffexon<-cbind(nasdiffexon[,-38], nasdevexon[,c("exonID", "ExonID", "nExons")])

# Annotate gene length

# Melt dataset
nasdiffexon<-melt(nasdiffexon, id.vars = c("geneID", "exonID", "ExonID", "nExons"))

## Are exons closer to one end more expressed than to the other?

ggplot(nasdiffexon, aes(x=nExons, y=value))+geom_smooth()+theme_bw()
ggplot(nasdiffexon, aes(x=nExons, y=value, col=variable))+geom_point(alpha="0.3")+geom_smooth()+theme_bw()


## Are exons closer to one end more likely to be facultative than to the other?

#### Outdated code

# # Compare number of exons per transcript in GFF and exon files
# table(str_extract(row.names(nasdevexon), pattern = "[[:alnum:]]*")[1:100])
# table(OGS2_GFF[1:100,"transcriptID"])
# # Each exon is reported one time for each isoform in the final expression file
# 
# # Compare exon expression in different transcripts in the same sample
# testset<-nasdevexon[1:5000,]
# testset$geneID<-str_extract(row.names(testset), pattern = "Nasvi2EG[[:digit:]]*[^t]")
# testset$exonID<-str_extract(row.names(testset), pattern = "[[:digit:]]*$")
# testset$ExonID<-paste(testset$geneID, testset$exonID, sep = "_")
# 
# testset<-melt(testset, id.vars = c("geneID","exonID","ExonID"))
# # testset$value<-ifelse(testset$value<=0,NA,testset$value)
# testset<-ddply(testset, .(ExonID, variable), summarize, discordance=max(value)-min(value), median=median(value), mad=mad(value), .progress = "text")
# densityplot(testset$discordance)
# densityplot(testset$median)
# densityplot(testset$mad)
