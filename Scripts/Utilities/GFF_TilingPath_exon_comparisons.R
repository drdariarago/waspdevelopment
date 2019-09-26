### GFF vs Tiling Path comparisons


library(stringr)
library(reshape2)
library(lattice)
library(ggplot2)

rm(list=ls())

### Load datasets

# Load exon expression data 
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
nasoniaE<-nasonia_devtesova.exon.r99[,c(1,grep("(adult|ovaries|testes).*[1-3]$", names(nasonia_devtesova.exon.r99)))]
rownames(nasoniaE)<-nasoniaE$EXON
nasoniaE<-as.matrix(nasoniaE[,-1])
nasoniaE<-as.data.frame(apply(nasoniaE, c(1,2), function(x){ifelse(x>0,x,0)})) # remove noise
nasoniaE$transcriptID<-str_extract(row.names(nasoniaE), pattern = "Nasvi[[:alnum:]]*")

# Load GFF genome annotation
OGS2_GFF <- read.table("./Input/nvit2_evigenes_pub11u.good.gff", header=FALSE, comment.char = "#", sep = "\t")
OGS2_GFF <- droplevels(OGS2_GFF[which(OGS2_GFF$V3=="exon"),])
OGS2_GFF$transcriptID<-str_extract(string = OGS2_GFF$V9, pattern = "Nasvi[[:alnum:]]*")

# Do both sets contain the same set of transcripts?
table(unique(nasoniaE$transcriptID)%in%unique(OGS2_GFF$transcriptID)) # 9446 transcripts are missing from the tiling path dataset
table(unique(OGS2_GFF$transcriptID)%in%unique(nasoniaE$transcriptID)) # 1159 transcripts are missing from the GFF dataset

# Filter only shared transcripts and sort them by transcript/exon ID (expression data) and transcript/start position (GFF)
nasoniaE<-nasoniaE[which(nasoniaE$transcriptID%in%OGS2_GFF$transcriptID),]
nasoniaE<-nasoniaE[order(row.names(nasoniaE)), ]
OGS2_GFF<-OGS2_GFF[which(OGS2_GFF$transcriptID%in%nasoniaE$transcriptID),]
OGS2_GFF<-OGS2_GFF[order(OGS2_GFF$transcriptID, OGS2_GFF$V4),]

# Do transcripts in tiling path and GFF files have the same number of total exons?
diffexons<-merge(as.data.frame(table(nasoniaE$transcriptID)),as.data.frame(table(OGS2_GFF$transcriptID)), by = "Var1")
diffexons<-diffexons$Freq.x-diffexons$Freq.y 
table(diffexons) # Some models differ in number of reported exons, GFF ones always have more
table(diffexons==0) # 2316 transcripts have different number of reported exons
prop.table(table(diffexons==0)) # 11% of the total shared transcripts in the two datasets

# Do transcripts in tiling path and GFF files have the same number of unique exons?
nasoniaE$exonID<-str_extract(string = row.names(nasoniaE), pattern = "[[:digit:]]*$")
nasoniaE$geneID<-str_extract(string = row.names(nasoniaE), pattern = "[^t]*")
nasoniaE$exonID<-paste(nasoniaE$geneID, nasoniaE$exonID, sep="e")
uniquenasoniaE<-nasoniaE[which(duplicated(nasoniaE$exonID)==F),]

OGS2_GFF$exonID<-paste(
  "s",str_extract(string = OGS2_GFF$V1, pattern = "[[:digit:]]"),
  "b",OGS2_GFF$V4, "e", OGS2_GFF$V5,OGS2_GFF$V7, sep=""
)
OGS2_GFF$geneID<-str_extract(string = OGS2_GFF$transcriptID, pattern = "[^t]*")
uniqueOGS2_GFF<-OGS2_GFF[which(duplicated(OGS2_GFF$exonID)==F),]

diffuniquexons<-merge(as.data.frame(table(uniquenasoniaE$geneID)),as.data.frame(table(uniqueOGS2_GFF$geneID)), by = "Var1")
diffuniquexons<-diffuniquexons$Freq.x-diffuniquexons$Freq.y
table(diffuniquexons) # both sets differ in numbers of unique exons
table(sign(diffuniquexons)) # most genes 
prop.table(table(sign(diffuniquexons)))
