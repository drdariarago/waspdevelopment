# Functional annotation of features
library(stringr)

# Import GFF
nvit2_GFF <- read.table("./Input/nvit2_evigenes_pub11u.good.gff", header=FALSE, sep = "\t", comment.char = "#",
                        colClasses=c("character", "character", "character", "integer",  "integer","character", "character", "character", "character"))
colnames(nvit2_GFF)<-c("seqname", "source", "feature", "start", "end","score", "strand", "frame", "attributes")

# extract gene and exon IDs
nvit2_GFF$geneID<-str_extract(string = nvit2_GFF$attributes, pattern = "Nasvi2EG[^t]*")
# extract exon ID (must check matching in main dataset, could be as simple as number from 5' end)

# select only [geneID, exonID, CDS, chromosome, start/end and strand]
nvit2_GFF<-nvit2_GFF[which(nvit2_GFF$feature=="CDS"),c("seqname","start","end","strand")]

# Import chromosome sequences

# Extract sequences for each exon annotated with its main gene and exon group

# Export fasta file of exons