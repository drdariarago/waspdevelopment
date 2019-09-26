## Exon reduction report

# Initialize script
rm(list=ls())
library(stringr)
library(ggplot2)
graphdir<-file.path(getwd(), "./Graphics/exon_to_transcript_clustering_report/")
dir.create(graphdir)

# Load dataset
load("./Output/exon_to_transcript_clustering_rankingonly_adjusted/workspace.R")

# How many genes are present in the total dataset?
length(unique(str_extract(string = nasonia_devtesova.exon.r99$EXON, pattern = "^[^t]*")))
# How many unique exons are present in total dataset?
length(unique(paste(str_extract(string = nasonia_devtesova.exon.r99$EXON, pattern = "^[^t]*"), str_extract(string = nasonia_devtesova.exon.r99$EXON, pattern = "[[:digit:]]*$"))))

# How many genes are present before clustering
length(unique(nasdevexon$geneID))
# How many unique exons are present before clustering (expressed in >3 samples)
length(unique(nasdevexon$ExonID))
# How may genes are present after clustering 
length(unique(eigenexons$geneID))
# How many exons are present after clustering 
length(unique(eigenexons$exonID))

# How many exons are retained after expression filtering?
length(unique(eigenexons$exonID))/length(unique(paste(str_extract(string = nasonia_devtesova.exon.r99$EXON, pattern = "^[^t]*"), str_extract(string = nasonia_devtesova.exon.r99$EXON, pattern = "[[:digit:]]*$"))))

# How many genes are retained after expression fitering?
length(unique(eigenexons$geneID))/length(unique(str_extract(string = nasonia_devtesova.exon.r99$EXON, pattern = "^[^t]*")))

# How many exons are retained of each class?
table(eigenexons$splicing_cathegory)
table(eigenexons$splicing_cathegory, eigenexons$constitutive)

# Which proportion of remaining exons is constitutive/facultative
table(eigenexons$constitutive)
prop.table(table(eigenexons$constitutive))

# How many nodes are present in the final dataset?
length(unique(eigenexons$eigenexonID))

# How many of those are constitutive/facultative?
table(grepl("_fac", eigenexons[-which(duplicated(eigenexons$eigenexonID)),]$eigenexonID))
prop.table(table(grepl("_fac", unique(eigenexons$eigenexonID))))

# How many nodes are present of each class?
table(grepl("_fac", eigenexons[-which(duplicated(eigenexons$eigenexonID)),]$eigenexonID), eigenexons[-which(duplicated(eigenexons$eigenexonID)),]$splicing_cathegory)


# How many genes have splicing nodes?
table(table(eigenexons$geneID, eigenexons$constitutive)[,2]>0)
prop.table(table(table(eigenexons$geneID, eigenexons$constitutive)[,2]>0))

# What's the distribution of facultative eigenexons amongst genes?
table(table(eigenexons$eigenexonID)-1)
prop.table(table(table(eigenexons$eigenexonID)-1))

# Plot distribution of facultative eigenexons amongst genes
fac_freqs<-data.frame(table(eigenexons$eigenexonID)-1)
# ggplot(data=fac_freqs, aes(x=Freq))+geom_bar()+theme_bw()+scale_y_log10()
ggplot(data=fac_freqs, aes(x=as.factor(Freq)))+geom_bar()+theme_bw()+scale_y_log10(breaks = c(1,5,10,25,50,100,250,500,1000,2000,4000,8000,16000,32000))+ggtitle("Distribution of facultative exon groups amongst genes\n")+xlab(label = "Number of facultative transcripts per gene")+ylab(label = "Log Counts")+theme(panel.grid.major.x = element_blank())

# What's the distribution of exons amongst eigenexons?
# Plot frequency of exons, split between con and fac
exon_freqs<-data.frame(table(eigenexons$eigenexonID))
exon_freqs$con<-grepl("_con",exon_freqs$Var1)
table(exon_freqs$Freq, exon_freqs$con)

ggplot(data = exon_freqs, aes(x=as.factor(Freq), fill=con))+geom_histogram(position="dodge")+theme_bw()+scale_y_log10(breaks = c(1,5,10,25,50,100,250,500,1000,2000,4000,8000,16000,32000))+ggtitle("Distribution of exons within transcripts\n")+xlab(label = "Number of exons per transcript")+ylab(label = "Log Counts")+theme(panel.grid.major.x = element_blank())+scale_fill_discrete(name = "Transcript type", labels=c("Facultative","Constitutive"))

# Plot exons per transcript vs exons in gene
exons_pergene<-data.frame(table(eigenexons$geneID))
names(exons_pergene)<-c("geneID","nExons")
exon_freqs$geneID<-str_extract(string = exon_freqs$Var1, pattern = "^[[:alnum:]]*")
exon_freqs<-merge(exon_freqs, exons_pergene, by="geneID", all.x=T)
exon_freqs$spliced<-exon_freqs$geneID%in%unique(exon_freqs[grepl("fac",exon_freqs$Var1),"geneID"])

ggplot(data = exon_freqs[which(exon_freqs$spliced==T),], aes(x=Freq/nExons, fill=con))+geom_histogram(position="dodge")+theme_bw()+ggtitle("Proportion of exons in each transcript\n")+xlab(label = "Proportion of exons per transcript")+ylab(label = "Log Counts")+theme(panel.grid.major.x = element_blank())+scale_fill_discrete(name = "Transcript type", labels=c("Facultative","Constitutive"))

# Do paralogs show more isoforms than single copy genes?
# Calculate number of eigenexons per gene
eigenexons_freqs<-data.frame(table(unique(eigenexons[,c("geneID","eigenexonID")])$geneID))
names(eigenexons_freqs)<-c("geneID","nEigenexons")
# load genome annotation
OGS2 <- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")
# merge with eigenexon annotation
eigenexons_freqs_OGS2<-merge(eigenexons_freqs, OGS2, by.x="geneID", by.y="geneID")
# plot
ggplot(data = eigenexons_freqs_OGS2, aes(x=as.factor(nEigenexons), fill=quality7))+geom_histogram(position="dodge")+theme_bw()+ggtitle("Number of nodes gene")+xlab(label = "Nodes per gene")+ylab(label = "Log Counts")+theme(panel.grid.major.x = element_blank())+scale_fill_discrete(name = "Homology class", labels=c("NA","Ortholog","Paralog"))+scale_y_log10(breaks = c(5,10,25,50,100,250,500,1000,1500,3000))

# plot per SMBE poster
pdf(file=file.path(graphdir,"nodes_per_gene.pdf"), width = 6, height = 4)
ggplot(data = eigenexons_freqs_OGS2, aes(x=as.factor(nEigenexons), fill=quality7))+geom_histogram(position="dodge")+theme_bw()+ggtitle("Number of nodes per gene\n")+xlab(label = "Nodes per gene")+ylab(label = "Log Counts")+theme(panel.grid.major.x = element_blank(), legend.position=c(0.8,0.8), axis.text.x = element_text(size=15))+scale_fill_discrete(name = "Homology class", labels=c("NA","Ortholog","Paralog"))+scale_y_log10(breaks = c(5,10,25,50,100,250,500,1000,1500,3000))+theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12))
dev.off()

# ggplot(data = eigenexons_freqs_OGS2, aes(x=as.factor(nEigenexons)))+geom_histogram(position="dodge")+theme_bw()+ggtitle("Number of detected eigenexons per gene")+xlab(label = "Eigenexons per gene")+ylab(label = "Log Counts")+theme(panel.grid.major.x = element_blank())+scale_y_log10()+facet_wrap(~quality7, ncol = 1)


# ##### Outdated code
# # Are those exons expressed in less than 3 samples?
# table(row.names(nasdevexon)%in%eigenexons$exonID,apply(nasdevexon[,1:30],1,function(x){sum(x>0)}))
# # YES
# summary(apply(eigenexons,2,as.factor))
# # Are there any exons expressed in less than 3 samples amongst spliced genes?
# exprExons<-as.data.frame(apply(nasdevexon[which(row.names(nasdevexon)%in%eigenexons$exonID),1:30],1,function(x){sum(x>0)}))
# exprExons<-as.data.frame(apply(nasdevexon[,1:30],1,function(x){sum(x>0)}))
# names(exprExons)<-"count"
# exprExons<-merge(exprExons, eigenexons, by.x="row.names", by.y="exonID", all.x=T)
# exprExons$splicing_cathegory<-as.factor(ifelse(is.na(exprExons$splicing_cathegory),"unexpressed_exon",exprExons$splicing_cathegory))
# table(exprExons$count, exprExons$splicing_cathegory)
# 
# # How many total exons are not expressed in our dataset? 
# table(exprExons$count>2)
# prop.table(table(exprExons$count>2))
# table(exprExons$count>2, exprExons$splicing_cathegory)
# 
# # How many genes do they belong to? How many final eigenexons? How many of those are constitutive and how many facultative? What's the distribution of facultatives vs constitutives?
# 
# 
# ### Eigenexon clustering report
# 
# # How many eigenexons are present? How many genes do they belong to?
# 
# # How many eigenexons are constitutive vs facultative?
# 
# # How many genes have at least 1 facultative eigenexon? What is the distribution of facultative eigenexons amongst genes with fac eigenexons?
