# Test for enrichment of splicing node proximity and splicing-based GO terms
library(stringr)
library(WGCNA)
# method: start from ajdacency matrix, filter only rows containing genes.
# load(file = "./Output/WGCNA_clustering_biweight/adj")

## create dummy dataset
eigenexon_evalues <- read.csv("~/My Box Files/Data/Nasonia_development/Output/exon_to_transcript_clustering_rankingonly/eigenexon_evalues.csv")
row.names(eigenexon_evalues)<-eigenexon_evalues$eigenexonID
adj<-corFast(t(eigenexon_evalues[sample.int(n = nrow(eigenexon_evalues), size = 5000),-1]), use="pairwise.complete.obs")
adjG<-adj[grep("con",row.names(adj)),]
# create two tables, one with the splicing nodes as columns, one with the gene nodes as columns
adjGS<-adjG[,grep("fac", colnames(adjG))]
adjGG<-adjG[,grep("con", colnames(adjG))]
# for every row, calculate minimum and median distance and create data.frame with genes and their min/median distance to other genes, splicing nodes and the ratio between the two

SplicingDist<-data.frame(
  geneID=str_extract(row.names(adjG),pattern = "^[[:alnum:]]*"),
  GSdist=apply(adjGS, 1, median, na.rm=T),
  GGdist=apply(adjGG, 1, median, na.rm=T))
SplicingDist$GSratio<-SplicingDist$GSdist/SplicingDist$GGdist

# perform correlation analysis between distance/ratio and GO term enrichment
# load Nasonia gene IDs, then screen for orthoDB ids and all associated GO terms.
OrthoDB7_Arthropoda_tabtext <- read.delim("~/My Box Files/Data/Nasonia_development/Input/OrthoDB7_Arthropoda_tabtext", dec=",")
NVITR_ODB7<-OrthoDB7_Arthropoda_tabtext[which(OrthoDB7_Arthropoda_tabtext$Organism=="Nasonia vitripennis"),]
# Attach ODB ID to each gene, alongside with INTERPRO domains
SplicingDist<-merge(SplicingDist, NVITR_ODB7[,c("ODB7_OG_ID","Gene_ID","Protein_ID","InterPro_domains")], by.x="geneID", by.y="Gene_ID", all.x = F)
# Import and add interpro GO terms (importing as txt, then parsing relevant terms)
interpro2GO<-scan(file="./Input/interpro2go_release_49.0", what = character(), skip=3, sep="\n")
interpro2GO<-data.frame(
  InterPro_domain=str_extract(interpro2GO, pattern = "IPR[[:alnum:]]*"),
  GO_term=str_extract(interpro2GO, pattern = "GO:[^;]*"),
  GO_ID=str_extract(interpro2GO, pattern="GO:[[:digit:]]{7}"))
# Attach distances and gene ID to each domain
DomainDist<-str_extract_all(SplicingDist$InterPro_domains, pattern = "IPR[[:alnum:]]*")
DomainDist<-ldply(lapply(X = DomainDist, rbind))
DomainDist<-cbind(SplicingDist, DomainDist)
head(melt(DomainDist, id.vars = c("geneID","ODB7_OG_ID","GSdist","GGdist","GSratio","Protein_ID","InterPro_domains")))
# DomainDist<-sapply(interpro2GO$InterPro_domain, )
# perform correlation analyses