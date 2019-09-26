## Variance partitioning of gene dataset
## load raw data, analyze with LIMMA and add Fdr and fdr for developmental time series and gonads vs adults

# data import and formatting
nasonia_devtesova.gene.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.gene.r99.score")
row.names(nasonia_devtesova.gene.r99)<-nasonia_devtesova.gene.r99$GENE
# filter only rows with samples (exclude means)
nasoniadevgene<-nasonia_devtesova.gene.r99[,grep("[.][123]",names(nasonia_devtesova.gene.r99))]
# Exclude gonads 
nasoniadevgene<-nasoniadevgene[,-grep("(ovaries|testes)", names(nasoniadevgene))]

library(vegan)
# Create dissimilarity matrix of samples by gene
tnasoniadevgene<-t(nasoniadevgene)
dissnasoniadevgene<-dist(nasoniadevgene)

# create data.frame with descriptors of samples
datagene<-data.frame(stage=factor(as.character(regmatches(row.names(tnasoniadevgene), gregexpr("^[[:alnum:]]+",row.names(tnasoniadevgene)))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult")), male=factor(grepl("[^[:alnum:]]male|testes", row.names(tnasoniadevgene))))
datagene$stage_ordered<-ordered(datagene$stage, levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))

# run model (computationally intensive)
adonis1<-adonis(dissnasoniadevgene~stage/male, data=datagene, permutations=1000)
adonis2<-adonis(dissnasoniadevgene~stage_ordered/male, data=datagene, permutations=1000)

# Save data
save.image(file="./Output/Adonis.RData")