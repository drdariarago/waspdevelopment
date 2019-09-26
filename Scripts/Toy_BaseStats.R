## this functions converts the logarithm of odds (lods) value from eBayes into probability of gene being DE
lodprobs<-function(y){
  x<-exp(y)
  x/(x+1)
}

#### generate toy dataset with a random selection of 1000 genes
# data import and formatting
nasonia_devtesova.gene.r99 <- read.delim("~/My Box Files/Data/Nasonia_development/Nasonia Score files/nasonia_devtesova.gene.r99.score")
row.names(nasonia_devtesova.gene.r99)<-nasonia_devtesova.gene.r99$GENE
# select 1000 random rows
nasonia_devtesova.gene.r99<-nasonia_devtesova.gene.r99[sample(1:nrow(nasonia_devtesova.gene.r99),10000, replace=F),]
# filter only rows with samples (exclude means)
nasoniadevgene<-nasonia_devtesova.gene.r99[,grep("[.][123]",names(nasonia_devtesova.gene.r99))]
# exclude gonads for now
nasoniadevgene<-nasoniadevgene[,-grep("(ovaries|testes)", names(nasoniadevgene))]
# and remove main non-toy dataset
rm(nasonia_devtesova.gene.r99)



####WGCNA based quality control
library(WGCNA)
allowWGCNAThreads()
library(ecodist)
library(ape)
# check genes with missing values or zero variance
goodSamplesGenes(t(nasoniadevgene), verbose=3)
# cluster samples to detect anomalies
h1<-flashClust(dist(t(nasoniadevgene)))
plot(as.phylo(h1), tip.color=grepl("female", h1$labels)+1, type="cladogram")


#### LIMMA data formatting 
library(limma)
## convert into EList object (list first, then create as EList)
nasoniadevlist<-list(E=as.matrix(nasoniadevgene), genes=as.character(row.names(nasoniadevgene)), targets=names(nasoniadevgene))
nasdevEList<-new("EList", nasoniadevlist)
rm(list=c("nasoniadevlist","nasoniadevgene"))


##### LIMMA modelling
## create (toy) unidimensional single factor design matrix
design<-model.matrix(~0+factor((grepl("female",nasdevEList$targets))))
colnames(design)<-c("female","male")
## run linear models
fit<-lmFit(nasdevEList, design)
## compute contrasts
contrast.matrix<-makeContrasts(female-male, levels=design)
fit2<-contrasts.fit(fit, contrast.matrix)
fit2<-eBayes(fit2)
## plot as venn diagram
vennDiagram(decideTests(fit2))

## create toy unidimensional multi-factor design matrix
design<-model.matrix(~0+factor(as.character(regmatches(nasdevEList$targets, gregexpr("^[[:alnum:]]+",nasdevEList$targets)))))
colnames(design)<-levels(factor(as.character(regmatches(nasdevEList$targets, gregexpr("^[[:alnum:]]+",nasdevEList$targets)))))
contrast.matrix<-makeContrasts(testes-ovaries, adult-ovaries, emb10-ovaries, levels=design)
# fit models
fit<-lmFit(nasdevEList, design)
fit2<-contrasts.fit(fit, contrast.matrix)
fit2<-eBayes(fit2)
# display results as table and venn diagram
topTable(fit2, coef=1)
vennDiagram(decideTests(fit2))

## create toy unidimensional multi-factor all-vs-all design matrix
design<-model.matrix(~0+factor(as.character(regmatches(nasdevEList$targets, gregexpr("^[[:alnum:]]+",nasdevEList$targets)))))
colnames(design)<-levels(factor(as.character(regmatches(nasdevEList$targets, gregexpr("^[[:alnum:]]+",nasdevEList$targets)))))

library(plyr) #compute all combinations of factors and paste for contrasts
allcontrasts<-laply(combn(levels(factor(as.character(regmatches(nasdevEList$targets, gregexpr("^[[:alnum:]]+",nasdevEList$targets))))),2, simplify=F), .fun=function(x)paste(x, collapse="-"))
contrast.matrix<-makeContrasts(contrasts=allcontrasts, levels=design)

# fit models
fit<-lmFit(nasdevEList, design)
fit2<-contrasts.fit(fit, contrast.matrix)
fit2<-eBayes(fit2)
# display results as table and venn diagram
topTable(fit2, coef=1)
vennDiagram(decideTests(fit2)) # not possible because of excessive factor levels
# the most DE gene in adults vs embryos is glucose dehidrogenase (check Bart, remember it's a subset)


## create toy bidimensional multi-factor all-vs-all design matrix
stage<-factor(as.character(regmatches(nasdevEList$targets, gregexpr("^[[:alnum:]]+",nasdevEList$targets))))
male<-factor(grepl("[^[:alnum:]]male|testes", nasdevEList$targets))
design<-model.matrix(~0+stage*male)

# fit models
fit<-lmFit(nasdevEList, design)
fit2<-eBayes(fit)

# display results as table and venn diagram
topTable(fit2, coef=1)
vennDiagram(decideTests(fit2)) # not possible because of excessive factor levels
