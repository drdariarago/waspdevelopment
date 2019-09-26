## Diagnostics of limma naming in exon dataset

# Script for detecting which exons have statistically significant alternative splicing across sex/development
# Main script identical to BaseStats applied to exon database, but with 1 modification:
# Load Basestats data, subtract gene estimates from exon scores (flattened by geneID, stage and sex, merge and subtract, then unflatten)
# run limmma and FDR

# Load LIMMA results
load(file="./Output/BaseStats_fit2.RData")
# Load exon data
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON
# filter only rows with samples (exclude means)
nasoniadevexon<-nasonia_devtesova.exon.r99[,grep("[.][123]",names(nasonia_devtesova.exon.r99))]
# Exclude gonads 
nasoniadevexon<-nasoniadevexon[,-grep("(ovaries|testes)", names(nasoniadevexon))]
# Sava dataset before analyses
write.csv(nasoniadevexon, file="./Output/nasoniadevexon.csv")

#### Subtract estimates from samples

# Annotate samples with male, stage and replicate
testData<-nasoniadevexon
testData$exonID<-row.names(testData)
testData$geneID<-strtrim(row.names(testData),14)
library(reshape)
testData<-melt(testData)
testData$Male<-grepl("_male",testData$variable)
testData$Stage<-factor(as.character(regmatches(testData$variable, gregexpr("^[[:alnum:]]+",testData$variable))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
testData$replicate<-factor(as.character(regmatches(testData$variable, gregexpr("[[:digit:]]$",testData$variable))))
testData<-subset(testData, select=-c(variable))
# Load and annotate coefficients wit stage and male specific
testCoefs<-fit2$coefficients[row.names(fit2$coefficients)%in%testData$geneID,]
testCoefs<-melt(testCoefs)
testCoefs$Male<-grepl("male",testCoefs$X2)
testCoefs$Stage<-factor(as.character(substring(as.character(regmatches(testCoefs$X2, gregexpr("^[[:alnum:]]*",testCoefs$X2))),6)),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
testCoefs<-subset(testCoefs, select=-c(X2))
names(testCoefs)<-c("geneID","Coef","Male","Stage")

# Subtract male coefficients to males
testData2<-merge(testData[which(testData$Male==T),], testCoefs, by=c("geneID", "Stage","Male"))
testData2$value<-testData2$value-testData2$Coef
testData2<-testData2[,c("exonID","geneID","value","Male","Stage","replicate")]

# Subtract stage specific coefficients
testData3<-rbind(testData2, testData[which(testData$Male==F),])
testData3<-merge(testData3, testCoefs[which(testCoefs$Male==F),-which(colnames(testCoefs)=="Male")], by=c("geneID", "Stage"), all.x=T)
testData3$value<-testData3$value-testData3$Coef
testData3<-testData3[,c("exonID","geneID","Stage","Male","replicate","value")]

# Put into exon by sample format
testData3$Male<-as.factor(ifelse(testData3$Male==T, "male","female"))
nasoniadevexon2<-cast(testData3, formula=exonID~Stage+Male+replicate)
row.names(nasoniadevexon2)<-nasoniadevexon2$exonID
nasoniadevexon2<-subset(nasoniadevexon2, select=-exonID)

### WARNING: gene Nasvi2EG025210 (transcript Nasvi2EG025210bt1.1) is absent from the tiling path data, it is therefore excluded from our analyses here
nasoniadevexon2<-na.omit(nasoniadevexon2)

# Save as CSV and clean workspace
write.csv(nasoniadevexon2, file="./Output/Nasonia_Dev_Exon_Corrected.csv")
# Cleaning disabled for diagnostics
# rm(list=(grep("test", ls(), value=T)))

## Create LIMMA object
library(limma)
# Reformat as list
nasoniadevlist_exon<-list(E=as.matrix(nasoniadevexon2), genes=as.character(row.names(nasoniadevexon2)), targets=names(nasoniadevexon2))
## BEAR does not add gene names and sample names to rows and column for some reason, adding them manually
row.names(nasoniadevlist_exon$E)<-nasoniadevlist_exon$genes
colnames(nasoniadevlist_exon$E)<-nasoniadevlist_exon$targets
# Converting to LIMMA EList format
nasdevEList_exon<-new("EList", nasoniadevlist_exon)
# disabled for diagnostics
# rm(list=c("nasoniadevlist_exon"))

# save workspace for diagnostics
save.image(file="Exon_LIMMA_diagnostics.RData")
