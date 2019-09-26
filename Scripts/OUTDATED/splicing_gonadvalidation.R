## Cross validation of splicing from gonads

## Initialize script
rm(list=ls())
library(plyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(vcd)
library(limma)
library(fdrtool)
source('./Scripts/multiplot.R', echo=F)
# initialize output path
newdir<-file.path(getwd(), "Output/splicing_gonadvalidation")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/splicing_gonadvalidation")
dir.create(graphdir)

## Load expression data
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/2015_05_19_Justin/nasonia_devtesova.exon.r99.score2")
nasoniaE<-nasonia_devtesova.exon.r99[,c(1,grep("(adult|ovaries|testes).*[1-3]$", names(nasonia_devtesova.exon.r99)))]
rownames(nasoniaE)<-nasoniaE$EXON
nasoniaE<-as.matrix(nasoniaE[,-1])
nasoniaE<-as.data.frame(apply(nasoniaE, c(1,2), function(x){ifelse(x>0,x,0)}))
rm(nasonia_devtesova.exon.r99)

## Load sex-specific genes
sexspec<- read.csv("./Output/sex_specific_nodes/sexspec_allstages.csv")[,-1]

## Load eigenexon assignments
eigenexons_assignments <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexons_assignments.csv")
eigenexons_assignments <- eigenexons_assignments[,c("exonID","eigenexonID","constitutive")]

## Merge exons with respective eigenexon ID
eigenexonsE<-merge(eigenexons_assignments, nasoniaE, by.x="exonID", by.y="row.names", all.x=F, all.y=F)
eigenexonsE_con<-eigenexonsE[which(eigenexonsE$constitutive=="constitutive"),]
eigenexonsE_fac<-eigenexonsE[which(eigenexonsE$constitutive=="facultative"),]

## Average expression in adults and gonads by eigenexons
eigenexonsE_con<-ddply(eigenexonsE_con, .(eigenexonID), numcolwise(median), .progress = "text")
eigenexonsE_fac<-ddply(eigenexonsE_fac, .(eigenexonID), numcolwise(median), .progress = "text")

## Divide spliced nodes by constitutive ones
eigenexonsE_fac2<-apply(eigenexonsE_fac, 1, function(x){
  as.numeric(x[-1])/eigenexonsE_con[grep(str_extract(x["eigenexonID"], "^[[:alnum:]]*."), eigenexonsE_con$eigenexonID),-1]
})
eigenexonsE_fac2<-ldply(eigenexonsE_fac2)
eigenexonsE_fac2<-cbind(eigenexonsE_fac[,1], eigenexonsE_fac2)
names(eigenexonsE_fac2)[1]<-"eigenexonID"

# set infinity scores to NA, set scores greater than 1 to 1
# reasons: higher than 1 result from noise on small-scale measurements, therefore possibly complete signal plus noise
# infinity results in constitutive exons ranking below cutoff in samples where nc exons are above cutoff, therefore possibly spurious signal

eigenexonsE_fac2[,sapply(eigenexonsE_fac2, is.numeric)]<-apply(eigenexonsE_fac2[,sapply(eigenexonsE_fac2, is.numeric)],c(1,2),function(x){as.numeric(ifelse(x==Inf, 0, x))})
eigenexonsE_fac2[,sapply(eigenexonsE_fac2, is.numeric)]<-apply(eigenexonsE_fac2[,sapply(eigenexonsE_fac2, is.numeric)],c(1,2),function(x){as.numeric(ifelse(x>1, 1, x))})

## stack data-frames
eigenexonsE<-rbind(eigenexonsE_con, eigenexonsE_fac2)
eigenexonsE<-eigenexonsE[order(eigenexonsE$eigenexonID),]

## Save as csv
write.csv(eigenexonsE, file=file.path(newdir, "gonadadult_eigenexonsE.csv"))
## Load data
eigenexonsE <- read.csv(file = "./Output/splicing_gonadvalidation/gonadadult_eigenexonsE.csv")[,-1]

# Remove sex-specific nodes
# eigenexonsE<- eigenexonsE[which(eigenexonsE$eigenexonID%in%sexspec[which(is.na(sexspec$spec==F)),"eigenexonID"]),]

# Convert to expression list
expr_list<-list(E=as.matrix(eigenexonsE[,-1]), genes=as.character(eigenexonsE[,1]), targets=names(eigenexonsE[,-1]))
# Converting to LIMMA EList format
expr_list<-new("EList", expr_list)

# Z transform separately transcription and splicing nodes, setting minimum to zero (rescales to a common distribution for easier comparisons)
z<-function(x){
  #   y<-x-mean(x, na.rm = T)
  y<-x-(min(x, na.rm = T)*sign(min(x, na.rm = T)))
  y<-y/sd(y, na.rm = T)
}

expr_list$E[grep("fac", expr_list$genes),]<-z(expr_list$E[grep("fac", expr_list$genes),])
expr_list$E[grep("con", expr_list$genes),]<-z(expr_list$E[grep("con", expr_list$genes),])

# Create factors
tissue<-as.factor(ifelse(grepl("adult", expr_list$targets), "Adult", "Gonads"))
sex<-as.factor(ifelse(grepl("female|ovaries", expr_list$targets), "Female", "Male"))
# Create design matrix
design<-model.matrix(~0+sex/tissue, contrasts.arg=list(tissue="contr.treatment",sex="contr.treatment"))
as.data.frame(cbind(expr_list$targets, as.character(tissue), as.character(sex), design))

glms<-lmFit(expr_list, design)
glms2<-eBayes(glms) ## detects zero variant genes (not expressed in either adults nor gonads)

# Save exon models and EList
save(file=file.path(newdir, "LIMMA_sexbias"), list=c("glms2"))

# Save abridged list with only expressed genes
glms_gonads <- data.frame(cbind(glms2[["genes"]],glms2[["coefficients"]],glms2[["p.value"]]))
glms_gonads <- na.exclude(glms_gonads)
names(glms_gonads) <- c("eigenexonID","coefFemale","coefMale","coefFemaleGonads","coefMaleGonads","pvalFemale","pvalMale","pvalFemaleGonads","pvalMaleGonads")

# This version applies FDR correction to all p-values
FDRfitall<-fdrtool(as.numeric(as.matrix(glms_gonads[,6:9])), statistic="pvalue", verbose=T) # strips vector columnwise

Fdrfitall<-matrix(FDRfitall$qval, 
                  ncol=4, nrow=nrow(glms_gonads), dimnames=list(glms_gonads$eigenexonID, paste(colnames(glms_gonads)[6:9],"_qval", sep=""))) # add global (tail area based)
fdrfitall<-matrix(FDRfitall$lfdr, 
                  ncol=4, nrow=nrow(glms_gonads), dimnames=list(glms_gonads$eigenexonID, paste(colnames(glms_gonads)[6:9],"_lfdr", sep=""))) # add global (tail area based)

# Save Fdr and fdr to base model
glms_gonads <- cbind(glms_gonads[,1:5],fdrfitall)

# and save as csv
write.csv(glms_gonads, file=file.path(newdir, "LIMMA_gonads_sexbias_fdr.csv"))

# Import compiled dataset
gonadannot <- read.csv(file = "./Output/splicing_gonadvalidation/LIMMA_gonads_sexbias_fdr.csv")[,-1]

# Annotate ovary and testes biased genes
testesannot <- ddply(.data = gonadannot, .variables = .(eigenexonID), .progress="text", summarize, 
  testes = as.factor(if (pvalMaleGonads_lfdr<0.05) {
    if (coefMaleGonads>0) {"Testes_biased"} else if (coefMaleGonads<0) {"Maleadult_biased"} else {"Unbiased"}
  } else {NA})
)
ovariesannot <- ddply(.data = gonadannot, .variables = .(eigenexonID), .progress="text", summarize, 
                     ovaries = as.factor(if (pvalFemaleGonads_lfdr<0.05) {
                       if (coefFemaleGonads>0) {"Ovaries_biased"} else if (coefFemaleGonads<0) {"Femaleadult_biased"} else {"Unbiased"}
                     } else {NA})
)
gonadannot <- merge(testesannot, ovariesannot)

summary(gonadannot)
table(gonadannot$ovaries, gonadannot$testes, useNA = "always")
margin.table(table(gonadannot$ovaries, gonadannot$testes, useNA = "always"), margin = 1)
margin.table(table(gonadannot$ovaries, gonadannot$testes, useNA = "always"), margin = 2)

table(gonadannot$ovaries, gonadannot$testes, grepl("con", gonadannot$eigenexonID), useNA = "always")
margin.table(table(gonadannot$ovaries, gonadannot$testes, grepl("con", gonadannot$eigenexonID), useNA = "always"), margin = c(1,3))
margin.table(table(gonadannot$ovaries, gonadannot$testes, grepl("con", gonadannot$eigenexonID), useNA = "always"), margin = c(2,3))

clusterdata <- read.csv(file = "./Output/Results_compiler/clusterdata_full.csv")[,-1]

table(clusterdata[grep("Ovari", clusterdata$gonadbias),"devsexbias"], clusterdata[grep("Ovari", clusterdata$gonadbias),"node_type"])
prop.table(table(grepl("f$", clusterdata[grep("Ovari", clusterdata$gonadbias),"devsexbias"]), clusterdata[grep("Ovari", clusterdata$gonadbias),"node_type"]), margin = 2)

table(clusterdata[grep("Ovaries_biasedNA", clusterdata$gonadbias),"devsexbias"], clusterdata[grep("Ovaries_biasedNA", clusterdata$gonadbias),"node_type"])

# How many genes have at least one isoform with emb10 fembias?
length(unique(clusterdata[grepl("^f", clusterdata$devsexbias),"geneID"]))
table(clusterdata[grepl("^f", clusterdata$devsexbias),"gonadbias"], clusterdata[grepl("^f", clusterdata$devsexbias),"node_type"], useNA = "always")
margin.table(table(clusterdata[grepl("^f", clusterdata$devsexbias),"gonadbias"], clusterdata[grepl("^f", clusterdata$devsexbias),"node_type"], useNA = "always"), margin = 2)
table(is.na(clusterdata[grepl("^f", clusterdata$devsexbias),]$gonadbias), clusterdata[grepl("^f", clusterdata$devsexbias),"node_type"])

# How many of them have at least one gonad or adultbiased isoform as well?
# earlyfemgenegonadbias <- apply(X = table(droplevels(clusterdata[grepl("^f", clusterdata$devsexbias),"geneID"]), droplevels(clusterdata[grepl("^f", clusterdata$devsexbias),"gonadbias"])), MARGIN = c(1,2), FUN = function(x){x>0})
# table(earlyfemgenegonadbias)

clusterdata[(grepl("^f", clusterdata$devsexbias)&grepl("Ovaries|Femaleadult", clusterdata$gonadbias)),]
clusterdata[(grepl("^f", clusterdata$devsexbias)&grepl("Maleadult|Testes", clusterdata$gonadbias)),]

# Plot as venn diagram
gonadannot_venn <- cbind(
  gonadannot$testes=="Maleadult_biased",
  gonadannot$testes=="Testes_biased",
  gonadannot$ovaries=="Femaleadult_biased",
  gonadannot$ovaries=="Ovaries_biased"
)
vennCounts(gonadannot_venn)
vennDiagram(gonadannot_venn, names = c("Adult Male","Testes","Adult Female","Ovaries"))

library(vcd)
mosaicplot(table(gonadannot$testes, gonadannot$ovaries, useNA = "always"), shade = T)
mosaicplot(table(gonadannot$testes, gonadannot$ovaries, useNA = "always"), shade = T, dir = c("h","v"))


# Compact gonadbias 
gonadannot$gonadbias <- as.factor(paste(gonadannot$ovaries, gonadannot$testes, sep=""))
gonadannot$gonadbias[which(gonadannot$gonadbias=="NANA")] <- NA
summary(gonadannot)
write.csv(gonadannot, file=file.path(newdir, "gonadannot.csv"))

##### Now with interaction design (first select gonad biased, then contrast males and females)

design <- model.matrix(~sex*tissue, contrasts.arg=list(tissue="contr.treatment",sex="contr.treatment"))
as.data.frame(cbind(expr_list$targets, as.character(tissue), as.character(sex), design))

glms<-lmFit(expr_list, design)
glms2<-eBayes(glms) ## detects zero variant genes (not expressed in either adults nor gonads)

# Save exon models and EList
save(file=file.path(newdir, "LIMMA_sexbias_interaction"), list=c("glms2"))

# Save abridged list with only expressed genes
glms_gonads <- data.frame(cbind(glms2[["genes"]],glms2[["coefficients"]],glms2[["p.value"]]))
glms_gonads <- na.exclude(glms_gonads)
names(glms_gonads) <- c("eigenexonID","coefIntercept","coefMale","coefGonads","coefMaleGonads","pvalIntercept","pvalMale","pvalGonads","pvalMaleGonads")

# This version applies FDR correction to all p-values
FDRfitall<-fdrtool(as.numeric(as.matrix(glms_gonads[,6:9])), statistic="pvalue", verbose=T) # strips vector columnwise

Fdrfitall<-matrix(FDRfitall$qval, 
                  ncol=4, nrow=nrow(glms_gonads), dimnames=list(glms_gonads$eigenexonID, paste(colnames(glms_gonads)[6:9],"_qval", sep=""))) # add global (tail area based)
fdrfitall<-matrix(FDRfitall$lfdr, 
                  ncol=4, nrow=nrow(glms_gonads), dimnames=list(glms_gonads$eigenexonID, paste(colnames(glms_gonads)[6:9],"_lfdr", sep=""))) # add global (tail area based)

# Save Fdr and fdr to base model
glms_gonads <- cbind(glms_gonads[,1:5],fdrfitall)

# and save as csv
write.csv(glms_gonads, file=file.path(newdir, "LIMMA_gonads_sexbias_interaction_fdr.csv"))

# Import compiled dataset
gonadannot <- read.csv(file = "./Output/splicing_gonadvalidation/LIMMA_gonads_sexbias_interaction_fdr.csv")[,-1]

# Annotate ovary and testes biased genes
sexannot <- ddply(.data = gonadannot, .variables = .(eigenexonID), .progress="text", summarize, 
                     adultsex = as.factor(if (pvalMale_lfdr<0.05) {
                       if (coefMale>0) {"Male_biased"} else if (coefMale<0) {"Female_biased"} else {"Unbiased"}
                     } else {NA})
)
testesannot <- ddply(.data = gonadannot, .variables = .(eigenexonID), .progress="text", summarize, 
                     sex = as.factor(if (pvalMaleGonads_lfdr<0.05) {
                       if (coefMaleGonads>0) {"Testes_biased"} else if (coefMaleGonads<0) {"Ovaries_biased"} else {"Unbiased"}
                     } else {NA})
)
tissueannot <- ddply(.data = gonadannot, .variables = .(eigenexonID), .progress="text", summarize, 
                      tissue = as.factor(if (pvalGonads_lfdr<0.05) {
                        if (coefGonads>0) {"Gonad_biased"} else if (coefGonads<0) {"Adult_biased"} else {"Unbiased"}
                      } else {NA})
)
gonadannot <- merge(testesannot, tissueannot)
gonadannot <- merge(gonadannot, sexannot)
summary(gonadannot)

table(gonadannot$sex, gonadannot$tissue, useNA = "always")

table(gonadannot$sex, grepl("con", gonadannot$eigenexonID), useNA = "always")
prop.table(table(gonadannot$sex, grepl("con", gonadannot$eigenexonID), useNA = "always"), margin = 1)
prop.table(table(gonadannot$sex, grepl("con", gonadannot$eigenexonID), useNA = "always"), margin = 2)

table(gonadannot$adultsex, gonadannot$sex, gonadannot$tissue, useNA = "ifany")
prop.table(table(gonadannot$adultsex, gonadannot$sex, gonadannot$tissue, useNA = "ifany"), margin = c(2,3))

gonadannot_interaction <- gonadannot[,c("eigenexonID", "sex", "tissue")]
names(gonadannot_interaction) <- c("eigenexonID", "gonadbias", "tissuebias")
write.csv(gonadannot_interaction, file=file.path(newdir, "gonadannot_interaction.csv"))

## Import into clusterdata and then import clusterdata

gonadannot <- read.csv(file = "./Output/splicing_gonadvalidation/gonadannot_interaction.csv")[,-1]
clusterdata <- read.csv(file = "./Output/Results_compiler/clusterdata_full.csv")[,-1]

# How many ovary-biased genes with adult malebias are somatic vs gonadal?
prop.table(table(clusterdata[which(clusterdata$node_type=="transcription"),]$gonadbias, grepl(pattern = "m$", clusterdata[which(clusterdata$node_type=="transcription"),]$devsexbias), clusterdata[which(clusterdata$node_type=="transcription"),]$tissuebias, useNA = "ifany"), margin = c(1,3))
# Same with isoforms
prop.table(table(clusterdata[which(clusterdata$node_type=="splicing"),]$gonadbias, grepl(pattern = "m$", clusterdata[which(clusterdata$node_type=="splicing"),]$devsexbias), clusterdata[which(clusterdata$node_type=="splicing"),]$tissuebias, useNA = "ifany"), margin = c(1,3))
# How many testes-biased genes with adult femalebias are somatic vs gonadal?
prop.table(table(clusterdata[which(clusterdata$node_type=="transcription"),]$gonadbias, grepl(pattern = "f$", clusterdata[which(clusterdata$node_type=="transcription"),]$devsexbias), clusterdata[which(clusterdata$node_type=="transcription"),]$tissuebias, useNA = "ifany"), margin = c(1,3))
# Same with isoforms
prop.table(table(clusterdata[which(clusterdata$node_type=="splicing"),]$gonadbias, grepl(pattern = "f$", clusterdata[which(clusterdata$node_type=="splicing"),]$devsexbias), clusterdata[which(clusterdata$node_type=="splicing"),]$tissuebias, useNA = "ifany"), margin = c(1,3))

prop.table(table(clusterdata[which(clusterdata$node_type=="transcription"),]$gonadbias=="Ovaries_biased"&grepl(pattern = "m$", clusterdata[which(clusterdata$node_type=="transcription"),"devsexbias"]), clusterdata[which(clusterdata$node_type=="transcription"),"tissuebias"], dnn = c("ovarian_malebias","tissue"), useNA = "ifany"), margin = 1)

prop.table(table(clusterdata[which(clusterdata$node_type=="transcription"),]$gonadbias=="Testes_biased"&grepl(pattern = "f$", clusterdata[which(clusterdata$node_type=="transcription"),"devsexbias"]), clusterdata[which(clusterdata$node_type=="transcription"),"tissuebias"], dnn = c("testes_femalebias","tissue"), useNA = "ifany"), margin = 1)


prop.table(table(clusterdata[which(clusterdata$gonadbias=="Ovaries_biased"&grepl("m$",clusterdata$devsexbias)),"tissuebias"], clusterdata[which(clusterdata$gonadbias=="Ovaries_biased"&grepl("m$",clusterdata$devsexbias)),"node_type"], useNA = "ifany"), margin = 2)
prop.table(table(clusterdata[which(clusterdata$gonadbias=="Testes_biased"&grepl("f$",clusterdata$devsexbias)),"tissuebias"], clusterdata[which(clusterdata$gonadbias=="Testes_biased"&grepl("f$",clusterdata$devsexbias)),"node_type"], useNA = "ifany"), margin = 2)


prop.table(table(grepl(pattern = "m.$",clusterdata$devsexbias), clusterdata$gonadbias, clusterdata$node_type, useNA = "ifany", dnn = c("pupal_malebias","gonadbias", "node type")), margin = c(2,3))
prop.table(table(grepl(pattern = "m$",clusterdata$devsexbias), clusterdata$gonadbias, clusterdata$node_type, useNA = "ifany", dnn = c("adult_malebias","gonadbias", "node type")), margin = c(2,3))


table(clusterdata$devsexbias, clusterdata$gonadbias, clusterdata$node_type)
table(clusterdata$devsexbias, clusterdata$gonadbias, clusterdata$node_type, useNA = "ifany")

table(ifelse(grepl(pattern = "f$", clusterdata$devsexbias), "fembias", "unbiased"), ifelse(is.na(clusterdata$gonadbias),"gonadbias","unbiased"), clusterdata$node_type, dnn = c("adult_fembias","gonadbias", "node_type"))

prop.table(table(ifelse(grepl(pattern = "f$", clusterdata$devsexbias), "fembias", "unbiased"), ifelse(is.na(clusterdata$gonadbias),"gonadbias","unbiased"), clusterdata$node_type, dnn = c("adult_fembias","gonadbias", "node_type")), margin = c(1,3))


# are genes with maleadult bias and ovary bias also gonad biased?
table(clusterdata$gonadbias, grepl(pattern = "m$", x = clusterdata$devsexbias), clusterdata$tissuebias, useNA = "ifany")

# # How many genes are true testes-specific-bias? Calculate by removing anything with adult female bias 
# table(grepl(pattern = "[^f]$", x = clusterdata$devsexbias)&(clusterdata$gonadbias=="Testes_biased"), clusterdata$node_type, useNA = "ifany")

#### Venn diagram of testes bias, ovary bias, gonad bias and adult bias (putting femalebiased testes genes in gonadbias?)
table(clusterdata$tissuebias, clusterdata$gonadbias, clusterdata$node_type)

pdf(file = file.path(graphdir, "venn_gonad_sexde_beforecorrection.pdf"), width = 8.24, height = 11.69/2)
par(mfrow = c(1,2))
# For Transcription
gonadannot_venn <- cbind(
  clusterdata[which(clusterdata$node_type=="transcription"),]$gonadbias=="Testes_biased",
  clusterdata[which(clusterdata$node_type=="transcription"),]$gonadbias=="Ovaries_biased",
  # clusterdata[which(clusterdata$node_type=="transcription"),]$tissuebias=="Adult_biased",
  clusterdata[which(clusterdata$node_type=="transcription"),]$tissuebias=="Gonad_biased"
)
gonadannot_venn <- apply(X = gonadannot_venn, MARGIN = c(1,2), FUN = function(x){ifelse(is.na(x),F,x)})
vennCounts(gonadannot_venn)
vennDiagram(gonadannot_venn, names = c("Testes","Ovaries","Gonads"), main = "Genes")
# For Splicing
gonadannot_venn <- cbind(
  clusterdata[which(clusterdata$node_type=="splicing"),]$gonadbias=="Testes_biased",
  clusterdata[which(clusterdata$node_type=="splicing"),]$gonadbias=="Ovaries_biased",
  # clusterdata[which(clusterdata$node_type=="splicing"),]$tissuebias=="Adult_biased",
  clusterdata[which(clusterdata$node_type=="splicing"),]$tissuebias=="Gonad_biased"
)
gonadannot_venn <- apply(X = gonadannot_venn, MARGIN = c(1,2), FUN = function(x){ifelse(is.na(x),F,x)})
vennCounts(gonadannot_venn)
vennDiagram(gonadannot_venn, names = c("Testes","Ovaries","Gonads"), main = "Isoforms")

dev.off()

## After correction (testes-biased adult-female genes added to gonads and removed from testes-biased)

pdf(file = file.path(graphdir, "venn_gonad_sexde_aftercorrection.pdf"), width = 8.24, height = 11.69/2)
par(mfrow = c(1,2))
## For transcription
gonadannot_venn <- cbind(
  clusterdata[which(clusterdata$node_type=="transcription"),]$gonadbias=="Testes_biased"&!grepl(pattern = "^[^f|m]*f$", x = clusterdata[which(clusterdata$node_type=="transcription"),]$devsexbias),
  clusterdata[which(clusterdata$node_type=="transcription"),]$gonadbias=="Ovaries_biased",
  clusterdata[which(clusterdata$node_type=="transcription"),]$tissuebias=="Gonad_biased"|(grepl(pattern = "^[^f|m]*f$", x = clusterdata[which(clusterdata$node_type=="transcription"),]$devsexbias)&clusterdata[which(clusterdata$node_type=="transcription"),]$gonadbias!="Ovaries_biased"))
gonadannot_venn <- apply(X = gonadannot_venn, MARGIN = c(1,2), FUN = function(x){ifelse(is.na(x),F,x)})
vennCounts(gonadannot_venn)
vennDiagram(gonadannot_venn, names = c("Testes","Ovaries","Gonads"), main = "Genes")
## For isoforms
gonadannot_venn <- cbind(
  clusterdata[which(clusterdata$node_type=="splicing"),]$gonadbias=="Testes_biased"&!grepl(pattern = "^[^f|m]*f$", x = clusterdata[which(clusterdata$node_type=="splicing"),]$devsexbias),  clusterdata[which(clusterdata$node_type=="splicing"),]$gonadbias=="Ovaries_biased",
  clusterdata[which(clusterdata$node_type=="splicing"),]$tissuebias=="Gonad_biased"|(grepl(pattern = "^[^f|m]*f$", x = clusterdata[which(clusterdata$node_type=="splicing"),]$devsexbias)&clusterdata[which(clusterdata$node_type=="splicing"),]$gonadbias!="Ovaries_biased"))
gonadannot_venn <- apply(X = gonadannot_venn, MARGIN = c(1,2), FUN = function(x){ifelse(is.na(x),F,x)})
vennCounts(gonadannot_venn)
vennDiagram(gonadannot_venn, names = c("Testes","Ovaries","Gonads"), main = "Isoforms")

dev.off()

######################

## Reshape as ID+stage+sex*expr
eigenexonsE_2<-melt(eigenexonsE)
eigenexonsE_2$stage<-ifelse(grepl("male",eigenexonsE_2$variable), "adult", "gonads")
eigenexonsE_2$sex<-ifelse(grepl("female|ovaries",eigenexonsE_2$variable), "female", "male")
eigenexonsE_2$con<-ifelse(grepl("con", eigenexonsE_2$eigenexonID), "constitutive", "facultative")

## Average the triplicates
eigenexonsE_2<-ddply(eigenexonsE_2, .variables = .(eigenexonID, stage, sex, con), summarize, E=median(value), .progress="text")
eigenexonsE_2$stage<-as.factor(eigenexonsE_2$stage)
eigenexonsE_2$sex<-as.factor(eigenexonsE_2$sex)
eigenexonsE_2$con<-as.factor(eigenexonsE_2$con)

## Reshape to enable direct adult/gonad comparisons
eigenexonsE_ratios<-recast(data = eigenexonsE_2, formula = eigenexonID+sex+con~stage)
eigenexonsE_ratios$meanE<-(eigenexonsE_ratios$adult+eigenexonsE_ratios$gonads)/2

## Plot Adult/Gonad expression, split by sex and splicing
ggplot(eigenexonsE_ratios, aes(x=adult, y=gonads))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+geom_point(alpha=0.1)+theme_bw()+geom_smooth(method="lm")

## add column for male and female expression
eigenexonsE_ratios_2 <- recast(data = eigenexonsE_2, formula = eigenexonID+con~stage+sex)
eigenexonsE_ratios_2 <- merge(eigenexonsE_ratios_2, gonadannot, by = "eigenexonID", all.x = T)
## Generate plot of male/testes vs female/gonads
ggplot(data = eigenexonsE_ratios_2, aes(x=adult_female/gonads_female, y=adult_male/gonads_male, col = tissuebias, alpha = 0.2))+geom_point()+scale_x_log10()+scale_y_log10()+geom_smooth(method="lm")+geom_abline(a=0, b=1)+geom_hline(yintercept=1)+geom_vline(xintercept=1)+facet_grid(.~con)
ggplot(data = eigenexonsE_ratios_2, aes(x=adult_female/gonads_female, y=adult_male/gonads_male, col = gonadbias, alpha = 0.2))+geom_point()+scale_x_log10()+scale_y_log10()+geom_smooth(method="lm")+geom_abline(a=0, b=1)+geom_hline(yintercept=1)+geom_vline(xintercept=1)+facet_grid(.~con)

pdf(file = file.path(graphdir, "maletestes_femalegonads_expr.pdf"), width = 20, height = 20)
ggplot(data = eigenexonsE_ratios_2, aes(x=adult_female/gonads_female, y=adult_male/gonads_male, col = gonadbias, alpha = 1-is.na(gonadbias)))+geom_point()+scale_x_log10()+scale_y_log10()+geom_smooth(method="lm")+geom_abline(a=0, b=1)+geom_hline(yintercept=1)+geom_vline(xintercept=1)+facet_grid(tissuebias~con)+theme_bw()
dev.off()

pdf(file = file.path(graphdir, "maletestes_femalegonads_expr.pdf"), width = 20, height = 20)
ggplot(data = eigenexonsE_ratios_2, aes(x=adult_female/gonads_female, y=adult_male/gonads_male, alpha=0.1))+geom_point()+scale_x_log10()+scale_y_log10()+geom_abline(a=0, b=1)+geom_hline(yintercept=1)+geom_vline(xintercept=1)+theme_bw()
dev.off()

## Compare distributions of gonadal bias
ggplot(eigenexonsE_ratios, aes(x=adult/gonads))+facet_grid(sex~con, scales = "fixed")+geom_density()+theme_bw()+scale_x_log10()

## Compare gonadal bias to average expression
ggplot(na.exclude(eigenexonsE_ratios), aes(x=meanE, y=gonads/adult))+facet_grid(sex~con, scales = "free_x")+geom_density2d()+geom_point(alpha=0.1)+theme_bw()+scale_y_log10()

## Calculate correlations between adults and gonads in m/f con/spl
adult_gonad_correlations<-ddply(eigenexonsE_ratios, .variables = .(sex, con), summarize, correlation_pearson=cor(adult, gonads, use = "pairwise",  method="pearson"), correlation_spearman=cor(adult, gonads, use = "pairwise",  method="spearman"), correlation_kendall=cor(adult, gonads, use = "pairwise",  method="kendall"), .progress = "text")
adult_gonad_correlations

write.csv(adult_gonad_correlations, file=file.path(newdir, "adult_gonad_correlations.csv"))

## Load developmental self-information scores
self.information<- read.csv("./Output/splicing_entropy/eigenexon_self-information.csv")

## Attach self-information to dataset
eigenexonsE_ratios_2<-merge(eigenexonsE_ratios, self.information[,-1], by="eigenexonID", all.x=F, all.y=F)
eigenexonsE_ratios_2$inf_coef_discr<-cut(x = eigenexonsE_ratios_2$inf_coef, breaks = 10)
eigenexonsE_ratios_2$group<-paste(eigenexonsE_ratios_2$sex, eigenexonsE_ratios_2$con, eigenexonsE_ratios_2$inf_coef)

## Plot correlations vs self-information (hypothesis, more specific nodes will correlate better between adults and gonads if the adult signal is gonadal, otherwise vice versa)
ggplot(na.exclude(eigenexonsE_ratios_2), aes(x=adult, y=gonads, alpha=inf_coef))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+geom_point()+theme_bw()

## Same, excluding nodes that are zero in either sample
eigenexonsE_ratios_2$adult1<-ifelse(eigenexonsE_ratios_2$adult>0,1,0)
eigenexonsE_ratios_2$gonad1<-ifelse(eigenexonsE_ratios_2$gonads>0,1,0)

inf_coef_cont1<-ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(x=adult, y=gonads, alpha=inf_coef))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+geom_point()+theme_bw()

ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(x=adult, y=gonads, alpha=0.5, col=inf_coef>0.45))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+geom_point()+theme_bw()

inf_coef_cont2<-ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(x=adult, y=gonads, alpha=inf_coef, col=inf_coef>0.45))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+theme_bw()

pdf(file = file.path(graphdir, "inf_coef_contours2.pdf"), height = 20, width = 20)
multiplot(inf_coef_cont1, inf_coef_cont2)
dev.off()

pdf(file = file.path(graphdir, "inf_coef_scatterplot.pdf"), height = 10, width = 10)
inf_coef_cont1
dev.off()

pdf(file = file.path(graphdir, "inf_coef_contours.pdf"), height = 10, width = 10)
inf_coef_cont2
dev.off()

# Clip to remove nonexpressed in either stage

genbias<-ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2&eigenexonsE_ratios_2$con=="constitutive"),]), aes(x=adult, y=gonads, alpha=inf_coef, col=inf_coef>0.55))+facet_wrap(~sex+con, scales = "free", nrow = 2)+geom_density2d()+theme_bw()+xlim(c(2,6.5))+ylim(c(2,6.5))+geom_smooth(method="lm")+scale_color_brewer(type = "qual", name="Specificity", labels=c("low","high"))+theme(legend.position="bottom")
splbias<-ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2&eigenexonsE_ratios_2$con=="facultative"),]), aes(x=adult, y=gonads, alpha=inf_coef, col=inf_coef>0.55))+facet_wrap(~sex+con, scales = "free", nrow = 2)+geom_density2d()+theme_bw()+xlim(c(0.2,1))+ylim(c(0.2,1))+geom_smooth(method="lm")+scale_color_brewer(type = "qual", name="Specificity", labels=c("low","high"))+theme(legend.position="bottom")

pdf(file=file.path(graphdir, "inf_coef_clipped_contours.pdf"), height = 10, width = 10)
multiplot(genbias, splbias, cols = 2)
dev.off()

# Plot gonadbias vs information coef
ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(y=abs(log10(gonads/adult)), x=inf_coef, col=sex, lty=con))+theme_bw()+scale_x_log10()+geom_boxplot(aes(group=group))+geom_smooth(method="lm")

# Split between gonad and adult biased
eigenexonsE_ratios_2$gonadbias<-(eigenexonsE_ratios_2$gonads/eigenexonsE_ratios_2$adult)<0
ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(y=log10(gonads/adult), x=inf_coef, col=sex, lty=con))+theme_bw()+scale_x_log10()+geom_boxplot(aes(group=group))+geom_smooth(method="lm")+facet_grid(.~gonadbias)

# Tabulate gonad/adult specific nodes
spectable<-table(eigenexonsE_ratios_2$adult1, eigenexonsE_ratios_2$gonad1, eigenexonsE_ratios_2$sex, eigenexonsE_ratios_2$con, dnn=c("adult", "gonad", "sex", "con"))
prop.table(table(eigenexonsE_ratios_2$adult1, eigenexonsE_ratios_2$gonad1, eigenexonsE_ratios_2$sex, eigenexonsE_ratios_2$con, dnn=c("adult", "gonad", "sex", "con")), margin = c(3,4))

# Plot tabulations
pdf(file = file.path(graphdir, "cotabplot_gonadspecific.pdf"))
cotabplot(x = spectable, formula = ~adult+con+gonad+sex, shade = T, panel = cotab_coindep, n=5000)
dev.off()

# Plot differences in self-information across different tabulations
eigenexonsE_ratios_2$express1<-as.factor(paste(eigenexonsE_ratios_2$adult1, eigenexonsE_ratios_2$gonad1))

ggplot(na.exclude(eigenexonsE_ratios_2), aes(x=express1, y=inf_coef))+geom_boxplot(varwidth = T, notch = T)+facet_grid(sex~con)+theme_bw()+scale_x_discrete(labels=c("Gonad \n specific", "Adult \n specific", "Aspecific"), name="")+scale_y_continuous(name="self_information")

eigenexonsE_ratios_2$group1<-paste(eigenexonsE_ratios_2$sex, eigenexonsE_ratios_2$con, eigenexonsE_ratios_2$express1)
ggplot(na.exclude(eigenexonsE_ratios_2), aes(x=group1, y=inf_coef, group=group1, col=sex, lty=con))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+scale_x_discrete(labels=c("Gonad \n specific", "Adult \n specific", "Aspecific"))

pdf(file=file.path(graphdir, "self_inf_gonadspecific.pdf"), width = 15, height = 10)
ggplot(na.exclude(eigenexonsE_ratios_2), aes(x=express1, y=inf_coef))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+scale_x_discrete(labels=c("Gonad \n specific", "Adult \n specific", "Aspecific"), name="")+scale_y_continuous(name="self_information")+facet_wrap(~con+sex, nrow=1)
dev.off()

cotabplot_data<-droplevels(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$express1%in%c("1 1", "0 1", "1 0")),])
cotabplot_data$express1<-mapvalues(cotabplot_data$express1, from = c("1 1", "0 1", "1 0"), to = c("Aspecific","Gonad\nSpecific", "Adult\nSpecific"))
cotabplot_data$inf_coef<-mapvalues(as.factor(cotabplot_data$inf_coef), from = levels(as.factor(cotabplot_data$inf_coef)), to = c(1:10))

pdf(file=file.path(graphdir, "cotabplot_express_selfinf.pdf"), width = 12, height = 12)
cotabplot(~inf_coef+express1|sex+con, data=na.exclude(cotabplot_data), split_vertical=T, shade = T, gp= shading_max, 
          labeling_args = list( rot_labels = c(top=0, right=0),
                                offset_labels = c(right = 0.7),
                                tl_labels=c(T,F),
                                varnames = c(T,F),
                                set_varnames = list(inf_coef="Self-Information Rank"))
)
dev.off()

pdf(file=file.path(graphdir, "cotabplot_express_selfinf_expected.pdf"), width = 12, height = 12)
cotabplot(~inf_coef+express1|sex+con, type="expected", data=na.exclude(cotabplot_data), split_vertical=T, shade = T, gp= shading_max, 
          labeling_args = list( rot_labels = c(top=0, right=0),
                                offset_labels = c(right = 0.7),
                                tl_labels=c(T,F),
                                varnames = c(T,F),
                                set_varnames = list(inf_coef="Self-Information Rank"))
)
dev.off()

pdf(file=file.path(graphdir, "cotabplot_gonadspec_testes.pdf"), width = 12, height = 12)
cotabplot(~con+sex, data=na.exclude(cotabplot_data[which(cotabplot_data$express1=="Gonad\nSpecific"),]), split_vertical=T, shade = T, gp= shading_max,
          labeling_args = list( rot_labels = c(top=0, right=0)),
                                offset_labels = c(right = 0.7),
                                tl_labels=c(T,T),
                                varnames = c(F,F)
                                )
dev.off()

## Plot density distributions
ggplot(na.exclude(eigenexonsE_ratios_2), aes(x=adult/gonads, col=inf_coef_discr))+facet_wrap(~sex+con, scales = "fixed")+theme_bw()+geom_density()+scale_x_log10()

## Save workspace
save.image(file = file.path(newdir, "splicing_gonadvalidation.R"))


# # Analyze tabulations as lm
# library(lme4)
# eigenexonsE_ratios_2$randomint<-as.factor(eigenexonsE_ratios_2$inf_coef>0.75)
# lm_infcoef<-glmer(droplevels(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$express1%in%c("0 1","1 0", "1 1")),]), 
#                formula = inf_coef>0.75~-1+sex*con*express1+(1|randomint), family=binomial)
# plot(lm_infcoef)

# ### Perform separate GLM testing for male adults vs testes and female adults vs ovaries
# 
# femaleigenexonsE <- eigenexonsE[,c(1,grep("female|ovaries", names(eigenexonsE)))]
# 
# # Convert to expression list
# f_expr_list<-list(E=as.matrix(femaleigenexonsE[,-1]), genes=as.character(femaleigenexonsE[,1]), targets=names(femaleigenexonsE[,-1]))
# # Converting to LIMMA EList format
# f_expr_list<-new("EList", f_expr_list)
# 
# # Normalize genes and isoforms separately
# f_expr_list$E[grep("fac", f_expr_list$genes),]<-z(f_expr_list$E[grep("fac", f_expr_list$genes),])
# f_expr_list$E[grep("con", f_expr_list$genes),]<-z(f_expr_list$E[grep("con", f_expr_list$genes),])
# 
# # Create factors
# tissue<-as.factor(ifelse(grepl("adult", f_expr_list$targets), "Adult", "Gonads"))
# # Create design matrix
# design<-model.matrix(~tissue, contrasts.arg=list(tissue="contr.treatment"))
# as.data.frame(cbind(f_expr_list$targets, as.character(tissue), design))
# 
# glms<-lmFit(f_expr_list, design)
# glms2<-eBayes(glms) ## detects zero variant genes (not expressed in either adults nor gonads)
# 
# # Save exon models and EList
# save(file=file.path(newdir, "LIMMA_sexbias_femonly.csv"), list=c("glms2"))
# 
# # Save abridged list with only expressed genes
# glms_ovaries <- data.frame(cbind(glms2[["genes"]],glms2[["coefficients"]],glms2[["p.value"]]))
# glms_ovaries <- na.exclude(glms_ovaries)
# names(glms_ovaries) <- c("eigenexonID","Intercept_coef","Ovaries_coef","Intercept_pval","Ovaries_pval")
# 
# # This version applies FDR correction to all p-values
# FDRfitall<-fdrtool(as.numeric(as.matrix(glms_ovaries[,4:5])), statistic="pvalue", verbose=T) # strips vector columnwise
# 
# Fdrfitall<-matrix(FDRfitall$qval, 
#                   ncol=2, nrow=nrow(glms_ovaries), dimnames=list(glms_ovaries$eigenexonID, paste(colnames(glms_ovaries)[4:5],"_qval", sep=""))) # add global (tail area based)
# fdrfitall<-matrix(FDRfitall$lfdr, 
#                   ncol=2, nrow=nrow(glms_ovaries), dimnames=list(glms_ovaries$eigenexonID, paste(colnames(glms_ovaries)[4:5],"_lfdr", sep=""))) # add global (tail area based)
# 
# # Save Fdr and fdr to base model
# glms_ovaries <- cbind(glms_ovaries[,1:3],fdrfitall)
# 
# # and save as csv
# write.csv(glms_ovaries, file=file.path(newdir, "LIMMA_ovaries_sexbias_fdr.csv"))

# ############# begin untested code 
# # Import compiled dataset
# ovariesannot <- read.csv(file = "./Output/splicing_gonadvalidation/LIMMA_ovaries_sexbias_fdr.csv")[,-1]
# 
# # Annotate ovary biased genes
# ovariesannot <- ddply(.data = ovariesannot, .variables = .(eigenexonID), .progress="text", summarize, 
#                       ovaries = as.factor(if (Ovaries_pval_lfdr<0.05) {
#                         if (Ovaries_coef>0) {"Ovaries_biased"} else if (Ovaries_coef<0) {"Femaleadult_biased"} else {"Unbiased"}
#                       } else {NA})
# )
# 
# #################### end untested code
