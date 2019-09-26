## Check whether methylation in adult females can be predicted from 
# Female bias in adults
# Female bias in any stage (nested in fembias)
# Expression level
# Female-specificity

# Initialize script
date()
# rm(list=setdiff(ls(),"ODB8_EukOGs_genes"))
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(limma)
library(fdrtool)
library(MuMIn)
rm(list = ls())
options(na.action = na.fail)
z<-function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
}
# initialize output path
newdir<-file.path(getwd(), "Output/methylation_devbias")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/methylation_devbias")
dir.create(graphdir)

# Load data
NVIT_OGS2_goodannotcomplete <- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")[,-1]
clusterdata <- read.csv("./Output/Results_compiler/clusterdata_full.csv")[,-1]

femmeth <- clusterdata[,c("eigenexonID","geneID","devsexbias","spec", "adult_female_meth_status", "node_type", "splicing_cathegory")]

expr_str <- NVIT_OGS2_goodannotcomplete[,c("geneID","quality2")]
femmeth <- merge(femmeth, expr_str, by = "geneID")

entropy <- read.csv(file = "./Output/splicing_entropy/eigenexon_self-information.csv")[,-1]
femmeth <- merge(femmeth, entropy, all.x = T, all.y = F)

dim(femmeth)
summary(femmeth)

# femmeth$meth<- as.factor(ifelse(femmeth$adult_female_meth_status=="Methylated", "Methylated", "Other"))
# femmeth$meth<- femmeth$adult_female_meth_status=="Methylated"
femmeth <- femmeth[-which(is.na(femmeth$adult_female_meth_status)),]
femmeth$meth <- femmeth$adult_female_meth_status
femmeth$meth <- relevel(femmeth$meth, ref = "Unmethylated")

femmeth$spec <- as.character(femmeth$spec)
femmeth$spec[which(is.na(femmeth$spec))] <- "Aspecific"
femmeth$spec <- as.factor(femmeth$spec)

# femmeth$adultfem <-grepl("f$", femmeth$devsexbias)
# femmeth$adultmal <-grepl("m$", femmeth$devsexbias)
# 
# femmeth$fem <- grepl("f", femmeth$devsexbias)
# femmeth$mal <- grepl("m", femmeth$devsexbias)

femmeth$bias <- as.factor(sapply(X = femmeth$devsexbias, FUN = function(x)
  if (grepl("f", x)) {"Fembias"} else if (grepl("m", x)) {"Malebias"} else {"Unbiased"}))
femmeth$bias <- relevel(femmeth$bias, ref = "Unbiased")

femmeth$adult <- as.factor(sapply(X = femmeth$devsexbias, FUN = function(x)
  if (grepl("f$", x)) {"Fembias"} else if (grepl("m$", x)) {"Malebias"} else {"Unbiased"}))
femmeth$adult <- relevel(femmeth$adult, ref = "Unbiased")

femmeth$pupa <- as.factor(sapply(X = femmeth$devsexbias, FUN = function(x)
  if (grepl("f.$", x)) {"Fembias"} else if (grepl("m.$", x)) {"Malebias"} else {"Unbiased"}))
femmeth$pupa <- relevel(femmeth$pupa, ref = "Unbiased")

femmeth$quality2 <- as.factor(ifelse(femmeth$quality2=="Express:Strong", "StrongExprSupport", "WeakExprSupport"))
femmeth$quality2 <- relevel(femmeth$quality2, ref = "WeakExprSupport")
femmeth <- femmeth[,-(which(names(femmeth)%in%c("devsexbias","adult_female_meth_status")))]
femmeth <- na.exclude(femmeth)

femmeth$inf_coef <- as.numeric(as.factor(rank(femmeth$inf_coef)))

row.names(femmeth) <- femmeth$eigenexonID

# Use glm
femmeth_glm<-glm(data = femmeth, formula = meth~spec+bias*inf_coef+adult+quality2+node_type+splicing_cathegory, family=binomial)
summary(femmeth_glm)
plot(femmeth_glm)

femmeth_d <- dredge(femmeth_glm)
model.sel(femmeth_d)
summary(model.avg(femmeth_d))

# Initial observations:
# Developmental specificity decreases probability of methylation in adult fem, expr support increases meth prob
# Sex-specificity increases meth prob, more in females (must re-test with nested model)
# Adult female biased expression increases meth prob, adult male biased expr decreases prob
# overall female bias DECREASES support, overall male bias has no support
# Transcription nodes are less likely to be methylated (must check for gene effect, splicing genes are more likely to be meth?)
# Total variance explained ~21%

#### Restrict to transcription nodes alone
femmeth_trans <- droplevels(femmeth[which(femmeth$node_type=="transcription"),])

femmeth_trans_glm <- glm(data = femmeth_trans, formula = meth~spec+bias*inf_coef+adult+quality2+splicing_cathegory, family=binomial)
summary(femmeth_trans_glm)
plot(femmeth_trans_glm)

femmeth_trans_d <- dredge(femmeth_trans_glm)
model.sel(femmeth_trans_d)
summary(model.avg(femmeth_trans_d))

prop.table(table(femmeth_trans$meth, femmeth_trans$spec), margin = 2)
# Methylation among spec genes is ~10 fold lower than genomewide

## Reduce to transcription, high expr support, no sex-specific

femmeth_slim <- droplevels(femmeth[which((femmeth$node_type=="transcription")&(femmeth$quality2=="StrongExprSupport")&(femmeth$spec=="Aspecific")),])
femmeth_slim <- femmeth_slim[,c("geneID","splicing_cathegory","inf_coef","bias","adult","pupa","meth")]
row.names(femmeth_slim) <- femmeth_slim$geneID
summary(femmeth_slim)

femmeth_slim_glm <- glm(data = femmeth_slim, formula = meth~bias*inf_coef+adult+I(pupa=="Malebias")+splicing_cathegory, family=binomial)
# femmeth_slim_glm <- glm(data = femmeth_slim, formula = meth~bias*inf_coef+adult+pupa+splicing_cathegory, family=binomial)
summary(femmeth_slim_glm)
plot(femmeth_slim_glm)

femmeth_slim_d <- dredge(femmeth_slim_glm)
model.sel(femmeth_slim_d)
summary(model.avg(femmeth_slim_d))

femmeth_slim_bestmodel <- glm(data = femmeth_slim, formula = meth~bias*inf_coef+adult+splicing_cathegory, family=binomial)
summary(femmeth_slim_bestmodel)
plot(femmeth_slim_bestmodel)

# ~30% variance explained with top model, best of set, DeltaAICc 1.93, second best includes male pupal expr (positive pred)

#### Graph trials
femmeth_slim_graph <- femmeth_slim
femmeth_slim_graph$meth <- as.numeric(femmeth_slim_graph$meth)
femmeth_slim_graph$grouper <- as.factor(paste(femmeth_slim_graph$adult, femmeth_slim_graph$bias))

# table(femmeth_slim_graph$meth)
femmeth_slim_graph <- ddply(femmeth_slim_graph, .(bias, adult ,meth), summarize, ngenes = length(unique(geneID)))

ggplot(data = femmeth_slim_graph[which(femmeth_slim_graph$meth==1),], aes(y=ngenes, x=bias, fill=adult))+geom_bar(stat = "identity", position="stack")+theme_bw()+scale_y_continuous("Number of Genes")

ggplot(data = femmeth_slim_graph, aes(y=ngenes, x=bias, fill=adult))+geom_bar(stat = "identity", position="stack")+theme_bw()+scale_y_continuous("Number of Methylated Genes")+facet_grid(.~meth)
ggplot(data = femmeth_slim_graph[which(femmeth_slim_graph$bias!="Unbiased"),], aes(y=ngenes, x=bias, fill=adult))+geom_bar(stat = "identity", position="dodge")+theme_bw()+scale_y_continuous("Number of Methylated Genes")+facet_grid(.~meth)
