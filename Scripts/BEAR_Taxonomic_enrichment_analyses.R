### Test different hypotheses on the evolution of sex-biased clusters
rm(list=ls())
library(lme4)
library(MuMIn)
library(reshape2)
library(stringr)
library(plyr)
options(na.action = "na.fail")
newdir<-file.path(getwd(), "Output/BEAR_taxonomic_enrichment_analyses")

# Load data
eigenexon_data <- read.csv("./Output/Results_compiler/eigenexon_data.csv")

# Reshape data for dWithin sexbias over developmental time of genes in each cluster, split by taxonomic depth and splicing
dWithin_Taxa<-eigenexon_data[,c("geneID","clusterID","eigenexonID","constitutive","Taxascan_max", grep("sexbias.*dWithin",names(eigenexon_data),value = T))]
dWithin_Taxa<-melt(dWithin_Taxa, id.vars = c("geneID","clusterID","eigenexonID","constitutive","Taxascan_max"))
dWithin_Taxa$stage<-str_extract(string = dWithin_Taxa$variable, pattern = "stage[[:alnum:]]*")
dWithin_Taxa$stage<-factor(dWithin_Taxa$stage, levels=c("stageemb10", "stageemb18", "stagelar51", "stagepupyel", "stageadult"))
dWithin_Taxa$value<-as.factor(sign(dWithin_Taxa$value))

# Collapse by cluster
dWithin_Taxa<-ddply(.data = dWithin_Taxa, .variables = .(clusterID, stage), .fun = summarize, value=value[1], Taxascan=mean(Taxascan_max), p.constitutive=sum(constitutive=="constitutive")/length(constitutive), log_clustersize=log10(length(geneID)), log_nGenes=log10(length(unique(geneID))), .progress = "text")

# Check for correlation between cluster size and proportion of constitutive nodes (absent)
ggplot(dWithin_Taxa, aes(x=p.constitutive, y=log_clustersize))+geom_density2d()+geom_point()+geom_smooth()+theme_bw()
# Check for distribution of mean taxonomic depth (controlling for clustersize)
ggplot(dWithin_Taxa, aes(x=Taxascan, y=log_clustersize))+geom_density2d()+theme_bw()+geom_rug()+geom_point()
# Check correlation between number of genes and cluster size (almost perfect log-log correlation)
ggplot(dWithin_Taxa, aes(x=log_nGenes, y=log_clustersize))+geom_density2d()+geom_point()+geom_smooth()+theme_bw()

# Code Z transform
Z<-function(x){(x-mean(x))/sd(x)}
logit<-function(x){
  y<-x+0.001 # adjusted to avoid -Inf on clusters composed only by splicing nodes
  log10(y/(1-y))
}

## Are sexbiased nodes more ancient? Does splicing play a role? Control for structure of cluster/gene and check for stage-specific age classes
GLM_dWithin_Taxa<-glmer(formula = (is.na(value)==F)~poly(as.integer(stage), degree = 4)*Z(Taxascan)*logit(p.constitutive)+(1|log_nGenes), data = dWithin_Taxa, family=binomial)
modelset_dWithin_Taxa<-dredge(GLM_dWithin_Taxa)
model.sel(modelset_dWithin_Taxa)
summary(model.avg(modelset_dWithin_Taxa))

bestGLM_dWithin_Taxa<-glmer(formula = (is.na(value)==F)~poly(as.integer(stage), degree=3)*logit(p.constitutive)+(1|log_nGenes), data = dWithin_Taxa, family=binomial)

# Reshape data for dOut sexbias over developmental time of genes in each cluster, split by taxonomic depth and splicing
dOut_Taxa<-eigenexon_data[,c("geneID","clusterID","eigenexonID","constitutive","Taxascan_max", grep("sexbias.*dOut",names(eigenexon_data),value = T))]
dOut_Taxa<-melt(dOut_Taxa, id.vars = c("geneID","clusterID","eigenexonID","constitutive","Taxascan_max"))
dOut_Taxa$stage<-str_extract(string = dOut_Taxa$variable, pattern = "stage[[:alnum:]]*")
dOut_Taxa$stage<-factor(dOut_Taxa$stage, levels=c("stageemb10", "stageemb18", "stagelar51", "stagepupyel", "stageadult"))
dOut_Taxa$value<-as.factor(sign(dOut_Taxa$value))

## Are sexbiased nodes more ancient? Does splicing play a role? Control for structure of cluster/gene and check for stage-specific age classes
GLM_dOut_Taxa<-glmer(formula = is.na(value)~-1+stage*Taxascan_max*constitutive+(1|clusterID/geneID), data = dOut_Taxa[which(is.na(dOut_Taxa$clusterID)==F),] , family=binomial)
modelset_dOut_Taxa<-dredge(GLM_dOut_Taxa)
model.sel(modelset_dOut_Taxa)
model.avg(modelset_dOut_Taxa)


# save workspace
save.image(file = file.path(newdir, BEAR_Taxonomic_enrichment_analyses.R))
