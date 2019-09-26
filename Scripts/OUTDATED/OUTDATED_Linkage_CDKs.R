###  Outdated linkage and other analyses

###### outdated linkage analyses
# load annotation
OGS2.0<- read.csv("../OGS2/GeneStats/NVIT_OGS2_goodannotcomplete.csv")

# how many annotated genes are represented on the array?
table(OGS2.0$geneID%in%fit2$genes)
prop.table(table(OGS2.0$geneID%in%fit2$genes)) #95%
# how many genes on the array are absent from the annotation?
fit2$genes[which(fit2$genes%in%OGS2.0$geneID==F)] # only one, excluded because of poor support

# Which CDKs are present in the OGS2?
OGS2.0[grep("yclin",OGS2.0$Name),"Name"]
# what are their gene names?
OGS2.0[grep("yclin",OGS2.0$Name),"geneID"]
# are they expressed in stage 0 vs the rest of development?
Fdrfitall[which(row.names(Fdrfitall)%in%OGS2.0[grep("yclin",OGS2.0$Name),"geneID"]),]
CDKexpr<-fit2$coefficients[which(row.names(fit2$coefficients)%in%OGS2.0[grep("yclin",OGS2.0$Name),"geneID"]),]

library(ggplot2)
library(reshape)
ggplot(melt(CDKexpr[,1:5]), aes(x=X2, y=value, colour=X1, group=X1))+geom_line()


# is there any linkage group enriched for male specific genes?
library(reshape)
# create dataset of DE genes 
sexDE<-(apply(Fdrfitall, c(1,2), function(x){x<0.0001}))
sexDE<-melt(sexDE)
names(sexDE)<-c("geneID","Factor","DE")
row.names(sexDE)<-paste(sexDE$geneID, sexDE$Factor, sep="_")
# add coefficients of change
sexcoef<-fit2$coefficients
sexcoef<-melt(sexcoef)
names(sexcoef)<-c("geneID","Factor","Coefficient")
row.names(sexcoef)<-paste(sexcoef$geneID, sexcoef$Factor, sep="_")

sexDE2<-merge(sexDE, sexcoef, by=0)
sexDE<-sexDE2[,c(2,3,4,7)]
names(sexDE)<-c("geneID","Factor","DE","Coefficient")
# add direction of change
sexDE$DE<-apply(sexDE,1,function(x){
  if(x["DE"]!=FALSE) {
    ifelse(as.numeric(x["Coefficient"])<0, yes="DN", no="UP")}
  else FALSE
})


# add linkage group
sexDE<-merge(sexDE,OGS2.0[,c("geneID","LinkageCluster")], by="geneID")
sexDE$LinkageCluster<-as.factor(sexDE$LinkageCluster)

# how many linkage clusters have DE genes on them?
prop.table(table(sexDEclusters$DEs==0)) # 20%

# plot raw data
library(plyr)
sexDEclusters<-ddply(sexDE[sexDE$Factor=="maleTRUE",], .(LinkageCluster), .fun=summarize, DEs=sum(DE!="FALSE"), Size=length(DE), DEp=sum(DE!="FALSE")/length(DE), Coefficient=median(Coefficient))
sexDEclusters$Chromosome<-as.character(lapply(strsplit(as.character(sexDEclusters$LinkageCluster),".", fixed=T),function(x){x[[1]]}))
ggplot(sexDEclusters, aes(x=LinkageCluster, y=DEp, col=Chromosome))+geom_point()+theme_bw()+scale_y_log10()

# add chromosome and total number of genes on lg to dataset
sexDE<-merge(sexDE, sexDEclusters[,c("LinkageCluster", "Size", "Chromosome")], all.x=T)


# specify model
library(lme4)
# linkage<-glm(DE!="FALSE"~Size+Chromosome/LinkageCluster-1, data=sexDE[sexDE$Factor=="maleTRUE",],family="quasibinomial")
linkage<-glm(DE!="FALSE"~Size+LinkageCluster-1, data=sexDE[sexDE$Factor=="maleTRUE",],family="quasibinomial")
# extract p values
linkageP<-summary(linkage)$coef[,"Pr(>|t|)"]
# apply FDR correction
library(fdrtool)
linkageq<-fdrtool(linkageP)$qval
# extract enrichment coefficients for each linkage group
linkageCoef<-summary(linkage)$coef[,"Estimate"]
# add coefficients, q-values and re-plot
SDC2<-cbind(linkageq, linkageCoef)
SDC2<-SDC2[-1,]
sexDEclusters2<-cbind(sexDEclusters[-(nrow(sexDEclusters)),], SDC2)
rm(SDC2)

ggplot(sexDEclusters2, aes(x=LinkageCluster, y=linkageq*sign(linkageCoef), col=as.factor(sign(Coefficient))))+geom_point()+theme_bw()

# select significant linkage clusters
sexDEclusters2[linkageq<0.001,]
# all of the clusters wich differ in mean count of sex DE genes have lower number of sex DE genes
# which is the intercept for all clusters?
summary(linkage)

# try with poisson distribution
glm(LinkageCluster~DE, data=droplevels(sexDE[sexDE$Factor=="maleTRUE",]),family="poisson")
# not working for strange reasons

# does any linkage group have an higher proportion of male DE genes?
linkageMF<-glm(DE=="UP"~Size+LinkageCluster-1, data=droplevels(sexDE[sexDE$Factor=="maleTRUE"&sexDE$DE!="FALSE",]),family="binomial")
step(linkageMF) # size (number of genes) is not significant in predicting the direction of sex bias in linkage clusters
linkageMF<-glm(DE=="UP"~LinkageCluster-1, data=droplevels(sexDE[sexDE$Factor=="maleTRUE"&sexDE$DE!="FALSE",]),family="binomial")

# extract p values
linkageMFP<-summary(linkageMF)$coef[,"Pr(>|z|)"]
# apply FDR correction
library(fdrtool)
linkageMFq<-fdrtool(linkageMFP)$qval
# extract enrichment coefficients for each linkage group
linkageMFCoef<-summary(linkageMF)$coef[,"Estimate"]
SDC3<-as.data.frame(cbind(linkageMFq, linkageMFCoef))
names(SDC3)<-c("linkageMFq", "linkageMFCoef")
SDC3$LinkageCluster<-substring(row.names(SDC3), 15,nchar(row.names(SDC3)))
# add coefficients, q-values and re-plot
sexDEclusters3<-merge(sexDEclusters2, SDC3, by="LinkageCluster", all.x=T)
rm(SDC3)

# and plot
ggplot(sexDEclusters3, aes(x=LinkageCluster, y=linkageMFq*sign(linkageMFCoef), col=as.factor(Chromosome)))+geom_point()+theme_bw()