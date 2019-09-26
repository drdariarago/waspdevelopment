### Validate splicing dynamics via visualization of dsx in male and female development

# Initialize script
rm(list=ls())
library(ggplot2)
library(plyr)
library(reshape2)
library(stringr)

## Import eigenexon Evalues, genome annotation, eigenexon assignments and raw score files
eigenexon_evalues <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexon_evalues.csv")
OGS2 <- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")
eigenexons_assignments <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexons_assignments.csv")
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/2015_05_19_Justin/nasonia_devtesova.exon.r99.score2")

##### Find DSX
grep("Doublesex", OGS2$Name, value = T)
dsx_geneID<-OGS2[grep("Doublesex", OGS2$Name, value = F)[1],"geneID"]

## Extract eigenexons corresponding to Doublesex
grep(dsx_geneID, eigenexon_evalues$eigenexonID, value=T)
## One isoform identified

## Check exon assignment file
grep(dsx_geneID, eigenexons_assignments$eigenexonID, value=T)
eigenexons_assignments[grep(dsx_geneID, eigenexons_assignments$eigenexonID, value=F),]

## Plot dynamics of eigenexons over time
dsx_eigenexons_evalues<-eigenexon_evalues[grep(dsx_geneID, eigenexon_evalues$eigenexonID),]
dsx_eigenexons_evalues<-melt(dsx_eigenexons_evalues)
dsx_eigenexons_evalues$stage<-as.factor(str_extract(string=dsx_eigenexons_evalues$variable, pattern="^[[:alnum:]]*"))
dsx_eigenexons_evalues$stage<-factor(dsx_eigenexons_evalues$stage, levels=c("emb10","emb18","lar51","pupyel","adult", "ovaries", "testes"))
dsx_eigenexons_evalues$sex<-as.factor(ifelse(grepl("female|ovaries", dsx_eigenexons_evalues$variable),"Female","Male"))
dsx_eigenexons_evalues$grouper<-as.factor(paste(dsx_eigenexons_evalues$eigenexonID, dsx_eigenexons_evalues$sex))

ggplot(dsx_eigenexons_evalues, aes(x=stage, y=value, lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~eigenexonID)

## Plot individual exons instead
dsx_exons_evalues<-nasonia_devtesova.exon.r99[grep(dsx_geneID, nasonia_devtesova.exon.r99$EXON, value=F),c(1,grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99)))]
dsx_exons_evalues<-melt(dsx_exons_evalues)

dsx_exons_evalues$stage<-as.factor(str_extract(string=dsx_exons_evalues$variable, pattern="^[[:alnum:]]*"))
dsx_exons_evalues$stage<-factor(dsx_exons_evalues$stage, levels=c("emb10","emb18","lar51","pupyel","adult", "ovaries", "testes"))
dsx_exons_evalues$sex<-as.factor(ifelse(grepl("female|ovaries", dsx_exons_evalues$variable),"Female","Male"))
dsx_exons_evalues$grouper<-as.factor(paste(dsx_exons_evalues$EXON, dsx_exons_evalues$sex))
dsx_exons_evalues<-merge(dsx_exons_evalues, eigenexons_assignments[,c("exonID","eigenexonID")], by.x="EXON", by.y="exonID")
dsx_exons_medians<-ddply(dsx_exons_evalues, .(variable), summarize, medianvalue=median(value))
dsx_exons_evalues<-merge(dsx_exons_evalues, dsx_exons_medians)
# dsx_exons_evalues<-merge(dsx_exons_evalues, dsx_eigenexons_evalues, by.x=c("variable","eigenexonID","stage","sex"), by.y=c("variable","eigenexonID","stage","sex"), suffixes = c("","_eigenexon"))

ggplot(dsx_exons_evalues, aes(x=stage, y=value, col=sex, group=grouper))+geom_point()+geom_smooth()+theme_bw()
ggplot(dsx_exons_evalues, aes(x=stage, y=value, col=eigenexonID, lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)
ggplot(dsx_exons_evalues, aes(x=stage, y=value-medianvalue, col=eigenexonID, lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)

# ggplot(dsx_exons_evalues, aes(x=stage, y=value-value_eigenexon, col=eigenexonID, lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)

## Alternative isoform is female specific in pupa/adult (low expression, possibly noise)

### Check for SXL
grep("Sex-lethal", OGS2$Name, value = T)
sxl_geneID<-OGS2[grep("Sex-lethal", OGS2$Name, value = F)[1],"geneID"]
## Extract eigenexons corresponding to Doublesex
grep(sxl_geneID, eigenexon_evalues$eigenexonID, value=T)

## Check exon assignment file
grep(sxl_geneID, eigenexons_assignments$eigenexonID, value=T)
eigenexons_assignments[grep(sxl_geneID, eigenexons_assignments$eigenexonID, value=F),]

## Plot dynamics of gene over time
sxl_eigenexons_evalues<-eigenexon_evalues[grep(sxl_geneID, eigenexon_evalues$eigenexonID),]
sxl_eigenexons_evalues<-melt(sxl_eigenexons_evalues)
sxl_eigenexons_evalues$stage<-as.factor(str_extract(string=sxl_eigenexons_evalues$variable, pattern="^[[:alnum:]]*"))
sxl_eigenexons_evalues$stage<-factor(sxl_eigenexons_evalues$stage, levels=c("emb10","emb18","lar51","pupyel","adult", "ovaries", "testes"))
sxl_eigenexons_evalues$sex<-as.factor(ifelse(grepl("female|ovaries", sxl_eigenexons_evalues$variable),"Female","Male"))
sxl_eigenexons_evalues$grouper<-as.factor(paste(sxl_eigenexons_evalues$eigenexonID, sxl_eigenexons_evalues$sex))

ggplot(sxl_eigenexons_evalues, aes(x=stage, y=value, lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~eigenexonID)


## Plot individual exons instead
sxl_exons_evalues<-nasonia_devtesova.exon.r99[grep(sxl_geneID, nasonia_devtesova.exon.r99$EXON, value=F),c(1,grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99)))]

sxl_exons_evalues<-melt(sxl_exons_evalues)
sxl_exons_evalues$value<-ifelse(sxl_exons_evalues$value>=0,sxl_exons_evalues$value,0)

sxl_exons_evalues$stage<-as.factor(str_extract(string=sxl_exons_evalues$variable, pattern="^[[:alnum:]]*"))
sxl_exons_evalues$stage<-factor(sxl_exons_evalues$stage, levels=c("emb10","emb18","lar51","pupyel","adult", "ovaries", "testes"))
sxl_exons_evalues$sex<-as.factor(ifelse(grepl("female|ovaries", sxl_exons_evalues$variable),"Female","Male"))
sxl_exons_evalues$grouper<-as.factor(paste(sxl_exons_evalues$EXON, sxl_exons_evalues$sex))
# sxl_exons_evalues$EXON<-paste(str_extract(sxl_exons_evalues$EXON, "[^t]*"),str_extract(sxl_exons_evalues$EXON, "[[:digit:]]*$"), sep = "t0.")
sxl_exons_evalues<-merge(sxl_exons_evalues, eigenexons_assignments[,c("exonID","eigenexonID")], by.x="EXON", by.y="exonID")
sxl_exons_medians<-ddply(sxl_exons_evalues, .(variable), summarize, medianvalue=median(value))
sxl_exons_evalues<-merge(sxl_exons_evalues, sxl_exons_medians)

ggplot(sxl_exons_evalues, aes(x=stage, y=value-medianvalue, col=eigenexonID, lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)

# ggplot(sxl_exons_evalues, aes(x=stage, y=value, col=sex, group=grouper))+geom_point()+geom_smooth(se=F)
ggplot(sxl_exons_evalues, aes(x=stage, y=value, col=eigenexonID, lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)

## Hclust exons
library(amap)
sxl_exons_evalues2<-as.matrix(nasonia_devtesova.exon.r99[grep("Nasvi2EG000104", nasonia_devtesova.exon.r99$EXON, value=F),grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99))])
row.names(sxl_exons_evalues2)<-nasonia_devtesova.exon.r99[grep("Nasvi2EG000104", nasonia_devtesova.exon.r99$EXON, value=F),"EXON"]
sxl_exons_evalues2<-apply(sxl_exons_evalues2, c(1,2), function(x){ifelse(x<=0, 0, x)})
sxl_exons_evalues2<-sxl_exons_evalues2[,-grep("testes|ovaries",colnames(sxl_exons_evalues2))]

sxl_hcluster<-hcluster(sxl_exons_evalues2, method = "correlation", link = "complete")
plot(sxl_hcluster)


# ### Check for TRA
# grep("transformer", OGS2$Name, value = T)
# tra_geneID<-OGS2[grep("transformer", OGS2$Name, value = F)[1],"geneID"]
# ## Extract eigenexons corresponding to Doublesex
# grep(tra_geneID, eigenexon_evalues$eigenexonID, value=T)
# 
# ## Two groups identified (one alternative isoform)
# 
# ## Check exon assignment file
# grep(tra_geneID, eigenexons_assignments$eigenexonID, value=T)
# eigenexons_assignments[grep(tra_geneID, eigenexons_assignments$eigenexonID, value=F),]
# 
# ## Plot dynamics of gene over time
# 
# ## Plot individual exons instead
# tra_exons_evalues<-nasonia_devtesova.exon.r99[grep(tra_geneID, nasonia_devtesova.exon.r99$EXON, value=F),c(1,grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99)))]
# 
# tra_exons_evalues<-melt(tra_exons_evalues)
# tra_exons_evalues$value<-ifelse(tra_exons_evalues$value>=0,tra_exons_evalues$value,0)
# 
# tra_exons_evalues$stage<-as.factor(str_extract(string=tra_exons_evalues$variable, pattern="^[[:alnum:]]*"))
# tra_exons_evalues$stage<-factor(tra_exons_evalues$stage, levels=c("emb10","emb18","lar51","pupyel","adult", "ovaries", "testes"))
# tra_exons_evalues$sex<-as.factor(ifelse(grepl("female|ovaries", tra_exons_evalues$variable),"Female","Male"))
# tra_exons_evalues$grouper<-as.factor(paste(tra_exons_evalues$EXON, tra_exons_evalues$sex))
# tra_exons_evalues<-merge(tra_exons_evalues, eigenexons_assignments[,c("exonID","eigenexonID")], by.x="EXON", by.y="exonID")
# 
# ggplot(tra_exons_evalues, aes(x=stage, y=value, col=sex, group=grouper))+geom_point()+geom_smooth(se=F)
# ggplot(tra_exons_evalues, aes(x=stage, y=value, col=eigenexonID, lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)
# 
# ## Hclust exons
# library(amap)
# tra_exons_evalues2<-as.matrix(nasonia_devtesova.exon.r99[grep(tra_geneID, nasonia_devtesova.exon.r99$EXON, value=F),grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99))])
# row.names(tra_exons_evalues2)<-nasonia_devtesova.exon.r99[grep(tra_geneID, nasonia_devtesova.exon.r99$EXON, value=F),"EXON"]
# tra_exons_evalues2<-apply(tra_exons_evalues2, c(1,2), function(x){ifelse(x<=0, 0, x)})
# tra_exons_evalues2<-tra_exons_evalues2[,-grep("testes|ovaries",colnames(tra_exons_evalues2))]
# 
# tra_hcluster<-hcluster(tra_exons_evalues2, method = "correlation", link = "complete")
# plot(tra_hcluster)
# tra_cutree<-data.frame(cluster=cutree(tra_hcluster, k=4))
# 
# tra_exons_evalues_cutree<-merge(tra_exons_evalues, tra_cutree, by.x="EXON", by.y=0)
# ggplot(tra_exons_evalues_cutree, aes(x=stage, y=value, col=as.factor(cluster), lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)

### Check for TRA2
grep("transformer", OGS2$Name, value = T)
tra2_geneID<-OGS2[grep("transformer", OGS2$Name, value = F)[1],"geneID"]
## Extract eigenexons corresponding to transformer2
grep(tra2_geneID, eigenexon_evalues$eigenexonID, value=T)

## Two groups identified (one alternative isoform)

## Check exon assignment file
grep(tra2_geneID, eigenexons_assignments$eigenexonID, value=T)
eigenexons_assignments[grep(tra2_geneID, eigenexons_assignments$eigenexonID, value=F),]

## Plot dynamics of gene over time

## Plot individual exons instead
tra2_exons_evalues<-nasonia_devtesova.exon.r99[grep(tra2_geneID, nasonia_devtesova.exon.r99$EXON, value=F),c(1,grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99)))]

tra2_exons_evalues<-melt(tra2_exons_evalues)
tra2_exons_evalues$value<-ifelse(tra2_exons_evalues$value>=0,tra2_exons_evalues$value,0)

tra2_exons_evalues$stage<-as.factor(str_extract(string=tra2_exons_evalues$variable, pattern="^[[:alnum:]]*"))
tra2_exons_evalues$stage<-factor(tra2_exons_evalues$stage, levels=c("emb10","emb18","lar51","pupyel","adult", "ovaries", "testes"))
tra2_exons_evalues$sex<-as.factor(ifelse(grepl("female|ovaries", tra2_exons_evalues$variable),"Female","Male"))
tra2_exons_evalues$grouper<-as.factor(paste(tra2_exons_evalues$EXON, tra2_exons_evalues$sex))
# tra2_exons_evalues$EXON<-paste(str_extract(tra2_exons_evalues$EXON, "[^t]*"),str_extract(tra2_exons_evalues$EXON, "[[:digit:]]*$"), sep = "t0.")
tra2_exons_evalues<-merge(tra2_exons_evalues, eigenexons_assignments[,c("exonID","eigenexonID")], by.x="EXON", by.y="exonID")

# ggplot(tra2_exons_evalues, aes(x=stage, y=value, col=sex, group=grouper))+geom_point()+geom_smooth(se=F)
ggplot(tra2_exons_evalues, aes(x=stage, y=value, col=eigenexonID, lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)

## Hclust exons
library(amap)
tra2_exons_evalues2<-as.matrix(nasonia_devtesova.exon.r99[grep(tra2_geneID, nasonia_devtesova.exon.r99$EXON, value=F),grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99))])
row.names(tra2_exons_evalues2)<-nasonia_devtesova.exon.r99[grep(tra2_geneID, nasonia_devtesova.exon.r99$EXON, value=F),"EXON"]
tra2_exons_evalues2<-apply(tra2_exons_evalues2, c(1,2), function(x){ifelse(x<=0, 0, x)})
tra2_exons_evalues2<-tra2_exons_evalues2[,-grep("testes|ovaries",colnames(tra2_exons_evalues2))]

tra2_hcluster<-hcluster(tra2_exons_evalues2, method = "correlation", link = "complete")
plot(tra2_hcluster)
tra2_cutree<-data.frame(cluster=cutree(tra2_hcluster, k=2))

tra2_exons_evalues_cutree<-merge(tra2_exons_evalues, tra2_cutree, by.x="EXON", by.y=0)
ggplot(tra2_exons_evalues_cutree, aes(x=stage, y=value, col=as.factor(cluster), lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)

### Check fruitless
grep("fruitless", OGS2$Name, value = T)
fru1_geneID<-OGS2[grep("fruitless", OGS2$Name, value = F)[1],"geneID"]
## Extract eigenexons corresponding to transformer2
grep(fru1_geneID, eigenexon_evalues$eigenexonID, value=T)

## Check exon assignment file
grep(fru1_geneID, eigenexons_assignments$eigenexonID, value=T)
eigenexons_assignments[grep(fru1_geneID, eigenexons_assignments$eigenexonID, value=F),]

## Plot dynamics of gene over time

## Plot individual exons instead
fru1_exons_evalues<-nasonia_devtesova.exon.r99[grep(fru1_geneID, nasonia_devtesova.exon.r99$EXON, value=F),c(1,grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99)))]

fru1_exons_evalues<-melt(fru1_exons_evalues)
fru1_exons_evalues$value<-ifelse(fru1_exons_evalues$value>=0,fru1_exons_evalues$value,0)

fru1_exons_evalues$stage<-as.factor(str_extract(string=fru1_exons_evalues$variable, pattern="^[[:alnum:]]*"))
fru1_exons_evalues$stage<-factor(fru1_exons_evalues$stage, levels=c("emb10","emb18","lar51","pupyel","adult", "ovaries", "testes"))
fru1_exons_evalues$sex<-as.factor(ifelse(grepl("female|ovaries", fru1_exons_evalues$variable),"Female","Male"))
fru1_exons_evalues$grouper<-as.factor(paste(fru1_exons_evalues$EXON, fru1_exons_evalues$sex))
fru1_exons_evalues<-merge(fru1_exons_evalues, eigenexons_assignments[,c("exonID","eigenexonID")], by.x="EXON", by.y="exonID")

# ggplot(fru1_exons_evalues, aes(x=stage, y=value, col=sex, group=grouper))+geom_point()+geom_smooth(se=F)
ggplot(fru1_exons_evalues, aes(x=stage, y=value, col=eigenexonID, lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)

## Hclust exons
library(amap)
fru1_exons_evalues2<-as.matrix(nasonia_devtesova.exon.r99[grep(fru1_geneID, nasonia_devtesova.exon.r99$EXON, value=F),grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99))])
row.names(fru1_exons_evalues2)<-nasonia_devtesova.exon.r99[grep(fru1_geneID, nasonia_devtesova.exon.r99$EXON, value=F),"EXON"]
fru1_exons_evalues2<-apply(fru1_exons_evalues2, c(1,2), function(x){ifelse(x<=0, 0, x)})
fru1_exons_evalues2<-fru1_exons_evalues2[,-grep("testes|ovaries",colnames(fru1_exons_evalues2))]

fru1_hcluster<-hcluster(fru1_exons_evalues2, method = "correlation", link = "complete")
plot(fru1_hcluster)
fru1_cutree<-data.frame(cluster=cutree(fru1_hcluster, k=2))

fru1_exons_evalues_cutree<-merge(fru1_exons_evalues, fru1_cutree, by.x="EXON", by.y=0)
ggplot(fru1_exons_evalues_cutree, aes(x=stage, y=value, col=as.factor(cluster), lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)

### Check fruitless
grep("fruitless", OGS2$Name, value = T)
fru2_geneID<-OGS2[grep("fruitless", OGS2$Name, value = F)[2],"geneID"]
## Extract eigenexons corresponding to transformer2
grep(fru1_geneID, eigenexon_evalues$eigenexonID, value=T)

## Check exon assignment file
grep(fru2_geneID, eigenexons_assignments$eigenexonID, value=T)
eigenexons_assignments[grep(fru2_geneID, eigenexons_assignments$eigenexonID, value=F),]

## Plot dynamics of gene over time

## Plot individual exons instead
fru2_exons_evalues<-nasonia_devtesova.exon.r99[grep(fru2_geneID, nasonia_devtesova.exon.r99$EXON, value=F),c(1,grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99)))]

fru2_exons_evalues<-melt(fru2_exons_evalues)
fru2_exons_evalues$value<-ifelse(fru2_exons_evalues$value>=0,fru2_exons_evalues$value,0)

fru2_exons_evalues$stage<-as.factor(str_extract(string=fru2_exons_evalues$variable, pattern="^[[:alnum:]]*"))
fru2_exons_evalues$stage<-factor(fru2_exons_evalues$stage, levels=c("emb10","emb18","lar51","pupyel","adult", "ovaries", "testes"))
fru2_exons_evalues$sex<-as.factor(ifelse(grepl("female|ovaries", fru2_exons_evalues$variable),"Female","Male"))
fru2_exons_evalues$grouper<-as.factor(paste(fru2_exons_evalues$EXON, fru2_exons_evalues$sex))
fru2_exons_evalues<-merge(fru2_exons_evalues, eigenexons_assignments[,c("exonID","eigenexonID")], by.x="EXON", by.y="exonID")

# ggplot(fru2_exons_evalues, aes(x=stage, y=value, col=sex, group=grouper))+geom_point()+geom_smooth(se=F)
ggplot(fru2_exons_evalues, aes(x=stage, y=value, col=eigenexonID, lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)

## Hclust exons
library(amap)
fru2_exons_evalues2<-as.matrix(nasonia_devtesova.exon.r99[grep(fru2_geneID, nasonia_devtesova.exon.r99$EXON, value=F),grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99))])
row.names(fru2_exons_evalues2)<-nasonia_devtesova.exon.r99[grep(fru2_geneID, nasonia_devtesova.exon.r99$EXON, value=F),"EXON"]
fru2_exons_evalues2<-apply(fru2_exons_evalues2, c(1,2), function(x){ifelse(x<=0, 0, x)})
fru2_exons_evalues2<-fru2_exons_evalues2[,-grep("testes|ovaries",colnames(fru2_exons_evalues2))]

fru2_hcluster<-hcluster(fru2_exons_evalues2, method = "correlation", link = "complete")
plot(fru2_hcluster)
fru2_cutree<-data.frame(cluster=cutree(fru2_hcluster, h=0.9))

fru2_exons_evalues_cutree<-merge(fru2_exons_evalues, fru2_cutree, by.x="EXON", by.y=0)
ggplot(fru2_exons_evalues_cutree, aes(x=stage, y=value, col=as.factor(cluster), lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)

## OGS1.2 TRA
traOGS12_ID<-OGS2[grep("NV12464-RA",OGS2$nasviOGS1),"geneID"]
## Extract eigenexons corresponding to transformer2
grep(traOGS12_ID, eigenexon_evalues$eigenexonID, value=T)

## Check exon assignment file
grep(traOGS12_ID, eigenexons_assignments$eigenexonID, value=T)
eigenexons_assignments[grep(traOGS12_ID, eigenexons_assignments$eigenexonID, value=F),]

## Plot dynamics of gene over time

## Plot individual exons instead
traOGS12_exons_evalues<-nasonia_devtesova.exon.r99[grep(traOGS12_ID, nasonia_devtesova.exon.r99$EXON, value=F),c(1,grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99)))]

traOGS12_exons_evalues<-melt(traOGS12_exons_evalues)
traOGS12_exons_evalues$value<-ifelse(traOGS12_exons_evalues$value>=0,traOGS12_exons_evalues$value,0)

traOGS12_exons_evalues$stage<-as.factor(str_extract(string=traOGS12_exons_evalues$variable, pattern="^[[:alnum:]]*"))
traOGS12_exons_evalues$stage<-factor(traOGS12_exons_evalues$stage, levels=c("emb10","emb18","lar51","pupyel","adult", "ovaries", "testes"))
traOGS12_exons_evalues$sex<-as.factor(ifelse(grepl("female|ovaries", traOGS12_exons_evalues$variable),"Female","Male"))
traOGS12_exons_evalues$grouper<-as.factor(paste(traOGS12_exons_evalues$EXON, traOGS12_exons_evalues$sex))
traOGS12_exons_evalues<-merge(traOGS12_exons_evalues, eigenexons_assignments[,c("exonID","eigenexonID")], by.x="EXON", by.y="exonID")

# ggplot(traOGS12_exons_evalues, aes(x=stage, y=value, col=sex, group=grouper))+geom_point()+geom_smooth(se=F)
ggplot(traOGS12_exons_evalues, aes(x=stage, y=value, col=eigenexonID, lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)

## Hclust exons
library(amap)
traOGS12_exons_evalues2<-as.matrix(nasonia_devtesova.exon.r99[grep(traOGS12_ID, nasonia_devtesova.exon.r99$EXON, value=F),grep("*.[[:digit:]]$",names(nasonia_devtesova.exon.r99))])
row.names(traOGS12_exons_evalues2)<-nasonia_devtesova.exon.r99[grep(traOGS12_ID, nasonia_devtesova.exon.r99$EXON, value=F),"EXON"]
traOGS12_exons_evalues2<-apply(traOGS12_exons_evalues2, c(1,2), function(x){ifelse(x<=0, 0, x)})
traOGS12_exons_evalues2<-traOGS12_exons_evalues2[,-grep("testes|ovaries",colnames(traOGS12_exons_evalues2))]

traOGS12_hcluster<-hcluster(traOGS12_exons_evalues2, method = "correlation", link = "complete")
plot(traOGS12_hcluster)
traOGS12_cutree<-data.frame(cluster=cutree(traOGS12_hcluster, k=3))

traOGS12_exons_evalues_cutree<-merge(traOGS12_exons_evalues, traOGS12_cutree, by.x="EXON", by.y=0)
ggplot(traOGS12_exons_evalues_cutree, aes(x=stage, y=value, col=as.factor(cluster), lty=sex, group=grouper))+geom_point()+geom_smooth(se=F)+geom_smooth(aes(group=sex))+theme_bw()+facet_wrap(~EXON)
  