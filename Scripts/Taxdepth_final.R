# Compare taxonomic depth of clusters
# Initialize script
library(reshape2)
library(plyr)
library(stringr)
library(vcd)
library(ggplot2)
library(lattice)
library(fdrtool)
library(MuMIn)
options(na.action="na.fail")
newdir<-file.path(getwd(),"Output/Taxdepth_final")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/Taxdepth_final")
dir.create(graphdir)

# Load dataset
read.csv(file = "./Output/Results_compiler/moduledata_short.csv")[,-1]

# Compare proprtion of nodes belonging to each stratum in each cluster, split by sexbias
taxdepth_boxplot<-recast(moduledata_short[,c("clusterID","signif","taxaNasonia_p","taxaHymenoptera_p","taxaInsecta_p","taxaArthropoda_p","taxaMetazoa_p")], formula = clusterID+signif+variable~.)
names(taxdepth_boxplot)<-c("clusterID","signif","stratum","proportion")
ggplot(taxdepth_boxplot, aes(x=signif, y=proportion))+geom_boxplot(notch = T, varwidth = T)+facet_grid(.~stratum)+theme_bw()

# Use linear model
taxdepth_boxplot$proportion<-ifelse(taxdepth_boxplot$proportion==0, 0.001, taxdepth_boxplot$proportion)
lm_taxa<-glm(formula = proportion~-1+stratum*signif, data = taxdepth_boxplot, family = Gamma)
plot(lm_taxa)
model.sel(dredge(lm_taxa))
lm_taxa_str<-glm(formula = proportion~-1+stratum, data = taxdepth_boxplot, family = Gamma)
plot(lm_taxa_str)
summary(lm_taxa_str)

# Restrict to clusters with non-zero proportions (rationale, clusters are taxonomically consistent, good hypothesis to check)
taxdepth_boxplot<-taxdepth_boxplot[which(taxdepth_boxplot$proportion>=0.01),]
lm_taxa<-glm(formula = proportion~-1+stratum*signif, data = taxdepth_boxplot, family = Gamma)
plot(lm_taxa)
model.sel(dredge(lm_taxa))
lm_taxa_str<-glm(formula = proportion~-1+stratum, data = taxdepth_boxplot, family = Gamma)
plot(lm_taxa_str)
summary(lm_taxa_str)

# Proportions are consistent across clusters, regardless of sexbias

# Compare proportion of sexbiased vs non-sexbiased nodes in each stratum
taxdepth_barplot<-recast(moduledata_short[,c("clusterID","signif","taxaNasonia_n","taxaHymenoptera_n","taxaInsecta_n","taxaArthropoda_n","taxaMetazoa_n")], formula = signif+variable~., fun.aggregate=sum)
names(taxdepth_barplot)<-c("signif","stratum","number")
# taxdepth_barplot<-ddply(taxdepth_barplot, .variables = .(signif, stratum), .fun = summarize, value=sum(value))
ggplot(taxdepth_barplot, aes(x=stratum, col=signif, y=number))+geom_point()+theme_bw()
taxdepth_barplot<-recast(taxdepth_barplot, formula = stratum~signif)
taxdepth_barplot<-rbind(taxdepth_barplot, c(NA,sum(taxdepth_barplot$"FALSE"),sum(taxdepth_barplot$"TRUE")))
taxdepth_barplot$stratum<-as.factor(c(as.character(taxdepth_barplot$stratum)[-6],"allnodes_n"))
taxdepth_barplot$prop_sexb<-taxdepth_barplot$"TRUE"/(taxdepth_barplot$"FALSE"+taxdepth_barplot$"TRUE")
taxdepth_barplot$prop_unb<-taxdepth_barplot$"FALSE"/(taxdepth_barplot$"FALSE"+taxdepth_barplot$"TRUE")
taxdepth_barplot<-melt(taxdepth_barplot[,c("stratum","prop_sexb","prop_unb")], variable.name = "sexbias")
ggplot(taxdepth_barplot, aes(x=stratum, col=sexbias, y=value))+geom_point()+geom_line(aes(group=sexbias))+theme_bw()
ggplot(taxdepth_barplot, aes(x=stratum, fill=sexbias, y=value))+geom_bar(stat="identity")+theme_bw()

taxdepth_mosaicplot<-recast(moduledata_short[,c("clusterID","signif","taxaNasonia_n","taxaHymenoptera_n","taxaInsecta_n","taxaArthropoda_n","taxaMetazoa_n")], 
                            formula = variable~signif, fun.aggregate=sum)
row.names(taxdepth_mosaicplot)<-taxdepth_mosaicplot$variable
taxdepth_mosaicplot<-taxdepth_mosaicplot[,-1]
mosaicplot(x = taxdepth_mosaicplot, shade = T)

## Conclusions
# Taxonomic composition of clusters is similar between sebiased and non-sexbiased ones
# Hymenopteran nodes belong to more sexbiased clusters than expected
