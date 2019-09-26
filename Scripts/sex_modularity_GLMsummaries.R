# Initialize script
date()
rm(list=ls())
newdir<-file.path(getwd(), "Output/sex_modularity_GLMsummaries")
graphicsdir<-file.path(getwd(), "Graphics/sex_modularity_GLMsummaries")
dir.create(newdir)
dir.create(graphicsdir)
library(plyr)
library(reshape)
library(ggplot2)
library(stringr)
library(vcd)
# options(na.action="na.fail")

######## NOTE ######
# Sexbias is expressed as number of female samples removed (0-3)
# Negative terms indicate that the parmeter decreases with removal of female samples (is female-specific)
# Positive terms indicate that the parameter increases with removal of female samples (is male-specific)
####################

### Load data (based on permutation analyses without resampling)
## Load changes in integration/pleiotropy
fdr_perm_dWithin <- read.csv("./Output/sex_modularity_GLMs/fdr_perm_dWithin.csv")
fdr_perm_dOut <- read.csv("./Output/sex_modularity_GLMs/fdr_perm_dOut.csv")
## Reorder stages
# dWithin
fdr_perm_dWithin$coefWithinFactornames<-str_sub(str_extract(fdr_perm_dWithin$coefWithinFactornames, pattern = "^[^:]*"), start = 6)
fdr_perm_dWithin$coefWithinFactornames<-factor(fdr_perm_dWithin$coefWithinFactornames, c("emb10","emb18", "lar51", "pupyel", "adult"))
# dOut
fdr_perm_dOut$coefOutFactornames<-str_sub(str_extract(fdr_perm_dOut$coefOutFactornames, pattern = "^[^:]*"), start = 6)
fdr_perm_dOut$coefOutFactornames<-factor(fdr_perm_dOut$coefOutFactornames, c("emb10","emb18", "lar51", "pupyel", "adult"))

## How many significant integration changes can we detect in total?
# From unadjusted lFDR
ggplot(fdr_perm_dWithin, aes(x=lFDR_individual))+geom_histogram(binwidth=0.01)+theme_bw()+scale_y_log10()+geom_vline(aes(xintercept=0.01))
table(fdr_perm_dWithin$lFDR_individual<0.001)
table(fdr_perm_dWithin$lFDR_individual<0.01) # selected threshold, expected 1 false discovery
table(fdr_perm_dWithin$lFDR_individual<0.05)
table(fdr_perm_dWithin$lFDR_individual<0.1) 
table(fdr_perm_dWithin$lFDR_individual<0.5)

prop.table(table(fdr_perm_dWithin$lFDR_individual<0.01)) # selected threshold, 7% of all contrasts is significant

# From adjusted lFDR
ggplot(fdr_perm_dWithin, aes(x=lFDR_adjusted))+geom_histogram(binwidth=0.01)+theme_bw()+scale_y_log10()+geom_vline(aes(xintercept=0.05))
table(fdr_perm_dWithin$lFDR_adjusted<0.001)
table(fdr_perm_dWithin$lFDR_adjusted<0.01) # selected threshold, expected 1 false discovery
table(fdr_perm_dWithin$lFDR_adjusted<0.05) 
table(fdr_perm_dWithin$lFDR_adjusted<0.1)

prop.table(table(fdr_perm_dWithin$lFDR_adjusted<0.01)) # selected threshold, 3% of all contrasts is significant

# What's the convergence between adjusted and unadjusted FDRs?
table(fdr_perm_dWithin$lFDR_adjusted<0.1, fdr_perm_dWithin$lFDR_individual<0.01)

## Pseudo-volcano plots
fdr_perm_dWithin$significantID<-fdr_perm_dWithin$.id
fdr_perm_dWithin$significantID[which(fdr_perm_dWithin$lFDR_individual>0.01&fdr_perm_dWithin$lFDR_adjusted>0.01)]<-NA
fdr_perm_dWithin$significantID<-droplevels(fdr_perm_dWithin$significantID)

ggplot(fdr_perm_dWithin, aes(y=lFDR_individual, x=Estimate, col=significantID))+geom_point()+theme_bw()+scale_y_log10()+geom_vline(aes(yintercept=0))+geom_hline(aes(yintercept=0.01))+facet_wrap(~coefWithinFactornames, nrow = 1)
ggplot(fdr_perm_dWithin, aes(y=lFDR_adjusted, x=Estimate, col=significantID))+geom_point()+theme_bw()+scale_y_log10()+geom_vline(aes(yintercept=0))+geom_hline(aes(yintercept=0.01))+facet_wrap(~coefWithinFactornames, nrow = 1)

# How many clusters show sex-biased changes in integration in how many stages?
table(table(fdr_perm_dWithin[which(fdr_perm_dWithin$lFDR_individual<0.01),]$.id))
table(table(fdr_perm_dWithin[which(fdr_perm_dWithin$lFDR_adjusted<0.01),]$.id))
# Only 1 cluster (thistle3) has significant differences in dWithin more than one stage, and only if we use unadjusted lFDR

# Tabulate clusters with significant sex-biased changes in integration, using adjusted lFDR
sign_dWithin<-fdr_perm_dWithin[which(fdr_perm_dWithin$lFDR_adjusted<0.01),]
# Tabulate numer of clusters per stage and direction
table(sign_dWithin$coefWithinFactornames,sign(sign_dWithin$Estimate))
margin.table(table(sign_dWithin$coefWithinFactornames,sign(sign_dWithin$Estimate)), margin = 2)
margin.table(table(sign_dWithin$coefWithinFactornames,sign(sign_dWithin$Estimate)), margin = 1)
# Greatest number of clusters in embryo18&pupa (surprising) lowest in lar51 (surprising)
# Male-specific integration (positive coefficients) predominates in all stages but emb10

# Save names of clusters with significant changes in integration, including stage and direction annotation
sign_dWithin2<-sign_dWithin[,c("coefWithinFactornames",".id","Estimate","lFDR_adjusted")]
sign_dWithin2$Sex<-as.factor(ifelse(sign_dWithin2$Estimate>0,"Male","Female"))
names(sign_dWithin2)<-c("Stage","clusterID","Estimate","lFDR_adjusted","Sex")
write.csv(sign_dWithin2, file=file.path(newdir, "dWithin_signClusters.csv"))

## How many significant constraint changes can we detect in total?
# From unadjusted lFDR
ggplot(fdr_perm_dOut, aes(x=lFDR_individual))+geom_histogram(binwidth=0.01)+theme_bw()+scale_y_log10()+geom_vline(aes(xintercept=0.05))
table(fdr_perm_dOut$lFDR_individual<0.001)
table(fdr_perm_dOut$lFDR_individual<0.01) 
table(fdr_perm_dOut$lFDR_individual<0.05) # selected threshold, expected 2 false discoveries
table(fdr_perm_dOut$lFDR_individual<0.1)

prop.table(table(fdr_perm_dOut$lFDR_individual<0.05)) # selected threshold, expected 3 false discoveries

# From adjusted lFDR
ggplot(fdr_perm_dOut, aes(x=lFDR_adjusted))+geom_histogram(binwidth=0.01)+theme_bw()+scale_y_log10()+geom_vline(aes(xintercept=0.05))
table(fdr_perm_dOut$lFDR_adjusted<0.001)
table(fdr_perm_dOut$lFDR_adjusted<0.01)
table(fdr_perm_dOut$lFDR_adjusted<0.05)
table(fdr_perm_dOut$lFDR_adjusted<0.1) # selected threshold, expected 2 false discoveries

prop.table(table(fdr_perm_dOut$lFDR_adjusted<0.1)) # selected threshold, expected 2 false discoveries

# What's the convergence between adjusted and unadjusted FDRs?
table(fdr_perm_dOut$lFDR_adjusted<0.1, fdr_perm_dOut$lFDR_individual<0.05)

## Pseudo-volcano plots
fdr_perm_dOut$significantID<-fdr_perm_dOut$.id
fdr_perm_dOut$significantID[which(fdr_perm_dOut$lFDR_individual>0.05&fdr_perm_dOut$lFDR_adjusted>0.1)]<-NA
fdr_perm_dOut$significantID<-droplevels(fdr_perm_dOut$significantID)

ggplot(fdr_perm_dOut, aes(y=lFDR_individual, x=Estimate, col=significantID))+geom_point()+theme_bw()+scale_y_log10()+geom_vline(aes(yintercept=0))+geom_hline(aes(yintercept=0.05))+facet_wrap(~coefOutFactornames, nrow = 1)
ggplot(fdr_perm_dOut, aes(y=lFDR_adjusted, x=Estimate, col=significantID))+geom_point()+theme_bw()+scale_y_log10()+geom_vline(aes(yintercept=0))+geom_hline(aes(yintercept=0.1))+facet_wrap(~coefOutFactornames, nrow = 1)

# How many clusters show sex-biased changes in constraint in how many stages?
table(table(fdr_perm_dOut[which(fdr_perm_dOut$lFDR_individual<0.05),]$.id))
table(table(fdr_perm_dOut[which(fdr_perm_dOut$lFDR_adjusted<0.1),]$.id))
# Only 3 clusters have significant differences in more than one stage, and only if we use unadjusted lFDR

# Tabulate clusters with significant sex-biased changes in constraint, using adjusted lFDR
sign_dOut<-fdr_perm_dOut[which(fdr_perm_dOut$lFDR_adjusted<0.1),]
# Tabulate numer of clusters per stage and direction
table(sign_dOut$coefOutFactornames,sign(sign_dOut$Estimate))
# Greatest number of clusters in adult (unsurprising) lowest in lar51 (surprising)
# Female-specific constraint (negative coefficients) predominates in all stages but emb18 and adult

### Test convergence between integration and pleiotropy
sign_clusters<-merge(sign_dOut, sign_dWithin, all=T, by.x=c("coefOutFactornames", ".id"), by.y=c("coefWithinFactornames", ".id"), suffixes = c("_dOut","_dWithin"))
table(sign(sign_clusters$Estimate_dOut),sign(sign_clusters$Estimate_dWithin))
# No cluster is significant in both

# Create compacted dataset with both dWithin and dOut 
names(fdr_perm_dWithin)[2]<-"Factornames"
names(fdr_perm_dOut)[2]<-"Factornames"
fdr_perm<-merge(fdr_perm_dOut[,-1], fdr_perm_dWithin[,-1], by=c(".id","Factornames"), suffixes = c("_dOut","_dWithin"))

# Which clusters are significant for changes in dWithin by stage?
dWithin_sign<-fdr_perm[which(fdr_perm$lFDR_adjusted_dWithin<0.05),c(".id","Estimate_dWithin" , "Factornames")]
dWithin_sign<-dWithin_sign[order(dWithin_sign$Factornames),]
write.csv(dWithin_sign, file=file.path(newdir, "dWithin_sign.csv"))
# Which clusters are significant for changes in dOut by stage?
dOut_sign<-fdr_perm[which(fdr_perm$lFDR_adjusted_dOut<0.1),c(".id","Estimate_dOut" , "Factornames")]
dOut_sign<-dOut_sign[order(dOut_sign$Factornames),]
write.csv(dOut_sign, file=file.path(newdir, "dOut_sign.csv"))

# Plot effect sizes of all clusters
ggplot(fdr_perm, aes(x=Estimate_dOut, y=Estimate_dWithin, xmax=Estimate_dOut+Std..Error_dOut, xmin=Estimate_dOut-Std..Error_dOut, ymax=Estimate_dWithin+Std..Error_dWithin, ymin=Estimate_dWithin-Std..Error_dWithin, col=Factornames))+geom_point(alpha="0.3")+theme_bw()+geom_errorbar(alpha="0.3")+geom_errorbarh(alpha="0.3")+geom_density2d(group="1")

ggplot(fdr_perm, aes(x=Estimate_dOut+1, y=Estimate_dWithin+1, xmax=Estimate_dOut+Std..Error_dOut+1, xmin=Estimate_dOut-Std..Error_dOut+1, ymax=Estimate_dWithin+Std..Error_dWithin+1, ymin=Estimate_dWithin-Std..Error_dWithin+1, col=Factornames))+geom_point()+theme_bw()+geom_errorbar()+geom_errorbarh()+scale_x_log10()+scale_y_log10()

# What is the correlation between correlation and parcellation in all clusters?
cor(fdr_perm$Estimate_dOut, fdr_perm$Estimate_dWithin)
cor(fdr_perm$Estimate_dOut, fdr_perm$Estimate_dWithin, method = "spearman")
cor(fdr_perm$Estimate_dOut, fdr_perm$Estimate_dWithin, method = "kendall")

cor_bystage<-ddply(.data = fdr_perm, .variables = .(Factornames), .fun = summarize, pearson=cor(Estimate_dOut, Estimate_dWithin), spearman=cor(Estimate_dOut, Estimate_dWithin, method = "spearman"), kendall=cor(Estimate_dOut, Estimate_dWithin, method = "kendall"))

ddply(.data = fdr_perm, .variables = .(Factornames), .fun = summarize, pearson=cor.test(Estimate_dOut, Estimate_dWithin)$p.value, spearman=cor.test(Estimate_dOut, Estimate_dWithin, method = "spearman")$p.value, kendall=cor.test(Estimate_dOut, Estimate_dWithin, method = "kendall")$p.value)

cor_bycluster<-ddply(.data = fdr_perm, .variables = .(.id), .fun = summarize, pearson=cor(Estimate_dOut, Estimate_dWithin), spearman=cor(Estimate_dOut, Estimate_dWithin, method = "spearman"), kendall=cor(Estimate_dOut, Estimate_dWithin, method = "kendall"))

ddply(.data = fdr_perm, .variables = .(.id), .fun = summarize, pearson=cor.test(Estimate_dOut, Estimate_dWithin)$p.value, spearman=cor.test(Estimate_dOut, Estimate_dWithin, method = "spearman")$p.value, kendall=cor.test(Estimate_dOut, Estimate_dWithin, method = "kendall")$p.value)

# What is the overall level of integration/parcellation of each cluster?
load(file="./Output/sex_modularity_GLMs/dWithinGLMs")
load(file="./Output/sex_modularity_GLMs/dOutGLMs")
library(lattice)

densityplot(unlist(lapply(dWithinGLMs, function(x){median(x$coefficients[1:5])})))
densityplot(unlist(lapply(dOutGLMs, function(x){median(x$coefficients[1:5])})))

# What are the average changes?
ggplot(data = fdr_perm, aes(x=Estimate_dOut,col=Factornames))+geom_density()+theme_bw()
ggplot(data = fdr_perm, aes(x=Estimate_dWithin,col=Factornames))+geom_density()+theme_bw()

ggplot(data = fdr_perm, aes(x=Estimate_dOut, y=Estimate_dWithin))+geom_density2d()+theme_bw()+facet_wrap(~Factornames)

ggplot(fdr_perm, aes(x=Estimate_dOut, y=Estimate_dWithin, xmax=Estimate_dOut+Std..Error_dOut, xmin=Estimate_dOut-Std..Error_dOut, ymax=Estimate_dWithin+Std..Error_dWithin, ymin=Estimate_dWithin-Std..Error_dWithin, alpha=2-(lFDR_individual_dWithin+lFDR_individual_dOut)))+geom_point()+theme_bw()+geom_errorbar()+geom_errorbarh()+geom_density2d()+facet_wrap(~Factornames)+theme(legend.position="bottom")

# fdr_perm$minFdr<-pmin(fdr_perm$lFDR_adjusted_dOut, fdr_perm$lFDR_adjusted_dWithin)
fdr_perm2<-fdr_perm
fdr_perm2$Factornames<-revalue(fdr_perm2$Factornames, c("emb10"="embryo 10", "emb18" = "embryo 18", "lar51" = "larva 51", "pupyel" = "pupa", "adult" = "adult"))
pdf(file = file.path(graphicsdir,"clusterbiasedchanges_2.pdf"), width = 8.605, height = 5.378, useDingbats=FALSE)
ggplot(fdr_perm2, aes(x=Estimate_dOut, y=Estimate_dWithin, xmax=Estimate_dOut+Std..Error_dOut, xmin=Estimate_dOut-Std..Error_dOut, ymax=Estimate_dWithin+Std..Error_dWithin, ymin=Estimate_dWithin-Std..Error_dWithin, alpha=1-pmin(lFDR_individual_dWithin,lFDR_individual_dOut)))+geom_point()+theme_bw()+geom_errorbar()+geom_errorbarh()+geom_density2d()+facet_wrap(~Factornames, nrow=2)+scale_alpha_continuous(name="Significance Score")+scale_x_continuous(name="Changes in Constraint \n per male transcriptome added")+scale_y_continuous(name="Changes in Integration \n per male transcriptome added")+coord_fixed()+ggtitle(label = "Sex-Specific changes in Network Parameters\n across Development\n")+theme(legend.position=c(0.8,0.2))
dev.off()


# Plot cross tabulation of clusters with significant changes in each stage
mosaic(formula = ~(Estimate_dOut>0)+(Estimate_dWithin>0)+(lFDR_adjusted_dOut<0.05)+(lFDR_adjusted_dWithin<0.05)|Factornames, data = fdr_perm)
mosaic(formula = ~(lFDR_adjusted_dOut<0.05)+(Estimate_dOut>0)|Factornames, data = fdr_perm)
mosaic(formula = ~(lFDR_adjusted_dWithin<0.05)+(Estimate_dWithin>0)|Factornames, data = fdr_perm)

# Load cluster annotation (size, count of splicing/expression nodes, count of paralog/ortholog/NA nodes)



######## Outdated code snippets

### Test convergence between integration and pleiotropy
# sign_clusters<-merge(sign_dOut, sign_dWithin, all=T, by.x=c("coefOutFactornames", ".id"), by.y=c("coefWithinFactornames", ".id"), suffixes = c("_dOut","_dWithin"))
# table(sign(sign_clusters$Estimate_dOut),sign(sign_clusters$Estimate_dWithin))
# No contrasts are shared between the two sets
# sign_clusters<-sign_clusters[which(is.na(sign_clusters$Estimate_dOut)==F&is.na(sign_clusters$Estimate_dWithin)==F),]
# # 9 out of 11 are male specific, 2 female specific, all cases have convergent signs
# table(sign(sign_clusters$Estimate_dOut), sign(sign_clusters$Estimate_dWithin))
# # Plot effect sizes of clusters with significant effect in both
# ggplot(sign_clusters, aes(x=Estimate_dOut, y=Estimate_dWithin, xmax=Estimate_dOut+Std..Error_dOut, xmin=Estimate_dOut-Std..Error_dOut, ymax=Estimate_dWithin+Std..Error_dWithin, ymin=Estimate_dWithin-Std..Error_dWithin, col=coefOutFactornames))+geom_point()+theme_bw()+geom_errorbar()+geom_errorbarh()
