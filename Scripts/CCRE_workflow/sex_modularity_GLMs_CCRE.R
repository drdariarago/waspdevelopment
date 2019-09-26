### Detect sex:stage related changes in module density using GLMs (merges sex_modularity_GLMs and sex_modularity_GLMsummaries)
date()
rm(list=ls())
library(plyr)
library(stringr)
library(WGCNA)
library(ggplot2)
library(fdrtool)
allowWGCNAThreads()
library(boot)
library(reshape2)
source(file = "./Scripts/multiplot.R")
# initialize output path
newdir<-file.path(getwd(), "Output/sex_modularity_GLMs_CCRE")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/sex_modularity_GLMs_CCRE")
dir.create(graphdir)

## import permuted cluster values
modulardensities <- read.csv("./Output/sex_modularity_biweight_CCRE/sex_modulewise_modulardensities.csv")
# reorder stages
modulardensities$stage<-factor(modulardensities$stage, c("emb10","emb18", "lar51", "pupyel", "adult"))

## QC
# Are there any densities greater than 1 or lesser than 0?
summary(apply(X = modulardensities[,c("dWithin", "dOut", "dMain")], MARGIN = 1, FUN = function(x){any(x>1)}))
summary(apply(X = modulardensities[,c("dWithin", "dOut", "dMain")], MARGIN = 1, FUN = function(x){any(x<0)}))
summary(apply(X = modulardensities[,c("dWithin", "dOut", "dMain")], MARGIN = 1, FUN = function(x){any(x==0)}))
# Which ones?
# modulardensities[which(apply(X = modulardensities[,c("dWithin", "dOut", "dMain")], MARGIN = 1, FUN = function(x){any(x>1)})),]
# modulardensities[which(apply(X = modulardensities[,c("dWithin", "dOut", "dMain")], MARGIN = 1, FUN = function(x){any(x<0)})),]

# is there any variance/mean correlation?
# ggplot(data=ddply(modulardensities, .(stage, sexbias, clusterID), summarize, dWithinM=mean(dWithin), dWithinS=sd(dWithin)), aes(x=dWithinM,y=dWithinS, col=clusterID))+geom_point()+geom_smooth()
pdf(file = file.path(graphdir, "mean_variance.pdf"))
with(ddply(modulardensities, .(stage, sexbias, clusterID), summarize, dWithinM=mean(dWithin), dWithinS=sd(dWithin)), cor.test(dWithinM, dWithinS, use = "pair"))
with(ddply(modulardensities, .(stage, sexbias, clusterID), summarize, dOutM=mean(dOut), dOutS=sd(dOut)), cor.test(dOutM, dOutS, use = "pair"))
dev.off()

## Remove cluster grey (has density zero in all permutations)
modulardensities <- modulardensities[-which(modulardensities$clusterID=="grey"),]

## Test via linear models
# For changes in integration coefficient
dWithinGLMs<-dlply(.data = modulardensities, 
                   .variables = .(clusterID),
                   .fun = glm, 
                   formula = dWithin~stage-1+stage:sexbias+dMain+dOut, 
                   family = Gamma(link="log"))
# Store diagnostic plots
pdf(file = file.path(graphdir, "dWithinGLM_diagnostics.pdf"))
par(mfrow=c(2,2))
lapply(X = dWithinGLMs, FUN = plot)
dev.off()
# And save model list
save(dWithinGLMs, file = file.path(newdir, "dWithinGLMs"))

# For changes in pleiotropy coefficient
dOutGLMs<-dlply(.data = modulardensities, 
                .variables = .(clusterID),
                .fun = glm, 
                formula = dOut~stage-1+stage:sexbias+dMain+dWithin, 
                family = Gamma(link="logit"))
# Store diagnostic plots
pdf(file = file.path(graphdir, "dOutGLM_diagnostics.pdf"))
par(mfrow=c(2,2))
lapply(X = dOutGLMs, FUN = plot)
dev.off()
# And save model list
save(dOutGLMs, file = file.path(newdir, "dOutGLMs"))

## Extract coefficients and p-values for integration coefficients
coefWithin<-llply(.data = dWithinGLMs, function(x){summary(x)$coefficient})
coefWithinFactornames<-unlist(llply(.data = coefWithin, function(x){row.names(x)}))
coefWithin<-ldply(coefWithin) 
coefWithin<-cbind(coefWithinFactornames,coefWithin)
coefWithin$qval<-fdrtool(coefWithin$"Pr(>|t|)", statistic = "pvalue")$qval

write.csv(coefWithin, file=file.path(newdir, "integration_coefs.csv"))

## Extract coefficients and p-values for pleiotropy coefficients
coefOut<-llply(.data = dOutGLMs, function(x){summary(x)$coefficient})
coefOutFactornames<-unlist(llply(.data = coefOut, function(x){row.names(x)}))
coefOut<-ldply(coefOut) 
coefOut<-cbind(coefOutFactornames,coefOut)
coefOut$qval<-fdrtool(coefOut$"Pr(>|t|)", statistic = "pvalue")$qval

write.csv(coefOut, file=file.path(newdir, "pleiotropy_coefs.csv"))

## Validate via bootstrapping (without replacement)

# Define number of bootstrap replicates
nperm<-1000
# Reduce dataset for testing
# modulardensities<-modulardensities[grep("white|red|orange",modulardensities$clusterID),]

#### Testing using bootstrapping, dWithin version without replacement

# Create index of permutations (one permutation per row)
indexes<-replicate(
  n = nperm,
  sample(x = 1:table(modulardensities$clusterID)[1], size = table(modulardensities$clusterID)[1], replace = F)
)
# Apply function to each module
dWithin_perm<-apply(X = indexes, MARGIN = 2, function(indexes){
  ddply(.data = modulardensities, 
        .variables = .(clusterID),
        .fun = function(cluster){
          newdata<-cluster
          newdata[,"sexbias"]<-newdata[indexes,"sexbias"] # randomly reassign sexbias values
          model<-glm(dWithin~stage-1+stage:sexbias+dMain+dOut, data = newdata, family=Gamma(link="logit"))
          summary(model)$coefficients[c(8:12),"Pr(>|t|)"] # save p-value estimates of stage:sex interactions
        }
  )
})

# Add observed p-values as last data.frame in list (will become last row of each data.frame after splitting by cluster)
dWithin_perm[[length(dWithin_perm)+1]]<-ddply(.data = modulardensities, 
                                              .variables = .(clusterID),
                                              .fun = function(cluster){
                                                model<-glm(dWithin~stage-1+stage:sexbias+dMain+dOut, data = cluster, family=Gamma(link="logit"))
                                                summary(model)$coefficients[c(8:12),"Pr(>|t|)"]
                                              }
)

# Split list of permutations by module
dWithin_perm<-dlply(ldply(dWithin_perm), .(clusterID))

# Calculate local FDR of p-values compared to null distribution of reshuffled labels
pdf(file = file.path(graphdir, "dWithin_perm_cluster_fdr.pdf"))
dWithin_perm_q.val_1<-ldply(lapply(X = dWithin_perm, FUN = function(cluster){
  sapply(X = 2:6, function(index){
    fdrtool(cluster[,index], statistic = "pvalue", plot = T)$lfdr[nrow(cluster)]
  })
}))
dev.off()
colnames(dWithin_perm_q.val_1)<-names(dWithin_perm[[1]])

# Adjust lFDR values for multiple hypothesis testing
pdf(file = file.path(graphdir, "dWithin_perm_Fdr.pdf"))
dWithin_perm_q.val_2<-matrix(fdrtool(c(unlist(dWithin_perm_q.val_1[,-1])), statistic = "pvalue", plot = T)$qval,
                             ncol = 5, 
                             dimnames = list(
                               c(dWithin_perm_q.val_1$clusterID),
                               c(colnames(dWithin_perm_q.val_1)[-1]))
)
dev.off()

# Attach to factor estimates
dWithin_perm_q.val_1<-melt(dWithin_perm_q.val_1)
names(dWithin_perm_q.val_1)<-c(".id","coefWithinFactornames", "lFDR_individual")
dWithin_perm_q.val_2<-melt(dWithin_perm_q.val_2)
names(dWithin_perm_q.val_2)<-c(".id","coefWithinFactornames", "lFDR_adjusted")
# factor_q.val_2$coefOutFactornames<-str_extract(factor_q.val_2$coefOutFactornames, pattern = "^[^:]*")
fdr_perm_dWithin<-merge(coefWithin, dWithin_perm_q.val_1)
fdr_perm_dWithin<-merge(fdr_perm_dWithin, dWithin_perm_q.val_2)

# And save as csv
write.csv(fdr_perm_dWithin, file=file.path(newdir, "fdr_perm_dWithin.csv"))

#### Testing using bootstrapping, dOut version without replacement
# Apply function to each module
dOut_perm<-apply(X = indexes, MARGIN = 2, function(indexes){
  ddply(.data = modulardensities, 
        .variables = .(clusterID),
        .fun = function(cluster){
          newdata<-cluster
          newdata[,"sexbias"]<-newdata[indexes,"sexbias"] # randomly reassign sexbias values
          model<-glm(dOut~stage-1+stage:sexbias+dMain+dWithin, data = newdata, family=Gamma(link="logit"))
          summary(model)$coefficients[c(8:12),"Pr(>|t|)"] # save p-value estimates of stage:sex interactions
        }
  )
})

# Add observed p-values as last data.frame in list (will become last row of each data.frame after splitting by cluster)
dOut_perm[[length(dOut_perm)+1]]<-ddply(.data = modulardensities, 
                                        .variables = .(clusterID),
                                        .fun = function(cluster){
                                          model<-glm(dOut~stage-1+stage:sexbias+dMain+dWithin, data = cluster, family=Gamma(link="logit"))
                                          summary(model)$coefficients[c(8:12),"Pr(>|t|)"]
                                        }
)

# Split list of permutations by module
dOut_perm<-dlply(ldply(dOut_perm), .(clusterID))

# Calculate local FDR of p-values compared to null distribution of reshuffled labels
pdf(file = file.path(graphdir, "dOut_perm_cluster_fdr.pdf"))
dOut_perm_q.val_1<-ldply(lapply(X = dOut_perm, FUN = function(cluster){
  sapply(X = 2:6, function(index){
    fdrtool(cluster[,index], statistic = "pvalue", plot = T)$lfdr[nrow(cluster)]
  })
}))
colnames(dOut_perm_q.val_1)<-names(dOut_perm[[1]])
dev.off()

# Adjust lFDR values for multiple hypothesis testing
pdf(file = file.path(graphdir, "dOut_perm_Fdr.pdf"))
dOut_perm_q.val_2<-matrix(fdrtool(c(unlist(dOut_perm_q.val_1[,-1])), statistic = "pvalue", plot = T)$qval,
                          ncol = 5, 
                          dimnames = list(
                            c(dOut_perm_q.val_1$clusterID),
                            c(colnames(dOut_perm_q.val_1)[-1]))
)
dev.off()

# Attach to factor estimates
dOut_perm_q.val_1<-melt(dOut_perm_q.val_1)
names(dOut_perm_q.val_1)<-c(".id","coefOutFactornames", "lFDR_individual")
dOut_perm_q.val_2<-melt(dOut_perm_q.val_2)
names(dOut_perm_q.val_2)<-c(".id","coefOutFactornames", "lFDR_adjusted")
# factor_q.val_2$coefOutFactornames<-str_extract(factor_q.val_2$coefOutFactornames, pattern = "^[^:]*")
fdr_perm_dOut<-merge(coefOut, dOut_perm_q.val_1)
fdr_perm_dOut<-merge(fdr_perm_dOut, dOut_perm_q.val_2)

# And save as csv
write.csv(fdr_perm_dOut, file=file.path(newdir, "fdr_perm_dOut.csv"))

# Save workspace
save.image(file = file.path(newdir, "GLM_workspace.RData"))

### Annotate significant contrasts and their direction
######## NOTE ######
# Sexbias is expressed as number of female samples removed (0-3)
# Negative terms indicate that the parmeter decreases with removal of female samples (is female-specific)
# Positive terms indicate that the parameter increases with removal of female samples (is male-specific)
####################

# load(file = "./Output/sex_modularity_GLMs_CCRE/GLM_workspace.RData")

## Reorder stages
# dWithin
fdr_perm_dWithin$coefWithinFactornames<-str_sub(str_extract(fdr_perm_dWithin$coefWithinFactornames, pattern = "^[^:]*"), start = 6)
fdr_perm_dWithin$coefWithinFactornames<-factor(fdr_perm_dWithin$coefWithinFactornames, c("emb10","emb18", "lar51", "pupyel", "adult"))
# dOut
fdr_perm_dOut$coefOutFactornames<-str_sub(str_extract(fdr_perm_dOut$coefOutFactornames, pattern = "^[^:]*"), start = 6)
fdr_perm_dOut$coefOutFactornames<-factor(fdr_perm_dOut$coefOutFactornames, c("emb10","emb18", "lar51", "pupyel", "adult"))

## How many significant integration changes can we detect in total?
# From unadjusted lFDR
rbind(table(fdr_perm_dWithin$lFDR_individual<0.001),
      table(fdr_perm_dWithin$lFDR_individual<0.005),
      table(fdr_perm_dWithin$lFDR_individual<0.01),
      table(fdr_perm_dWithin$lFDR_individual<0.05), # selected threshold, expected >2 false discoveries
      table(fdr_perm_dWithin$lFDR_individual<0.1), 
      table(fdr_perm_dWithin$lFDR_individual<0.5))
prop.table(table(fdr_perm_dWithin$lFDR_individual<0.05)) # selected threshold, 7% of all contrasts is significant

# From adjusted lFDR
rbind(table(fdr_perm_dWithin$lFDR_adjusted<0.001),
      table(fdr_perm_dWithin$lFDR_adjusted<0.005),
      table(fdr_perm_dWithin$lFDR_adjusted<0.01), 
      table(fdr_perm_dWithin$lFDR_adjusted<0.05), 
      table(fdr_perm_dWithin$lFDR_adjusted<0.1), # selected threshold, expected >2 false discoveries
      table(fdr_perm_dWithin$lFDR_adjusted<0.5)) 
prop.table(table(fdr_perm_dWithin$lFDR_adjusted<0.1)) # selected threshold, ~3% of all contrasts is significant

# What's the convergence between adjusted and unadjusted FDRs (dWithin)?
table(fdr_perm_dWithin$lFDR_adjusted<0.1, fdr_perm_dWithin$lFDR_individual<0.05, dnn = c("lFDR_adj","lFDR_ind"))
# Perfect match, but adjusted loses 17 contrasts

## Pseudo-volcano plots
fdr_perm_dWithin$significantID<-fdr_perm_dWithin$.id
fdr_perm_dWithin$significantID[which(fdr_perm_dWithin$lFDR_individual>0.05&fdr_perm_dWithin$lFDR_adjusted>0.1)]<-NA
multiplot(cols = 2,
  ggplot(fdr_perm_dWithin, aes(y=lFDR_individual, x=Estimate, col=significantID))+geom_point()+theme_bw()+scale_y_log10(limits=c(1e-12,1))+geom_vline(aes(xintercept=0))+geom_hline(aes(yintercept=0.05))+facet_wrap(~coefWithinFactornames, nrow = 1),
  ggplot(fdr_perm_dWithin, aes(y=lFDR_adjusted, x=Estimate, col=significantID))+geom_point()+theme_bw()+scale_y_log10(limits=c(1e-12,1))+geom_vline(aes(xintercept=0))+geom_hline(aes(yintercept=0.1))+facet_wrap(~coefWithinFactornames, nrow = 1)
)
# How many clusters show sex-biased changes in integration in how many stages?
table(table(fdr_perm_dWithin[which(fdr_perm_dWithin$lFDR_individual<0.05),]$.id))
table(table(fdr_perm_dWithin[which(fdr_perm_dWithin$lFDR_adjusted<0.1),]$.id))
# One, only with unadjusted

# Tabulate clusters with significant sex-biased changes in integration, using adjusted lFDR
sign_dWithin<-fdr_perm_dWithin[which(fdr_perm_dWithin$lFDR_adjusted<0.1),]
# Tabulate numer of clusters per stage and direction
table(sign_dWithin$coefWithinFactornames,sign(sign_dWithin$Estimate))
margin.table(table(sign_dWithin$coefWithinFactornames,sign(sign_dWithin$Estimate)), margin = 2)
margin.table(table(sign_dWithin$coefWithinFactornames,sign(sign_dWithin$Estimate)), margin = 1)
# Greatest number of DI clusters in adult (unsurprising), DI present in early and late embryo but not larva and pupa (surprising!)
# Sexes balanced but not evenly partitioned between stages: female DI prevails in early embryo, male in late embryo, balanced in adults

# Save names of clusters with significant changes in integration, including stage and direction annotation
sign_dWithin2<-sign_dWithin[,c("coefWithinFactornames",".id","Estimate","lFDR_adjusted")]
sign_dWithin2$Sex<-as.factor(ifelse(sign_dWithin2$Estimate>0,"Male","Female"))
names(sign_dWithin2)<-c("Stage","clusterID","Estimate","lFDR_adjusted","Sex")
write.csv(sign_dWithin2, file = file.path(newdir, "dWithin_signClusters.csv"))

## How many significant constraint changes can we detect in total?
# From unadjusted lFDR
rbind(table(fdr_perm_dOut$lFDR_individual<0.001),
      table(fdr_perm_dOut$lFDR_individual<0.005),
      table(fdr_perm_dOut$lFDR_individual<0.01),
      table(fdr_perm_dOut$lFDR_individual<0.05), 
      table(fdr_perm_dOut$lFDR_individual<0.1), # selected threshold, expected 3 false discoveries
      table(fdr_perm_dOut$lFDR_individual<0.5)
)
prop.table(table(fdr_perm_dOut$lFDR_individual<0.1)) # selected threshold, 8% total clusters DC

# From adjusted lFDR
rbind(table(fdr_perm_dOut$lFDR_adjusted<0.001),
      table(fdr_perm_dOut$lFDR_adjusted<0.005),
      table(fdr_perm_dOut$lFDR_adjusted<0.01),
      table(fdr_perm_dOut$lFDR_adjusted<0.05),
      table(fdr_perm_dOut$lFDR_adjusted<0.1), # selected threshold, expected 1 false discovery
      table(fdr_perm_dOut$lFDR_adjusted<0.5)
)
prop.table(table(fdr_perm_dOut$lFDR_adjusted<0.1)) # selected threshold, 2% of total clusters DC

# What's the convergence between adjusted and unadjusted FDRs?
table(fdr_perm_dOut$lFDR_adjusted<0.1, fdr_perm_dOut$lFDR_individual<0.1)
# Perfect convergence, adjusted removes 23 contrasts

## Pseudo-volcano plots
fdr_perm_dOut$significantID<-fdr_perm_dOut$.id
fdr_perm_dOut$significantID[which(fdr_perm_dOut$lFDR_individual>0.1&fdr_perm_dOut$lFDR_adjusted>0.1)]<-NA
# fdr_perm_dOut$significantID<-droplevels(fdr_perm_dOut$significantID)
multiplot(cols = 2,
          ggplot(fdr_perm_dOut, aes(y=lFDR_individual, x=Estimate, col=significantID))+geom_point()+theme_bw()+scale_y_log10(limits=c(1e-5,1))+geom_vline(aes(xintercept=0))+geom_hline(aes(yintercept=0.1))+facet_wrap(~coefOutFactornames, nrow = 1),
          ggplot(fdr_perm_dOut, aes(y=lFDR_adjusted, x=Estimate, col=significantID))+geom_point()+theme_bw()+scale_y_log10(limits=c(1e-5,1))+geom_vline(aes(xintercept=0))+geom_hline(aes(yintercept=0.1))+facet_wrap(~coefOutFactornames, nrow = 1)
)

# How many clusters show sex-biased changes in constraint in how many stages?
table(table(fdr_perm_dOut[which(fdr_perm_dOut$lFDR_individual<0.1),]$.id))
table(table(fdr_perm_dOut[which(fdr_perm_dOut$lFDR_adjusted<0.1),]$.id))
# None

# Tabulate clusters with significant sex-biased changes in constraint, using adjusted lFDR
sign_dOut<-fdr_perm_dOut[which(fdr_perm_dOut$lFDR_adjusted<0.1),]
# Tabulate numer of clusters per stage and direction
table(sign_dOut$coefOutFactornames,sign(sign_dOut$Estimate))
# Greatest number of clusters in pupae and adults, absent from late embryo only
# Male-specific constraint (positive coefficients) predominates in all stages but early embryo

# Save names of clusters with significant changes in constraint, including stage and direction annotation
sign_dOut2<-sign_dOut[,c("coefOutFactornames",".id","Estimate","lFDR_adjusted")]
sign_dOut2$Sex<-as.factor(ifelse(sign_dOut2$Estimate>0,"Male","Female"))
names(sign_dOut2)<-c("Stage","clusterID","Estimate","lFDR_adjusted","Sex")
write.csv(sign_dOut2, file=file.path(newdir, "dOut_signClusters.csv"))

## Merge sign dWithin and dOut tests in single table
dOut <- data.frame( 
  clusterID = sign_dOut2$clusterID,
  DiffConstrained = paste(sign_dOut2$Sex,sign_dOut2$Stage, sep = "_")
)
dWithin <- data.frame( 
  clusterID = sign_dWithin2$clusterID,
  DiffIntegrated = paste(sign_dWithin2$Sex, sign_dWithin2$Stage, sep = "_")
)
signcontrasts <- merge(dWithin, dOut, all=T)
write.csv(signcontrasts, file = file.path(newdir, "DIDC_signclusters.csv"))
