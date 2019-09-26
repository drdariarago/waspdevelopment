# Sex-specific analyses via GLMs
# clean space and load packages
date()
rm(list=ls())
library(plyr)
library(stringr)
library(WGCNA)
library(ggplot2)
library(fdrtool)
allowWGCNAThreads()
# initialize output path
newdir<-file.path(getwd(), "Output/stage_modularity_GLMs")
dir.create(newdir)

# import permutated cluster values
modulardensities <- read.csv("./Output/stage_modularity_biweight/modulardensities.csv")

## Fit GLMs for each module, validating the effect of stage and sex interaction on network parameters
# reorder stages
modulardensities$stage<-factor(modulardensities$stage, c("emb10","emb18", "lar51", "pupyel", "adult"))

# For changes in integration coefficient
dWithinGLMs<-dlply(.data = modulardensities, 
                   .variables = .(clusterID),
                   .fun = glm, 
                   formula = dWithin~stage1+stage2+stage3+stage4+stage5, 
                   family = gaussian(link="logit"),
                   .progress = "text")
# Store diagnostic plots
pdf(file = file.path(newdir, "dWithinGLM_diagnostics.pdf"))
par(mfrow=c(2,2))
lapply(X = dWithinGLMs, FUN = plot)
dev.off()

# For changes in pleiotropy coefficient
dOutGLMs<-dlply(.data = modulardensities, 
                .variables = .(clusterID),
                .fun = glm, 
                formula = dOut~stage1+stage2+stage3+stage4+stage5, 
                family = gaussian(link="logit"),
                .progress = "text")
# Store diagnostic plots
pdf(file = file.path(newdir, "dOutGLM_diagnostics.pdf"))
par(mfrow=c(2,2))
lapply(X = dOutGLMs, FUN = plot)
dev.off()


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