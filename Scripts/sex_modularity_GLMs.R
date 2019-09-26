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
library(boot)
library(reshape2)
# initialize output path
newdir<-file.path(getwd(), "Output/sex_modularity_GLMs")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/sex_modularity_GLMs")
dir.create(graphdir)

# import permutated cluster values
modulardensities <- read.csv("./Output/sex_modularity_biweight_2_01/sex_modulewise_modulardensities.csv")

## Fit GLMs for each module, validating the effect of stage and sex interaction on network parameters
# reorder stages

modulardensities$stage<-factor(modulardensities$stage, c("emb10","emb18", "lar51", "pupyel", "adult"))

# Are there any densities greater than 1 or lesser than 0?
summary(apply(X = modulardensities[,c("dWithin", "dOut", "dMain")], MARGIN = 1, FUN = function(x){any(x>1)}))
summary(apply(X = modulardensities[,c("dWithin", "dOut", "dMain")], MARGIN = 1, FUN = function(x){any(x<0)}))

# # # Which ones?
# modulardensities[which(apply(X = modulardensities[,c("dWithin", "dOut", "dMain")], MARGIN = 1, FUN = function(x){any(x>1)})),]
# modulardensities[which(apply(X = modulardensities[,c("dWithin", "dOut", "dMain")], MARGIN = 1, FUN = function(x){any(x<0)})),]

# is there any variance/mean correlation?
# ggplot(data=ddply(modulardensities, .(stage, sexbias, clusterID), summarize, dWithinM=mean(dWithin), dWithinS=sd(dWithin)), aes(x=dWithinM,y=dWithinS, col=clusterID))+geom_point()+geom_smooth()
with(ddply(modulardensities, .(stage, sexbias, clusterID), summarize, dWithinM=mean(dWithin), dWithinS=sd(dWithin)), cor(dWithinM, dWithinS, use = "pair"))
with(ddply(modulardensities, .(stage, sexbias, clusterID), summarize, dOutM=mean(dOut), dOutS=sd(dOut)), cor(dOutM, dOutS, use = "pair"))
# no in dWithin, yes in dOut

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

#### Testing using bootstrapping

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

# # Create index of permutations (one permutation per row) (recycling from dWithin)
# indexes<-replicate(
#   n = nperm,
#   sample(x = 1:table(modulardensities$clusterID)[1], size = table(modulardensities$clusterID)[1], replace = F)
# )

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

#### OUTDATED code snippets

# #### Testing using bootstrapping, dWithin version with replacement
# 
# # Enable error logging
# # options(error=recover)
# 
# ## Create bootstrapping function
# bf_dWithin<-function(modulardensities,index){
#   newdata<-modulardensities
#   newdata[,"sexbias"]<-newdata[index,"sexbias"] # randomly reassign sexbias values
#   #   newdata<-data.frame(modulardensities[,-4],sexbias=modulardensities[index,4]) # alternative version for swapping sexbias indexes
#   model<-glm(dWithin~stage-1+stage:sexbias+dMain+dOut, data = newdata, family=Gamma(link="logit"))
#   #   coef(model) # save coefficient estimates
#   summary(model)$coefficients[c(8:12),"Pr(>|t|)"] # save p-value estimates of stage:sex interactions
# }
# 
# # Apply bootstrapping function to every module, generate list of boot items
# 
# boot_dWithin<-dlply(.data = modulardensities, .variables = .(clusterID), .fun = boot, statistic=bf_dWithin, sim="ordinary", R=nperm)
# 
# # Calculate local FDR of p-values compared to null distribution of reshuffled labels
# 
# pdf(file = file.path(graphdir, "dWithin_boot_cluster_fdr.pdf"))
# dWithin_boot_q.val_1<-ldply(lapply(X = boot_dWithin, FUN = function(cluster){
#   sapply(X = 1:5, function(index){
#     fdrtool((c(cluster$t[,index],cluster$t0[index])), statistic = "pvalue", plot = T)$lfdr[nperm+1]
#   })
# }))
# dev.off()
# colnames(dWithin_boot_q.val_1)<-c("clusterID", names(boot_dWithin[[1]]$t0))
# 
# # Adjust lFDR values for multiple hypothesis testing
# pdf(file = file.path(graphdir, "dWithin_boot_Fdr.pdf"))
# dWithin_boot_q.val_2<-matrix(fdrtool(c(unlist(dWithin_boot_q.val_1[,-1])), statistic = "pvalue", plot = T)$qval,
#                              ncol = 5, 
#                              dimnames = list(
#                                c(dWithin_boot_q.val_1$clusterID),
#                                c(colnames(dWithin_boot_q.val_1)[-1]))
# )
# dev.off()
# 
# # Attach to factor estimates
# dWithin_boot_q.val_1<-melt(dWithin_boot_q.val_1)
# names(dWithin_boot_q.val_1)<-c(".id","coefWithinFactornames", "lFDR_individual")
# dWithin_boot_q.val_2<-melt(dWithin_boot_q.val_2)
# names(dWithin_boot_q.val_2)<-c(".id","coefWithinFactornames", "lFDR_adjusted")
# # factor_q.val_2$coefOutFactornames<-str_extract(factor_q.val_2$coefOutFactornames, pattern = "^[^:]*")
# fdr_boot_dWithin<-merge(coefWithin, dWithin_boot_q.val_1)
# fdr_boot_dWithin<-merge(fdr_boot_dWithin, dWithin_boot_q.val_2)
# 
# # And save as csv
# write.csv(fdr_boot_dWithin, file=file.path(newdir, "fdr_boot_dWithin.csv"))
# 
# #### Testing using bootstrapping, dOut version with replacement (not working, wrong indexes in function)
# 
# ## Create bootstrapping function
# bf_dOut<-function(modulardensities,index){
#   newdata<-modulardensities
#   newdata[,"sexbias"]<-newdata[index,"sexbias"] # randomly reassign sexbias values
#   #   newdata<-data.frame(modulardensities[,-4],sexbias=modulardensities[index,4]) # alternative version for swapping sexbias indexes
#   model<-glm(dOut~stage-1+stage:sexbias+dMain+dWithin, data = newdata, family=Gamma(link="logit"))
#   #   coef(model) # save coefficient estimates
#   summary(model)$coefficients[c(8:12),"Pr(>|t|)"] # save p-value estimates of stage:sex interactions
# }
# 
# # Apply bootstrapping function to every module, generate list of boot items, check to avoid resampling
# 
# boot_dOut<-dlply(.data = modulardensities, .variables = .(clusterID), .fun = boot, statistic=bf_dOut, sim="ordinary", R=nperm)
# 
# # Calculate local FDR of p-values compared to null distribution of reshuffled labels
# 
# pdf(file = file.path(graphdir, "dOut_boot_cluster_fdr.pdf"))
# dOut_boot_q.val_1<-ldply(lapply(X = boot_dOut, FUN = function(cluster){
#   sapply(X = 1:5, function(index){
#     fdrtool((c(cluster$t[,index],cluster$t0[index])), statistic = "pvalue", plot = T)$lfdr[nperm+1]
#   })
# }))
# colnames(dOut_boot_q.val_1)<-c("clusterID", names(boot_dOut[[1]]$t0))
# dev.off()
# 
# # Adjust lFDR values for multiple hypothesis testing
# pdf(file = file.path(graphdir, "dOut_boot_fdr.pdf"))
# dOut_boot_q.val_2<-matrix(fdrtool(c(unlist(dOut_boot_q.val_1[,-1])), statistic = "pvalue", plot = T)$qval,
#                           ncol = 5, 
#                           dimnames = list(
#                             c(dOut_boot_q.val_1$clusterID),
#                             c(colnames(dOut_boot_q.val_1)[-1]))
# )
# dev.off()
# 
# # Attach to factor estimates
# dOut_boot_q.val_1<-melt(dOut_boot_q.val_1)
# names(dOut_boot_q.val_1)<-c(".id","coefOutFactornames", "lFDR_individual")
# dOut_boot_q.val_2<-melt(dOut_boot_q.val_2)
# names(dOut_boot_q.val_2)<-c(".id","coefOutFactornames", "lFDR_adjusted")
# # factor_q.val_2$coefOutFactornames<-str_extract(factor_q.val_2$coefOutFactornames, pattern = "^[^:]*")
# fdr_boot_dOut<-merge(coefOut, dOut_boot_q.val_1)
# fdr_boot_dOut<-merge(fdr_boot_dOut, dOut_boot_q.val_2)
# 
# # And save as csv
# write.csv(fdr_boot_dOut, file=file.path(newdir, "fdr_boot_dOut.csv"))


# # plot variation in network parameters for each cluster
# ggplot(modularconnectivities, aes(y=kOut/kTotal, x=sexbias, col=clusterID))+geom_point(alpha="0.1")+geom_line()+facet_wrap(~clusterID)+scale_y_log10()
# ggplot(modularconnectivities, aes(y=kWithin/kTotal, x=sexbias, col=clusterID))+geom_point(alpha="0.1")+geom_line()+facet_wrap(~clusterID)+scale_y_log10()
# ggplot(modularconnectivities, aes(y=kWithin/kOut, x=sexbias, col=clusterID))+geom_point(alpha="0.1")+geom_line()+facet_wrap(~clusterID)
# 
# # summarize module data
# a<-ddply(modularconnectivities[,c(1:7,11,12)], .variables = .(clusterID, sexbias, stage), .fun = numcolwise(median))
# names(a)[grep("k", names(a))]<-paste0(names(a)[grep("k", names(a))], "_median")
# b<-ddply(modularconnectivities[,c(1:7,11,12)], .variables = .(clusterID, sexbias), .fun = numcolwise(mad))
# names(b)[grep("k", names(b))]<-paste0(names(b)[grep("k", names(b))], "_mad")
# c<-merge(a,b)
# ggplot(c, aes(x=sexbias, y=kWithin_median, col=clusterID))+facet_grid(.~stage)+geom_ribbon(aes(ymin=(kWithin_median)-(kWithin_mad), ymax=(kWithin_median/kTotal_median)+(kWithin_mad)))+geom_line()+theme_bw()
# 
# ggplot(c, aes(x=sexbias, y=kWithin_median, col=clusterID))+facet_grid(.~stage)+geom_line()+theme_bw()+geom_pointrange(aes(ymin=kWithin_median-kWithin_mad, ymax=kWithin_median+kWithin_mad))
# 

# ##
# 
# # BOOT: Generate table of factor by confidence interval (95%), check appropriate namings, apply to whole list of boot models
# 
# b<-t(sapply(X = c(8:12), function(x){
#   bootCI<-boot.ci(boot_dWithin, index = x, type = "bca")
#   bootCI$bca[,4:5]
# }))
# 
# # Check whether p-values fall within the random interval of p-values (95% confidence intervals)

## Compare null distribution to p-vaues using single sample t.test

# t-test of boot$t0 vs boot$t (sample vs distribution, logit transformed p-values to normalize)
# logit<-function(x){log10(x/(1-x))}
# 
# factor_p.val<-ldply(lapply(X = boot_dWithin, FUN = function(cluster){
#   sapply(X = 1:5, function(index){
#     wilcox.test(x = cluster$t[,index], mu = cluster$t0[index], alternative = "greater")$p.value
#   })
# }))
# colnames(factor_p.val)<-c("clusterID", names(boot_dWithin[[1]]$t0))

#  Annotate also mean/sd
# factor_p.val<-ldply(lapply(X = boot_dWithin, FUN = function(cluster){
#   sapply(X = 1:5, function(index){
#     cbind(
#       cluster$t0[index],
#       t.test(x = logit(cluster$t[,index]), mu = logit(cluster$t0[index]), alternative = "greater")$p.value,
#       mean(cluster$t[,index]),
#       sd(cluster$t[,index], na.rm = T)
#     )
#   })
# }))
# factor_p.val$parameter<-c("p","delta.p","mean.null.p","sd.null.p")


### Testing Bootstrapping functions (with replacement)
# Single cluster version for testing
# boot_dWithin<-boot(data=modulardensities[which(modulardensities$clusterID=="aliceblue"),], statistic=bf_dWithin, sim="ordinary", R=500)

# # Restricted to few clusters for testing on list
# boot_dWithin<-dlply(.data = modulardensities[grep("white|red|blue",modulardensities$clusterID),], .variables = .(clusterID), 
#                     .fun = boot, statistic=bf_dWithin, sim="ordinary", R=nperm)
