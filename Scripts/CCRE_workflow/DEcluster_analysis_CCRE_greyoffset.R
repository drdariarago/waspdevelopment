### DE clusters in CCRE network
# Using grey eigenexon as offset in glms
# Initialize script
date()
rm(list=ls())
library(ggplot2)
library(reshape2)
library(plyr)
library(stringr)
library(fdrtool)
source(file = "./Scripts/multiplot.R")

z <- function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}

# initialize output path
newdir<-file.path(getwd(), "Output/DEcluster_analysis_CCRE_greynormalized")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/DEcluster_analysis_CCRE_greynormalized")
dir.create(graphdir)

# Load datasets
load(file = "./Output/WGCNA_CCRE_clustering_biweight/CCRE_MEList.RData")
CCRE_nasdevgeneexon <- read.csv(file = "./Output/CCREs/CCRE_nasdevgeneexon.csv")
ME_expr <- MEList$eigengenes
row.names(ME_expr) <- colnames(CCRE_nasdevgeneexon)[-1]

# Store grey ME
MEgrey <- ME_expr[,"MEgrey"]
# ME_expr <- ME_expr[,-which(colnames(ME_expr)=="MEgrey")]
# # subtract and re-normalize from other clusters
# ME_expr <- apply(ME_expr, 2, function(x){x-MEgrey})
# ME_expr <- apply(ME_expr, 2, z)

# Convert to data.frame
ME_expr <- as.data.frame(ME_expr)
# Store sample stage
ME_expr$stage <- factor(as.character(str_extract(row.names(ME_expr), pattern = "^[^_]*")) ,levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
# Store sample sex
ME_expr$sex <- as.factor(ifelse(grepl("female", row.names(ME_expr)), "Female", "Male"))
# Store sample IDs
ME_expr$sampleID <- row.names(ME_expr)

# split into list of data.frames (one for each cluster)
ME_expr_list <- melt(data = ME_expr, id.vars = c("stage","sex","sampleID"), variable.name = "ClusterID", value.name = "E")
ME_expr_list <- dlply(.data = ME_expr_list, .variables = .(ClusterID))

### Analyze as GLMs
ME_LM_list <- lapply(X = names(ME_expr_list), FUN = function(x) {lm},
                     data = ME_expr_list[[x]],
                     formula = E ~ 0+stage/sex, 
                     contrasts = list(stage="contr.Sum",sex="contr.treatment"), offset = MEgrey)

# Perform GLM diagnostics
pdf(file = file.path(graphdir, "QCplots.pdf"))
lapply(X = ME_LM_list, FUN = plot)
dev.off()

## Save LMs to file
save(ME_LM_list, file = file.path(newdir, "ME_LMs.RData"))

## Extract coefficients
coefitall <- lapply(ME_LM_list, function(model){
  model$coefficients
})
coefmatrix <- as.matrix(ldply(.data = coefitall)[,-1])
row.names(coefmatrix) <- names(ME_LM_list)

## Extract p-values
pvals <- lapply(X = ME_LM_list, FUN = function(x){
  summary(x)[[4]][,4]
})
pvalmatrix <- as.matrix(ldply(.data = pvals)[,-1])
row.names(pvalmatrix) <- names(ME_LM_list)
## Apply fdr correction
FDRfitall<-fdrtool(as.vector(pvalmatrix), statistic="pvalue", verbose=T)
Fdrfitall<-matrix(FDRfitall$qval, 
                  ncol=ncol(pvalmatrix), nrow=nrow(pvalmatrix), dimnames=list(row.names(pvalmatrix), colnames(pvalmatrix))) # add global (tail area based)
fdrfitall<-matrix(FDRfitall$lfdr, 
                  ncol=ncol(pvalmatrix), nrow=nrow(pvalmatrix), dimnames=list(row.names(pvalmatrix), colnames(pvalmatrix))) # local fdr (density based)

# Check for clusters above threshold and annotate sign (sex only)
threshold<-5*10^-2
sexbiased_clusters_short <- apply(fdrfitall[,6:10], MARGIN = c(1,2), function(x){x<threshold})*apply(coefmatrix[,6:10], MARGIN = c(1,2), sign)
sexbiased_clusters_short <- apply(X = sexbiased_clusters_short, MARGIN = c(1,2), function(x){
  if (x==0) {"."} else if (x>0) {"m"} else if (x<0) {"f"}
})
sexbiased_clusters_short <- apply(X = sexbiased_clusters_short, 1, paste, collapse="")
sexbiased_clusters_short <- data.frame(names(ME_LM_list), sexbias = sexbiased_clusters_short)
sum(grepl("m|f", sexbiased_clusters_short$sexbias))
# Save as csv
names(sexbiased_clusters_short) <- c("clusterID", "cluster_devsexbias")
sexbiased_clusters_short$clusterID <- str_extract(string = sexbiased_clusters_short$clusterID, pattern = "[^ME].*")
write.csv(x = sexbiased_clusters_short, file = file.path(newdir, "sexbiased_clusters_short_5E-2.csv"))

# Check for clusters above threshold and annotate sign (stage only)
stagebiased_clusters_short <- apply(fdrfitall[,1:5], MARGIN = c(1,2), function(x){x<threshold})*apply(coefmatrix[,1:5], MARGIN = c(1,2), sign)
stagebiased_clusters_short <- apply(X = stagebiased_clusters_short, MARGIN = c(1,2), function(x){
  if (x==0) {"."} else if (x>0) {"+"} else if (x<0) {"-"}
})
stagebiased_clusters_short <- apply(X = stagebiased_clusters_short, 1, paste, collapse="")
stagebiased_clusters_short <- data.frame(names(ME_LM_list), stagebias = stagebiased_clusters_short)
sum(grepl("+|-", stagebiased_clusters_short$stagebias))
# Save as csv
names(stagebiased_clusters_short) <- c("clusterID", "cluster_stagebiased")
stagebiased_clusters_short$clusterID <- str_extract(string = stagebiased_clusters_short$clusterID, pattern = "[^ME].*")
write.csv(x = stagebiased_clusters_short, file = file.path(newdir, "stagebiased_clusters_short_5E-2.csv"))

# Merge stage and sexbias and save as csv
biased_clusters <- merge(sexbiased_clusters_short, stagebiased_clusters_short, by = "clusterID")
write.csv(biased_clusters, file = file.path(newdir, "biased_clusters.csv"))

# Plot individual eigengenes
# ggplot(data = ME_expr_list[[49]], aes(x=stage, y=E, col=sex))+geom_point()+theme_bw()

# Plot all clusters with modelled lms
ME_E_data <- melt(ME_expr, value.name = "E", variable.name = "moduleID")
pdf(file = file.path(graphdir,"cluster_lm.pdf"))
ggplot(data = ME_E_data, aes(x=stage, y=E, col=sex, group = sex))+geom_point()+geom_smooth(method="lm")+theme_bw()+facet_wrap(~moduleID)
dev.off()

# cp with bias without greynormalization
biased_nonormal <- read.csv(file = "./Output/DEcluster_analysis_CCRE/biased_clusters.csv")[,-1]
merge(biased_nonormal, biased_clusters, by = "clusterID", suffixes = c("_base", "_normal"), all.x = T)
grey_correction <- data.frame(
  row.names = biased_clusters[,1],
  devsexbias_base = as.character(biased_nonormal[-52,2]),
  devsexbias_corr = as.character(biased_clusters[,2]),
  stagebias_base = as.character(biased_nonormal[-52,3]),
  stagebias_corr = as.character(biased_clusters[,3]))
write.csv(x = grey_correction, file = file.path(newdir, "grey_correction.csv"))
