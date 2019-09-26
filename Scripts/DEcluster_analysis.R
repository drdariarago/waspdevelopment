# DEcluster analyses
# Initialize script
date()
rm(list=ls())
library(ggplot2)
library(reshape2)
library(plyr)
library(stringr)
library(fdrtool)
source(file = "./Scripts/multiplot.R")

# Load dataset
load(file = "./Output/WGCNA_clustering_biweight/WGCNA_clustering_biweight.RData")

# Copy sample by cluster matrix
ME_expr <- MElist2$newMEs
# Add sample names
row.names(ME_expr) <- row.names(t_nasdevgeneexon)
# Convert to data.frame
ME_expr <- as.data.frame(ME_expr)
# Store sample stage
ME_expr$stage <- factor(as.character(str_extract(row.names(ME_expr), pattern = "^[^_]*")) ,levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
# Store sample sex
ME_expr$sex <- as.factor(ifelse(grepl("female", row.names(ME_expr)), "Female", "Male"))
# Store sample IDs
ME_expr$sampleID <- row.names(ME_expr)
# Remove other objects
rm(list=setdiff(ls(),c("ME_expr","multiplot")))

# initialize output path
newdir<-file.path(getwd(), "Output/DEcluster_analysis")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/DEcluster_analysis")
dir.create(graphdir)

# split into list of data.frames (one for each cluster)
ME_expr_list <- melt(data = ME_expr, id.vars = c("stage","sex","sampleID"), variable.name = "ClusterID", value.name = "E")
ME_expr_list <- dlply(.data = ME_expr_list, .variables = .(ClusterID))

### Analyze as GLMs
ME_LM_list <- lapply(X = ME_expr_list, FUN = lm, formula = E ~ 0+stage/sex, contrasts = list(stage="contr.Sum",sex="contr.treatment"))

# Perform GLM diagnostics
pdf(file = file.path(graphdir, "QCplots"))
lapply(X = ME_LM_list, FUN = plot)
dev.off()

## Save LMs to file

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

# Check for clusters above threshold and annotate sign
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

## Outdated code
# # Extract p-values (version w only 1 value per factor)
# pvals<-lapply(X = ME_LM_list, FUN = function(model){
#   pval <- anova(model)[5]
#   pval[1:3,]
# })
# pvals <- ldply(pvals)
# pvalmatrix <- as.matrix(pvals[2:4])
# row.names(pvalmatrix) <- pvals$.id
# names(pvalmatrix) <- c("clusterID","Stage","Sex","Sex:Stage")
# Apply multiple hypothesis correction