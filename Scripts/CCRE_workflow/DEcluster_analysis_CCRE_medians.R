### DE clusters in CCRE network based on mean expression values
# Initialize script
date()
rm(list=ls())
library(ggplot2)
library(reshape2)
library(plyr)
library(stringr)
library(fdrtool)
source(file = "./Scripts/multiplot.R")

# initialize output path
newdir<-file.path(getwd(), "Output/DEcluster_analysis_CCRE")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/DEcluster_analysis_CCRE")
dir.create(graphdir)

# Load datasets
load(file = "./Output/WGCNA_CCRE_clustering_biweight/CCRE_MEList.RData")
CCRE_nasdevgeneexon <- read.csv(file = "./Output/CCREs/CCRE_nasdevgeneexon.csv")
ME_expr <- MEList$eigengenes
row.names(ME_expr) <- colnames(CCRE_nasdevgeneexon)[-1]

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

#### Compute based on average normalized expression instead of module eigengene (sexbias could be in components other than the first one) 
# low overall proportion of variance explained by first PCs
dim(MEList$varExplained)
varExplained <- as.data.frame(t(MEList$varExplained))
varExplained$clusterID <- colnames(MEList$eigengenes)
varExplained <- melt(varExplained)
varExplained$variable <- as.integer(varExplained$variable)
ggplot(data = varExplained, aes(y=value, x=variable, col=clusterID))+geom_line()+scale_x_discrete()+scale_y_log10()+theme_bw()+xlab(label = "Number of Principal Component") + ylab(label = "Proportion of Variance Explained")+ggtitle(label = "Proportion of Variance Explained by each Eigengene vs other PCs")
## Compute based on gene averages
# Convert to data.frame
MA_expr <- as.data.frame(MEList$averageExpr)
row.names(MA_expr) <- colnames(CCRE_nasdevgeneexon)[-1]
# Store sample stage
MA_expr$stage <- factor(as.character(str_extract(row.names(MA_expr), pattern = "^[^_]*")) ,levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
# Store sample sex
MA_expr$sex <- as.factor(ifelse(grepl("female", row.names(MA_expr)), "Female", "Male"))
# Store sample IDs
MA_expr$sampleID <- row.names(MA_expr)

# split into list of data.frames (one for each cluster)
MA_expr_list <- melt(data = MA_expr, id.vars = c("stage","sex","sampleID"), variable.name = "ClusterID", value.name = "E")
MA_expr_list <- dlply(.data = MA_expr_list, .variables = .(ClusterID))

### Analyze as GLMs
MA_LM_list <- lapply(X = MA_expr_list, FUN = lm, formula = E ~ 0+stage/sex, contrasts = list(stage="contr.Sum",sex="contr.treatment"))

# Perform GLM diagnostics
pdf(file = file.path(graphdir, "QCplots_MA"))
lapply(X = MA_LM_list, FUN = plot)
dev.off()

## Save LMs to file
save(MA_LM_list, file = file.path(newdir, "MA_LMs.RData"))

## Extract coefficients
coefitall <- lapply(MA_LM_list, function(model){
  model$coefficients
})
coefmatrix <- as.matrix(ldply(.data = coefitall)[,-1])
row.names(coefmatrix) <- names(MA_LM_list)

## Extract p-values
pvals <- lapply(X = MA_LM_list, FUN = function(x){
  summary(x)[[4]][,4]
})
pvalmatrix <- as.matrix(ldply(.data = pvals)[,-1])
row.names(pvalmatrix) <- names(MA_LM_list)
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
sexbiased_clusters_short$clusterID <- str_extract(string = sexbiased_clusters_short$clusterID, pattern = "[^AE].*")
write.csv(x = sexbiased_clusters_short, file = file.path(newdir, "MA_sexbiased_clusters_short_5E-2.csv"))

# Only difference cluster yellow, which has male bias only in module eigengene testing adults