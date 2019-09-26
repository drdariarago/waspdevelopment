### DE clusters in CCRE network
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

### Analyze as GLMs
ME_LM_list <- lapply(X = ME_expr_list, FUN = lm, formula = E ~ 0+stage/sex, contrasts = list(stage="contr.Sum",sex="contr.treatment"))

# Perform GLM diagnostics
pdf(file = file.path(graphdir, "QCplots"))
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

## Save coefficients and pvalues for each cluster
m_coefmatrix <- melt(coefmatrix)
names(m_coefmatrix) <- c("clusterID", "contrast", "coeff")
m_fdrfitall <- melt(fdrfitall)
names(m_fdrfitall) <- c("clusterID", "contrast", "fdr")
clusterGLMsummaries <- merge(m_coefmatrix, m_fdrfitall)
clusterGLMsummaries$stage <- substring(text = str_extract(string = clusterGLMsummaries$contrast, pattern = "^[^:]*"), first = 6)
clusterGLMsummaries$stage <- factor(x = clusterGLMsummaries$stage, levels = c("emb10","emb18","lar51","pupyel","adult"))
clusterGLMsummaries$contrast <- factor(ifelse(grepl(pattern = "Male", x = clusterGLMsummaries$contrast), "Sex", "Stage")) 
clusterGLMsummaries <- clusterGLMsummaries[order(clusterGLMsummaries$clusterID, clusterGLMsummaries$stage, clusterGLMsummaries$contrast),c("clusterID","stage","contrast","fdr","coeff")]
write.csv(x = clusterGLMsummaries, file = file.path(newdir, "clusterGLMsummaries.csv"))

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