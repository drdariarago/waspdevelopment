### Explore results from Motif Enrichment Tool
# Initialize script
date()
# rm(list=ls())
newdir<-file.path(getwd(), "Output/MET_overview")
graphicsdir<-file.path(getwd(), "Graphics/MET_overview")
dir.create(newdir)
dir.create(graphicsdir)
library(plyr)
library(ggplot2)
library(stringr)
library(fdrtool)

# Load module metadata
moduledata <- read.csv(file = "./Output/Results_compiler/moduledata_full.csv")

# Only need to run once after MET run
# Import data
MET_Results <- read.delim("./Input/MET_output/contrasts_networkgenes_500targets_005.txt", comment.char="#")
# Import clean set
MET_Results <- read.csv(file = "./Output/MET_overview/333clusters_allgenes_1000targets_005_slim.csv")


# # Import individual datasets for each stage:sex combination
# datalist = list.files(path = "./Input/MET_output/Compiled_dI_clusters/", pattern = "ale_")
# datalist = lapply(X = datalist, function(x) read.delim(file = file.path("./Input/MET_output/Compiled_dI_clusters/",x)))
# names(datalist) <- str_extract(string = list.files(path = "./Input/MET_output/Compiled_dI_clusters/", pattern = "ale_"), pattern = "^[^_]*")
# # Merge into single data.frame
# MET_Results <- ldply(.data = datalist, .id = "contrast")

# Remove extra title rows
MET_Results <- MET_Results[-which(MET_Results$UserSet=="UserSet"),c("UserSet","Normalization","RegRegion","Motif","User_ct","Intersection","Significance")]
# Export slimmed dataset
write.csv(x = MET_Results, file = file.path(newdir, "contrasts_networkgenes_500targets_005.csv"))
# Re-import slimmed dataset
MET_Results<-read.csv(file = "./Output/MET_overview/contrasts_networkgenes_500targets_005.csv")[,-1]
summary(MET_Results)

# # # Subselect modules with interesting sexbias dynamics
# sex_modules <- vector()
# sex_modules <- which(apply(moduledata[,grep("enriched", names(moduledata))], 1, function(x){any(x=="Enriched")}))
# sex_modules <- union(sex_modules, which(is.na(moduledata$dWithinBias)==F))
# sex_modules <- moduledata[sex_modules, "clusterID"]
# length(sex_modules)
# dim(MET_Results)
# MET_Results <- MET_Results[which(MET_Results$UserSet%in%sex_modules),]
# dim(MET_Results)

# Distribution of significance across all motifs in different regions with different normalization procedures
ggplot(data = MET_Results, aes(x=Significance, col=RegRegion, lty=Normalization))+geom_density()+theme_bw()+scale_x_log10()
ggplot(data = MET_Results, aes(y=Significance, x=RegRegion, col=RegRegion, lty=Normalization))+geom_boxplot(notch = T, notchwidth = T)+theme_bw()+scale_y_log10()

# Distribution of significance across significant motifs in different regions with different normalization procedures
pdf(file = file.path(graphicsdir, "boxplot_sign_regions.pdf"), width = 10)
ggplot(data = MET_Results, aes(y=Significance, x=RegRegion, col=RegRegion, lty=Normalization))+geom_boxplot(notch = T, notchwidth = T)+theme_bw()+scale_y_log10(limits=c(1E-16,0.05))+ggtitle(label = "Significance of significant motifs (p<0.05) \nin different regions and normalization procedures")+scale_x_discrete(labels=NULL, name="Regulatory Region Type")+scale_color_discrete(name="Regulatory Region Type")
dev.off()

# Subset to region of interest
MET_Results <- MET_Results[which(MET_Results$Normalization=="normalized by GC content"&MET_Results$RegRegion=="1kb upstream of TSS"),]

## Version from raw p-values
# Check number of hypotheses
nrow(MET_Results) # 2002, 2 exp false positives with p-value 0.001
table(MET_Results$Significance<0.001) # 3 positives

# Convert to fdrs
MET_fdr<-fdrtool(x = MET_Results$Significance, statistic = "pvalue")
min(MET_fdr$qval)
min(MET_fdr$lfdr)

MET_Results$qval <- MET_fdr$qval
MET_Results$lfdr <- MET_fdr$lfdr

threshold<-0.5
table(MET_Results$qval<threshold)
table(MET_Results$qval<threshold)*threshold 

## Create results table. For each cluster, report scoring motifs with region and probability. Split by normalization procedure.
MET_summary<-MET_Results[which(MET_Results$qval<threshold),c("Normalization","contrast","Motif","lfdr","qval")]
MET_summary<- MET_summary[order(MET_summary$Normalization, MET_summary$UserSet, MET_summary$RegRegion,MET_summary$Significance),]
# MET_summary<- MET_summary[order(MET_summary$Normalization, MET_summary$UserSet, MET_summary$Motif ,MET_summary$Significance),]
# And write to CSV
write.csv(MET_summary, file=file.path(newdir, "MET_summary.csv"))


############# Outdated code snippets

# # Tabulate expected false discoveries per threshold, per region (using standard normalization)
# sapply(X = levels(MET_Results$RegRegion), function(x){
#   sapply(X = c(0.1,0.05,0.01,0.05), FUN = function(y)table(MET_Results$Significance[which(MET_Results$RegRegion==x&MET_Results$Normalization=="standard normalization")]<y))
# })
# 
# sapply(X = levels(MET_Results$RegRegion), function(x){
#   sapply(X = c(0.1,0.05,0.01,0.05), FUN = function(y)table(MET_Results$Significance[which(MET_Results$RegRegion==x&MET_Results$Normalization=="standard normalization")]<y)*y)
# })
# 
# # Tabulate expected false discoveries per threshold, per region (using GC normalization)
# sapply(X = levels(MET_Results$RegRegion), function(x){
#   sapply(X = c(0.1,0.05,0.01,0.05), FUN = function(y)table(MET_Results$Significance[which(MET_Results$RegRegion==x&MET_Results$Normalization=="normalized by GC content")]<y))
# })
# 
# sapply(X = levels(MET_Results$RegRegion), function(x){
#   sapply(X = c(0.1,0.05,0.01,0.05), FUN = function(y)table(MET_Results$Significance[which(MET_Results$RegRegion==x&MET_Results$Normalization=="normalized by GC content")]<y)*y)
# })
# 
# # 0.01 seems an appropriate threshold in most cases (<1 FD), maybe too stringent

# ## Which motifs are specific of a given region rather than shared (in case of doubt 1KB up seems to capture most)
# # Split dataset by Userset(cluster), Motif and Normalization, Then check if duplicates are present amongst significant motifs
# RegSpec<-ddply(.data = MET_Results[which(MET_Results$Significance<0.05),], .variables = .(Normalization, contrast, Motif), summarize, RegSpec=length(Significance))
# prop.table(table(RegSpec$RegSpec==1, RegSpec$Normalization), margin = 1)
