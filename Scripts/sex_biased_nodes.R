## Define sex-biased nodes
# Initialize script
date()
rm(list=setdiff(ls(),"ODB8_EukOGs_genes"))
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(limma)
library(fdrtool)
# initialize output path
newdir<-file.path(getwd(), "Output/sex_biased_nodes")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/sex_biased_nodes")
dir.create(graphdir)

# load datasets
eigenexon_evalues <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexon_evalues.csv")
sexspec<- read.csv("./Output/sex_specific_nodes/sexspec_allstages.csv")[,-1]

# Remove sex-specific nodes
eigenexon_evalues<-eigenexon_evalues[which(eigenexon_evalues$eigenexonID%in%sexspec[which(is.na(sexspec$spec==F)),"eigenexonID"]),]

# Convert to expression list
expr_list<-list(E=as.matrix(eigenexon_evalues[,-1]), genes=as.character(eigenexon_evalues[,1]), targets=names(eigenexon_evalues[,-1]))
# Converting to LIMMA EList format
expr_list<-new("EList", expr_list)

# Z transform separately transcription and splicing nodes, setting minimum to zero (rescales to a common distribution for easier comparisons)
z<-function(x){
#   y<-x-mean(x, na.rm = T)
  y<-x-(min(x, na.rm = T)*sign(min(x, na.rm = T)))
  y<-y/sd(y, na.rm = T)
}

expr_list$E[grep("fac", expr_list$genes),]<-z(expr_list$E[grep("fac", expr_list$genes),])
expr_list$E[grep("con", expr_list$genes),]<-z(expr_list$E[grep("con", expr_list$genes),])

# Create factors
stage<-factor(as.character(str_extract(expr_list$targets, pattern = "^[^_]*")) ,levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
stage_ordered<-ordered(stage, levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
sex<-as.factor(ifelse(grepl("female", expr_list$targets), "Female", "Male"))
# Create design matrix
design<-model.matrix(~0+stage/sex, contrasts.arg=list(stage="contr.Sum",sex="contr.treatment"))
# design<-model.matrix(~0+stage/sex, contrasts.arg=list(stage="contr.treatment",sex="contr.treatment"))
# check if data matrix is correct
as.data.frame(cbind(expr_list$targets, as.character(stage), as.character(sex), design))

# Fit linear models
glms<-lmFit(expr_list, design)
glms2<-eBayes(glms)

# Save exon models and EList
save(file=file.path(newdir, "LIMMA_sexbias"), list=c("glms2"))

# This version applies FDR correction to all p-values
FDRfitall<-fdrtool(as.vector(glms2$p.value), statistic="pvalue", verbose=T)
Fdrfitall<-matrix(FDRfitall$qval, 
                  ncol=ncol(glms2$p.value), nrow=nrow(glms2$p.value), dimnames=list(glms2$genes, colnames(glms2$p.value))) # add global (tail area based)
fdrfitall<-matrix(FDRfitall$lfdr, 
                  ncol=ncol(glms2$p.value), nrow=nrow(glms2$p.value), dimnames=list(glms2$genes, colnames(glms2$p.value))) # local fdr (density based)

# Save Fdr and fdr to base model
glms2$Fdr<-Fdrfitall
glms2$fdr<-fdrfitall

# and save as csv
write.csv(fdrfitall, file=file.path(newdir, "LIMMA_sexbias_fdr"))

# Summarize, number of nodes significantly sexbiased and number of expected false discoveries per fdr threshold
thresholds<-alply(c(5*10^-(1:5)), .margins = 1)
discoveries<-sapply(X = thresholds, FUN =  function(x){sum(c(glms2$fdr)<x)})
false_disoveries<-discoveries*unlist(thresholds)
fdr_table<-rbind(discoveries, false_disoveries)
fdr_table # Choosing lfdr 0.05, expected 7703 false discoveries (individual prob for each gene 0.05)
# Percent proportions of all contrasts
(fdr_table/length(expr_list$E))*100

# Set lfdr threshold
threshold<-5*10^-2
# Save coefficients and p-values
glms2_coefs<-data.frame(glms2$coefficients, row.names = glms2$genes)
names(glms2_coefs)<-colnames(fdrfitall)
glms2_fdr_coef<-merge(fdrfitall,glms2_coefs, by="row.names", suffixes = c("_fdr", "_coef"))
names(glms2_fdr_coef)[1]<-"node_ID"
names(glms2_fdr_coef)[-1]<-str_sub(names(glms2_fdr_coef)[-1], start = 6)
names(glms2_fdr_coef)[-1]<-str_replace(names(glms2_fdr_coef)[-1], pattern = ":", replacement = "_")
names(glms2_fdr_coef)[-1]<-str_replace(names(glms2_fdr_coef)[-1], pattern = "sex", replacement = "")

write.csv(glms2_fdr_coef, file=file.path(newdir, "glms2_fdr_coef.csv"))

# Save table of nodes with compressed sexbias
sexbiased_nodes_short <- apply(glms2_fdr_coef[,7:11], MARGIN = c(1,2), function(x){x<threshold})*apply(glms2_fdr_coef[,17:21], MARGIN = c(1,2), sign)
sexbiased_nodes_short <- apply(X = sexbiased_nodes_short, MARGIN = c(1,2), function(x){
  if (x==0) {"."} else if (x>0) {"m"} else if (x<0) {"f"}
})
sexbiased_nodes_short <- apply(X = sexbiased_nodes_short, 1, paste, collapse="")
sexbiased_nodes_short <- data.frame(glms2_fdr_coef$node_ID, sexbias = sexbiased_nodes_short)
sum(grepl("m|f", sexbiased_nodes_short$sexbias))
write.csv(x = sexbiased_nodes_short, file = file.path(newdir, "sexbiased_nodes_short_5E-2.csv"))

# Reshape dataset for summary statistics
melt_glms2_fc<-melt(glms2_fdr_coef)
melt_glms2_fc$stage<-factor(as.character(str_extract(melt_glms2_fc$variable, pattern = "^[^_]*")) ,levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
melt_glms2_fc$sex<-ifelse(grepl("Male", melt_glms2_fc$variable), "Sex", "Stage")
melt_glms2_fc$type<-factor(str_extract(string = melt_glms2_fc$variable, pattern = "[^_]*$"))
melt_glms2_fc$spl<-as.factor(ifelse(grepl("fac", melt_glms2_fc$node_ID), "Splicing", "Transcription"))
glms2_fc<-dcast(data = melt_glms2_fc, formula = node_ID+stage+sex+spl~type)


# # Remove conflict nodes and annotate sexbiased ones
# sxbnodes<-sxbnodes[which(!(row.names(sxbnodes)%in%row.names(sxbconflict))),]
# Save as csv
sxbnodes2<-glms2_fc[which(glms2_fc$sex=="Sex"),c("node_ID", "stage", "coef", "fdr")]
sxbnodes2$sex<-ifelse(sign(sxbnodes2$coef)>0,"Male","Female")
sxbnodes2$sex[which(sxbnodes2$fdr>=threshold)]<-NA
sxbnodes2$gene_ID<-str_extract(string = sxbnodes2$node_ID, pattern = "^[^_]*")
sxbnodes2<-sxbnodes2[,c("gene_ID","node_ID","sex")]
names(sxbnodes2)<-c("gene_ID", "eigenexonID", "sexbias")
sxbnodes2<-unique(sxbnodes2)
write.csv(sxbnodes2, file = file.path(newdir, "sexbiased_nodes.csv"))

# save list of male and female biased genes
write(as.character(sxbnodes2[which(sxbnodes2$sexbias=="Male"),]$gene_ID), file=file.path(newdir, "male_bias.txt"))
write(as.character(sxbnodes2[which(sxbnodes2$sexbias=="Female"),]$gene_ID), file=file.path(newdir, "female_bias.txt"))