### Analyze gonad data to determine which genes are gonad specific and which ones are male vs female gonad specific rather than shared

library(stringr)
library(magrittr)
# library(reshape2)
library(plyr)
library(limma)
library(fdrtool)

## Create output directories
newdir<-file.path(getwd(), "Output/gonads_DE")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/gonads_DE")
dir.create(graphdir)


## Load datasets
gonad_expr <- read.csv(file = "./Output/gonads_cluster/gonad_expr_CCRE.csv", row.names = 1)

## Female soma vs gonads

# Convert to expression list
expr_list<-list(
  E = as.matrix(gonad_expr), 
  genes = as.character(gonad_expr),
  targets=names(gonad_expr))

# Converting to LIMMA EList format
expr_list<-new("EList", expr_list)

# Create factors
soma <- factor(x = ifelse(grepl("testes|ovaries", colnames(gonad_expr)), "Gonads", "Soma"), 
  levels = c("Soma", "Gonads")) 
sex <- factor(x = ifelse(grepl("adult|ovaries", colnames(gonad_expr)), "Female", "Male"), 
  levels = c("Female", "Male")) 

# Create design matrix
# Check male and female baseline separately
# Contrast all gonads vs all soma
# Contrast male gonads vs female gonads

design <- model.matrix(
  ~ 0 + sex * soma, 
  contrasts.arg=list(
    sex = "contr.sum",
    soma = "contr.treatment")
)

data.frame(sex, soma, design)

# Fit linear models
glms<-lmFit(expr_list, design)
glms2<-eBayes(glms)
# Save exon models and EList
save(file=file.path(newdir, "LIMMA_sexbias_gonads"), list=c("glms2"))

## Annotate contrast names in human-readable format

contrast_names <- c("Female_Upregulated", "Male_Upregulated", "Gonads_vs_Soma", "Ovaries_vs_Testes")

# Apply fdr correction
FDRfitall<-fdrtool(as.vector(glms2$p.value), statistic="pvalue", verbose=T)

Fdrfitall<-matrix(FDRfitall$qval, 
  ncol=ncol(glms2$p.value), nrow=nrow(glms2$p.value), 
  dimnames=list(row.names(glms2$p.value), contrast_names)) 

fdrfitall<-matrix(FDRfitall$lfdr, 
  ncol=ncol(glms2$p.value), nrow=nrow(glms2$p.value), 
  dimnames=list(row.names(glms2$p.value), contrast_names)) 

# Save Fdr and fdr to base model
glms2$Fdr<-Fdrfitall
glms2$fdr<-fdrfitall

# and save as csv
write.csv(fdrfitall, file=file.path(newdir, "LIMMA_sexbias_gonads_fdr"))


## Summarize discoveries
thresholds<-alply(c(5*10^-(1:5)), .margins = 1)
names(thresholds) <- thresholds
discoveries<-sapply(X = thresholds, FUN =  function(x){sum(c(glms2$fdr)<x)})
false_disoveries<-discoveries*unlist(thresholds)
fdr_table<-rbind(discoveries, false_disoveries)
fdr_table # Choosing lfdr 0.05, expected 7703 false discoveries (individual prob for each gene 0.05)
# Percent proportions of all contrasts
(fdr_table/length(expr_list$E))*100

# Set lfdr threshold
threshold<-5*10^-2


## Save table of nodes, gonad vs soma bias and ovaries vs testes bias
gonads_summary_signif <- glms2$fdr[,c(3,4)] %>%
  apply(., c(1,2), function(x){
    ifelse(x < threshold, 1, 0)
  })

gonads_summary_direction <- glms2$coefficients[,c(3,4)] %>%
  apply(., c(1,2), sign)

gonads_summary <- as.data.frame(gonads_summary_signif * gonads_summary_direction)

# Filter only testes vs ovaries contrasts with significant gonad enrichment
gonads_summary$Ovaries_vs_Testes <- factor(
  x = gonads_summary$Ovaries_vs_Testes * (gonads_summary$Gonads_vs_Soma == 1),
  levels = c(0,1,-1),
  labels = c("Unbiased", "Female", "Male")
)

gonads_summary$Gonads_vs_Soma <- factor(
  x = gonads_summary$Gonads_vs_Soma, 
  levels = c(0,1,-1), 
  labels = c("Unbiased", "Gonads", "Soma")
)

table(gonads_summary$Gonads_vs_Soma, gonads_summary$Ovaries_vs_Testes)

write.csv(gonads_summary, file = file.path(newdir, "gonads_DE_summary.csv"))
