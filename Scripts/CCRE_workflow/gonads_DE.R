### Analyze gonad data to determine which genes are gonad specific and which ones are male vs female gonad specific rather than shared

library(stringr)
library(magrittr)
# library(reshape2)
# library(plyr)
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
    sex = "contr.treatment",
    soma = "contr.sum")
)
design
