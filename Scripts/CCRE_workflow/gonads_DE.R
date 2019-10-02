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

female_expr <- gonad_expr[, grepl("adult|ovaries", colnames(gonad_expr))]

# Convert to expression list
female_expr_list<-list(
  E = as.matrix(female_expr), 
  genes = as.character(female_expr),
  targets=names(female_expr))

# Converting to LIMMA EList format
female_expr_list<-new("EList", female_expr_list)

# Create factors
soma <- as.factor(ifelse(grepl("testes|ovaries", colnames(female_expr)), "Gonads", "Soma")) 

# Create design matrix
design <- model.matrix(
  ~ soma, 
  contrasts.arg=list(
    soma="contr.treatment")
)

