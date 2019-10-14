### Create table of expressed vs sex-biased genes and transcripts at each stage
## Annotate proportion of genes which also have testes/ovary expression

# Initialize output path
newdir<-file.path(getwd(),"Output/GeneExpressionTable")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/GeneExpressionTable")
dir.create(graphdir)

# Load libraries
library(tidyverse)
library(magrittr)

## Import normalized dataset of gene expression per stage
nasdevexpression <- 
  read_csv(file = "./Output/CCREs/nasdevgeneexon.csv") %>%
  column_to_rownames(., "X1")

## Annotate transcript and gene IDs
nasdevexpression$transcriptID <- row.names(nasdevexpression)
nasdevexpression$geneID <- 
  row.names(nasdevexpression) %>%
  str_extract(pattern = "^[^_]*")

# Annotate stage, sex and replicate
tidyexpression <- pivot_longer(
  data = nasdevexpression, 
  cols = tidyselect::matches(match = "*[0-3]$"), 
  names_to = c("Stage", "Sex", "Replicate"), 
  names_pattern = "(.*)_(.*)\\.([0-3])",
  names_ptypes = list(
    Sex = factor(levels = c("female", "male")),
    Stage = factor(levels = c("emb10", "emb18", "lar51", "pupyel", "adult"))
  ),
  values_to = "Expression"
)

## Calculate how many transcripts have expression >0 for 2/3 replicates



## Bin by stage (at least one of the two sexes is expressed)

## Bin by gene (at least one of their transcripts are expressed)


# nodedata <- read.csv(
#   file = "./Output/Results_compiler_CCRE/nodedata_full.csv", 
#   row.names = 2)[,-1]
