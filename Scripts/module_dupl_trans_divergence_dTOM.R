# Initialize script
date()
rm(list=ls())
library(plyr)
library(stringr)
# initialize output path
newdir<-file.path(getwd(), "Output/dupl_trans_divergence")
dir.create(newdir)
## Load dTOM
load(file = "./Output/dTOMmat/dTOMmat_test")
dTOMmat<-dTOMmat_test
# create gene list
genelist<-unique(str_extract(string = colnames(dTOMmat), pattern = "[[:alnum:]]*"))
# # create data.frame to store within and between gene topological dissimilarities
trans_dTOM<-data.frame(geneID = genelist, pval = NA, m_within = NA, m_between = NA)
# for every gene, calculate p-value of within/between transcript distace and store median log distances
for (x in genelist){
  print(x)
  pos<-grep()
}
############ outdated code

# load orthologs and genes

# create list with gene names for every OG

# calculate median of log10 pairwise distances within OGs
# use lapply over list of OG/genes

# calculate median of log10 pairwise distances within genes
# create list of genes/transcripts

# calculate median of log10 pairwise distances within random groups of genes