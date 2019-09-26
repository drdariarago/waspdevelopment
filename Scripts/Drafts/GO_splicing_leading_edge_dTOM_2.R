## GO testing for enrichment in splicing terms by closeness to splicing nodes
# Initialize script
date()
rm(list=ls())
library(stringr)
# initialize output path
newdir<-file.path(getwd(), "Output/GO_splicing")
dir.create(newdir)
## Load dTOM
load(file = "./Output/dTOMmat/dTOMmat")
## Remove rows containing constitutive genes, remove columns containing facultative transcripts
dTOMmat2<-dTOMmat[grep("_fac", row.names(dTOMmat)),grep("_con", colnames(dTOMmat))]
## Rename columns to get gene names
colnames(dTOMmat)<-str_extract(string = colnames(dTOMmat), pattern = "^[[:alnum:]]*")

## Collapse columns by median of log dTOM
splicing_dTOM<-apply(X = dTOMmat2, 2, function(x){median(log10(x))})
## Export gene list for leading edge enrichment analysis
write.csv(splicing_dTOM, file=file.path(newdir, "splicing_dTOM_median.csv"))

## Export first and last percentiles for standard GO enrichment analysis
write(row.names(splicing_dTOM)[which(splicing_dTOM<quantile(x = splicing_dTOM, probs = 0.01))], file = file.path(newdir, "splicing_dTOM_0.01.txt"))
write(row.names(splicing_dTOM)[which(splicing_dTOM<quantile(x = splicing_dTOM, probs = 0.001))], file = file.path(newdir, "splicing_dTOM_0.001.txt"))
write(row.names(splicing_dTOM)[which(splicing_dTOM<quantile(x = splicing_dTOM, probs = 0.99))], file = file.path(newdir, "splicing_dTOM_0.99.txt"))
write(row.names(splicing_dTOM)[which(splicing_dTOM<quantile(x = splicing_dTOM, probs = 0.999))], file = file.path(newdir, "splicing_dTOM_0.999.txt"))

## Export ranked version of splicing dTOM
ranked_splicing_dTOM<-rank(x = splicing_dTOM, na.last = NA, ties.method = "average")
write.csv(ranked_splicing_dTOM, file=file.path(newdir, "ranked_splicing_dTOM.csv"))

## Collapse columns by median dissimilarity to closest 100 splicing nodes
splicing_dTOM_closest<-apply(X = dTOMmat2, 2, function(x){
  y<-x[order(x, decreasing=F, na.last = NA)][1:100]
  median(y)
  })
## Export gene list for leading edge enrichment analysis
write.csv(splicing_dTOM_closest, file=file.path(newdir, "splicing_dTOM_closest.csv"))

## Collapse columns by 0.1st percentile  of dTOM
splicing_dTOM_quantile<-apply(X = dTOMmat2, 2, function(x){quantile(x = x, probs = 0.001)})
## Export gene list for leading edge enrichment analysis
write.csv(splicing_dTOM_quantile, file=file.path(newdir, "splicing_dTOM_quantile.csv"))

