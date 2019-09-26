# Eigenexon based method for deconvoluting transcripts
# Basic concept: transcripts are independent components of a gene's expression. We can collapse exons within a gene in PCs up to a certain threshold of explained variance. Loadings reflect a transcripts' expression levels, while scores the different representation of exons within transcripts
# Important considerations: must collapse negatives to zeroes and choose appropriate correlation measure (kendall or spearmann)
# Must save two files: loadings as new expression data, to be normalized to total gene expression values, while scores as annotation of transcripts of our experiment

# Load packages
date()
library(plyr)
library(stringr)
library(smart)

# Import exon dataset
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON
# filter only rows with samples (exclude means)
rm(list=setdiff(ls(),"nasonia_devtesova.exon.r99"))
# nasdevexon<-nasonia_devtesova.exon.r99[14420:19420,grep("[.][123]",names(nasonia_devtesova.exon.r99))]
nasdevexon<-nasonia_devtesova.exon.r99[grep("Nasvi2EG016484",row.names(nasonia_devtesova.exon.r99)),grep("[.][123]",names(nasonia_devtesova.exon.r99))]
# Exclude gonads 
nasdevexon<-nasdevexon[,-grep("(ovaries|testes)", colnames(nasdevexon))]
# add gene IDs
nasdevexon$geneID<-str_extract(row.names(nasdevexon), "^[^t]*t") # including the t postscript to disambiguate *b and non b gene names
# add exon ID
nasdevexon$exonID<-str_extract(row.names(nasdevexon), "[^t]*$")
# Create matrix of gene expression data
nasdevexonE<-as.matrix(nasdevexon[,1:30])
# Remove negative values
nasdevexonE<-as.data.frame(apply(nasdevexonE, c(1,2), function(x){ifelse(x<0,0,x)}))

# calculate explained variance
$stdev^2
# and proportion of cumulative explained variance
plot(cumsum((a$sdev)^2) / sum(a$sdev^2))

# for every gene apply PCA analysis and save first n principal components in list
# rotation=loadings (combinations of variables used to form PCs), sum ofsquare loadings for each PC=1

library(nsprcomp)
pcobject<-prcomp(nasdevexonE, subset = grep("Nasvi2EG000002t.*", row.names(nasdevexonE)), center = T, tol = 0, scale. = T)

pcobject<-prcomp(t(nasdevexonE[grep("Nasvi2EG000002t.*", row.names(nasdevexonE)),]), center = T, tol = 0, scale. = F)

eigenexons<-lapply(row.names(nasdevexonE), function(x){prcomp(nasdevexonE, subset = grep(x, row.names(nasdevexonE)), center = T, tol = 0)})


# script for rescaling PCs to original data, from http://tinyurl.com/mu9437f to be used with prcomp

pc.use <- 2 # explains 93% of variance
trunc <- pcobject$x[,1:pc.use] %*% t(pcobject$rotation[,1:pc.use])

#and add the center (and re-scale) back to data
if(pcobject$scale != FALSE){
  trunc <- scale(trunc, center = FALSE , scale=1/pcobject$scale)
}
if(is.logical(pcobject$center)==F){
  trunc <- scale(trunc, center = -1 * pcobject$center, scale=FALSE)
}
dim(trunc)

# re-scale each PC in turn, then collapse by averaging exon scores

pcobject<-prcomp(t(nasdevexonE[grep("Nasvi2EG000025t.*", row.names(nasdevexonE)),]), center = T, tol = 0, scale. = F)
pcobject<-nscumcomp(t(nasdevexonE), center = T, tol = 0, scale. = F)
screeplot(pcobject)

pc.use <- 3 # explains 93% of variance
trunc <- sapply(1:pc.use, function(z){
  pc<-pcobject$x[,z] %*% t(pcobject$rotation[,z]) # select PC of interest
  pc<-scale(pc, center = -1 * pcobject$center, scale=FALSE) # scale back exon means
  pc<-if(pcobject$scale != FALSE){trunc <- scale(trunc, center = FALSE , scale=1/pcobject$scale)}else{pc} # remove variance scaling
  pc<-apply(pc, 1, median) # average exon values
  pc
})

# how different are the two scores once averaged?
apply(trunc, 1, median)/apply(nasdevexonE, 2, median)
apply(trunc, 1, sum)/apply(nasdevexonE, 2, sum)
apply(trunc, 1, sum)/apply(nasdevexonE, 2, median)
# how different is each PC versus the genewide average?
trunc/apply(nasdevexonE, 2, median)
# express as proportion of total PC scores
trunc/apply(trunc, 1, sum)
(trunc/apply(trunc, 1, sum))/(1/pc.use)

# how much does every PC increase total expression?
plot(cumsum(apply(trunc, 2, median)))

#and add the center (and re-scale) back to data
if(pcobject$scale != FALSE){
  trunc <- scale(trunc, center = FALSE , scale=1/pcobject$scale)
}
if(is.logical(pcobject$center)==F){
  trunc <- scale(trunc, center = -1 * pcobject$center, scale=FALSE)
}
dim(trunc)
