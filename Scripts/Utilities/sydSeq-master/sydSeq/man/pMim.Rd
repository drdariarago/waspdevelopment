\name{pMim}
\Rdversion{1}
\alias{pMim}
\title{Infer function of miRNA in RNA-Seg and miRNA-Seq data}
\description{
pMim jointly ranks miRNA and pathways with respect to a condition of interest.
}
\usage{
pMim(DataMi, DataG, classes, targets, pathways, Zmi=NULL, Zg=NULL, stat = 'cor', corP = 'Fisher', nInt = 2, randomTarg=FALSE, randomPath=FALSE, verbose == TRUE)}
\arguments{
\item{DataMi}{A matrix of normalised miRNA expression values.}
\item{DataG}{A matrix of normalised gene expression values.}
  \item{classes}{A vector containing two class labels for each sample. The names of which are the sample names.}
\item{targets}{A named list of the genes targeted by the named miRNA.}
\item{pathways}{A named list of the genes in each named pathway.}
\item{Zmi}{A vector of z-scores from a DE test on miRNA.}
\item{Zg}{A vector of z-scores from a DE test on genes.}
\item{stat}{A string describing which gene measure to use, 'cor' or 'de'.}
\item{nInt}{The minimum number of genes in each mir-pathway.}
\item{randomTarg}{Logical of whether gene target information should be randomised}
\item{randomPath}{Logical of whether pathway information should be randomised}
\item{verbose}{Should warnings be printed}
}
  \details{
Method ranks mir-pathways with respect to a condition of interest.
}
\value{
A list containing a formatted results matrix and the unformatted scores for each mir-pathway. Small scores provide evidence that the mir-pathway is interesting.
}
\author{ Ellis Patrick }
\examples{
### Generate Data
nMir = 10
nGene = 1000
nPath = 100
n = 20
DataMi = matrix(rnorm(nMir*n,1000,sqrt(1000)),nMir,n)
rownames(DataMi) = paste('mir',1:nMir,sep = '-')
DataG = matrix(rnorm(nGene*n,1000,sqrt(1000)),nGene,n)
rownames(DataG) = paste('gene',1:nGene,sep = '-')

targets = sample(rownames(DataG),nMir*100,replace = TRUE)
targets = split(targets,rep(rownames(DataMi),rep(100,nMir)))
targets = lapply(targets,unique)

pathways = sample(rownames(DataG),nPath*100,replace = TRUE)
pathways = split(pathways,rep(paste('path',1:nPath,sep = '-'),rep(100,nPath)))
pathways = lapply(pathways,unique)

classes = rep(0:1,c(10,10))
names(classes) = colnames(DataMi) = colnames(DataG) = paste('sample',1:n,sep = '-')

### Run pMim
output = pMim(DataMi, DataG, classes, targets, pathways)
}
