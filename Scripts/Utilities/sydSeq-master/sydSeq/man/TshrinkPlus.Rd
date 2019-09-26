\name{TshrinkPlus}
\Rdversion{1}
\alias{TshrinkPlus}
\title{Perform a moderated t-test using additional information}
\description{
TshrinkPlus is a performs a moderated t-test using additional information for normalised RNA-Seq counts.
}
\usage{
TshrinkPlus(Data, classes, meanPlus = NULL, variancePlus = NULL, verbose = TRUE,nn = 0.7)}
\arguments{
\item{Data}{A matrix of normlised counts with genes as row names and samples as column names.}
  \item{classes}{A vector containing two class labels for each sample. The names of which are the sample sample names.}
\item{meanPlus}{A vector of average gene expression from the additional data set named accordingly.}
\item{variancePlus}{A vector of gene variances from the additional data set named accordingly.}
\item{verbose}{Should warnings be printed}
\item{nn}{A parameter for kernal smoothing}
}
  \details{
Method performs a moderated t-test. The method uses additional information to fit the common variance curve if either \code{variancePlus} is not NULL or if both \code{variancePlus} and \code{meanPlus} are not NULL.
}
\value{
A list containg the pvalues, degrees of freedom and shrinkage parameter for the tests.
}
\author{ Ellis Patrick }
\examples{

set.seed(123)
# Generate true mean and variance
m = exp(rgamma(10000,5,1))
v = rnorm(10000)
v = pmax(m+1/40*m^2*exp(v),m)
# simulate "normalised counts"
Data = matrix(rnorm(10000*8,m,sqrt(v)),10000,8)
colnames(Data) = paste('s',1:8,sep ='')
rownames(Data) = paste('g',1:10000,sep ='')

# simulate additional data

Additional.Data = matrix(rnorm(10000*8,m,sqrt(pmax(m+1/40*m^2*exp(v),m))),10000,8)
rownames(Additional.Data) = paste('g',1:10000,sep ='')

# check mean vs variance
plot(apply(Data,1,mean),apply(Data,1,var),log = 'xy',xlab = 'Mean',ylab = 'Variance')
abline(0,1)

# calculate additional mean and variance vector
meanPlus = apply(Additional.Data,1,mean)
variancePlus = apply(Additional.Data,1,var)

# filter out genes with small means
variancePlus = variancePlus[meanPlus>5]
meanPlus = meanPlus[meanPlus>5]
Data = Data[rowMeans(Data)>10,]

# make class vector
classes = rep(1:2,c(4,4))
names(classes) = colnames(Data)

# test using tshrink and tshrinkplus
tshrinkplus = TshrinkPlus(Data = Data,classes = classes ,meanPlus,variancePlus)
tshrink = TshrinkPlus(Data,classes = classes)
}
