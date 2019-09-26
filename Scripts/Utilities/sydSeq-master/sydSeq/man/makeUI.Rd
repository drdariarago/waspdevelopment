\name{makeUI}
\Rdversion{1}
\alias{makeUI}
\title{Make UI annotation}
\description{
Make a UI annotation from GRanges object.
}
\usage{
makeUI(simpleGRanges, gFeature)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{simpleGRanges}{An object of type \code{GRanges} which has been simplified}
  \item{gFeature}{a genome annotation of type \code{gFeature}}
}
  \details{

Method makes UI definition.
}
\value{
An object of type \code{GRanges} containing UI constitutive exons.
}
\author{ Ellis Patrick }
\examples{
\dontrun{
hsapiensEnsembl <-makeTranscriptDbFromBiomart(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
saveFeatures(hsapiensEnsembl, file="hsapiensEnsembl.sqlite")
simpleEnsembl <- simplifyGFeatures(hsapiensEnsembl)
annoUI <- makeUI(simpleEnsembl,hsapiensEnsembl)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
