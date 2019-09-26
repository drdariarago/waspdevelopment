\name{simplifyGFeatures}
\Rdversion{1}
\alias{simplifyGFeatures}
\title{Simplifies a GFeatures object}
\description{
Simplifies a GFeatures object removing regions that overlap multiple genes and splits exons into unique regions.
}
\usage{
simplifyGFeatures(gFeature)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{gFeature}{An object of type \code{gFeature} }
}
\details{
.
}
\value{
An object of type \code{GRanges} }
\author{ Ellis Patrick }
\examples{
\dontrun{
hsapiensEnsembl <-makeTranscriptDbFromBiomart(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
saveFeatures(hsapiensEnsembl, file="hsapiensEnsembl.sqlite")
simpleEnsembl = simplifyGFeatures(hsapiensEnsembl)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
