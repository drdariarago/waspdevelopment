\name{exClust}
\Rdversion{1}
\alias{exClust}
\title{Estimates constitutive exons using both annotation and data}
\description{
Use clustering to estimate constitutive exons.
}
\usage{
exClust(simpleGRanges, cut = 1.5)
}
\arguments{
\item{simpleGRanges}{An object of type \code{GRanges} which has been simplified and had a counts matrix merged.}
  \item{cut}{The level to cut clustering tree.}
}
  \details{
Method clusters exons and cuts clustering tree at value \code{cut}.
}
\value{
An object of type \code{GRanges} containing the data specific constitutive exons.
}
\author{ Ellis Patrick }
\examples{
\dontrun{
hsapiensEnsembl <-makeTranscriptDbFromBiomart(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
saveFeatures(hsapiensEnsembl, file="hsapiensEnsembl.sqlite")
simpleEnsembl = simplifyGFeatures(hsapiensEnsembl)

###  Read in Data from bam files
bamDir = '.'
bamFls <- list.files(bamDir, "bam$", full=TRUE)
names(bamFls) <- sub("\\..*", "", basename(bamFls))
counter <- function(fl, gnModel)
{
    aln <- readGappedAlignments(fl)
	seqlevels(aln) <- sub("chr","", seqlevels(aln))
    strand(aln) <- "*" # for strand-blind sample prep protocol
    hits <- countOverlaps(aln, gnModel)
    counts <- countOverlaps(gnModel, aln[hits==1])
    names(counts) <- names(gnModel)
    counts
}
counts <- sapply(bamFls, counter, simpleEnsembl)

simpleEnsembl= mergeCounts(counts,simpleEnsembl)
save(simpleEnsembl,file = 'simpleEnsembl.RData')
exClustAnno = exClust(simpleEnsembl,cut = 2)
}
}
