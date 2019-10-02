## Cluster gonad expression by eigenexon, normalize and then cluster by CCRE

library(FESTA)
library(magrittr)
library(stringr)

## Create output directories
newdir<-file.path(getwd(), "Output/gonads_cluster")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/gonads_cluster")
dir.create(graphdir)

# Load eigenexon assignments
eigenexons_assignments <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexons_assignments.csv")
names(eigenexons_assignments)[7] <- "transcriptID"
eigenexons_assignments <- eigenexons_assignments[,-4]

## Load CCRE assignments
CCREs <- read.csv(file = "./Output/CCREs/CCRE_assignments.csv")[,-1]
names(CCREs) <- c("eigenexonID", "CCRE_ID", "Representative_Node")

## Load filtered gonad expression data
gonad_expr <- read.csv(file = "./Output/gonads_filter/gonad_exon_expr_postfiltering.csv", row.names = 1)

## Filter only transcripts present in main experiment
gonad_expr <- gonad_expr[which(row.names(gonad_expr)%in%eigenexons_assignments$exonID),]
## Filter only transcripts that are present in gonads/soma
eigenexons_assignments <- eigenexons_assignments[which(eigenexons_assignments$exonID%in%row.names(gonad_expr)),]

## Remove isoforms that have no main gene expression
no_constitutive <- setdiff(
  eigenexons_assignments[eigenexons_assignments$constitutive=='facultative', "geneID"],
  eigenexons_assignments[eigenexons_assignments$constitutive=='constitutive', "geneID"]
)

no_constitutive_exons <- str_extract(row.names(gonad_expr), "^[^t]*") %>%
  is_in(., no_constitutive) %>%
  which(. == T) 
gonad_expr <- gonad_expr[-no_constitutive_exons, ]

eigenexons_assignments <- eigenexons_assignments[-which(eigenexons_assignments$geneID %in% no_constitutive),]


## Bin eigenexons
debugonce(FESTA::AverageExons)
gonad_expr_eigenexons <- AverageExons(
  data = gonad_expr, 
  spliceID = eigenexons_assignments,
  splicingRatios = T,
  NAcorrection = T
)

## Bin CCREs


## Perform DE between ovaries and adult females


## Perform DE between testes and male pupae


## Categorize in: male gonads, female gonads, general gonads


## Intersect with other DE data from development


## Categorize proportion of Development DE from gonad expression