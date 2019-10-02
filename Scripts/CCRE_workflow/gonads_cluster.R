## Cluster gonad expression by eigenexon, normalize and then cluster by CCRE

## Create output directories
newdir<-file.path(getwd(), "Output/gonads_cluster")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/gonads_cluster")
dir.create(graphdir)


# Load eigenexon assignments
eigenexons_assignments <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexons_assignments.csv")

## Load CCRE assignments
CCREs <- read.csv(file = "./Output/CCREs/CCRE_assignments.csv")[,-1]
names(CCREs) <- c("eigenexonID", "CCRE_ID", "Representative_Node")

## Load filtered gonad expression data
gonad_expr <- read.csv(file = "./Output/gonads_filter/gonad_exon_expr_postfiltering.csv", row.names = 1)

## Bin eigenexons


## Bin CCREs


## Perform DE between ovaries and adult females


## Perform DE between testes and male pupae


## Categorize in: male gonads, female gonads, general gonads


## Intersect with other DE data from development


## Categorize proportion of Development DE from gonad expression