## Cluster gonad expression by eigenexon, normalize and then cluster by CCRE

library(FESTA)
library(magrittr)
library(stringr)
library(WGCNA)

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
gonad_expr_eigenexons <- AverageExons(
  data = gonad_expr, 
  spliceID = eigenexons_assignments,
  splicingRatios = T,
  NAcorrection = T
)

rownames(gonad_expr_eigenexons) <- gonad_expr_eigenexons$transcriptID

## Normalize dataset
# Rescale costitutive nodes to 0-1 space by dividing by largest expression value
max_expr <- 
  gonad_expr_eigenexons[,grep(pattern = "[0-3]$", x = colnames(gonad_expr_eigenexons))] %>%
  c(.) %>% unlist(.) %>%
  max(.)

gonad_expr_eigenexons_rescaled_con <- 
  gonad_expr_eigenexons[grep(pattern = "_con$", x = gonad_expr_eigenexons$transcriptID), -(1:2)] %>%
  apply(., c(1,2), function(x){ x / max_expr }) 

gonad_expr_eigenexons_rescaled <- 
  gonad_expr_eigenexons[grep(pattern = "_fac$", x = gonad_expr_eigenexons$transcriptID), -(1:2)] %>%
  rbind(gonad_expr_eigenexons_rescaled_con, .)

gonad_expr_eigenexons_rescaled <- 
  gonad_expr_eigenexons_rescaled[sort(row.names(gonad_expr_eigenexons_rescaled)),]

write.csv(gonad_expr_eigenexons_rescaled, file=file.path(newdir,"gonad_expr_eigenexons.csv"), row.names=F)

# QC filtering of genes with low variance or low expression
# Check for zero variance nodes and nodes with too many missing values
goodCols<-goodSamplesGenes(
  datExpr = transposeBigData(gonad_expr_eigenexons_rescaled),
  minNSamples = 3, minFraction = 1/10)
lapply(goodCols, function(x) {table(x)})
lapply(goodCols, function(x) {prop.table(table(x))})

## Bin CCREs

# Subset only CCRE transcripts
eigenexon_CCRE <- gonad_expr_eigenexons_rescaled[which(row.names(gonad_expr_eigenexons_rescaled) %in% CCREs$eigenexonID),]

# Create match table between row IDs and groups
CCRE_lookup <- merge(
  eigenexon_CCRE[,1 , drop=F], CCREs, 
  by.x = 'row.names', by.y = 'eigenexonID',
  all.x = T, all.y = F)[,-2]

# Create expression matrix for eigenexons that need to be converted to CCREs
row.names(eigenexon_CCRE) <- eigenexon_CCRE$transcriptID
eigenexon_CCRE_expression <- eigenexon_CCRE[,grep(pattern = "[0-3]$", x = colnames(eigenexon_CCRE))]
eigenexon_CCRE_expression <- as.matrix(as.numeric(eigenexon_CCRE_expression))

# Collapse CCREs
CCRE_list <- collapseRows(
  datET = eigenexon_CCRE, 
  rowGroup = CCRE_lookup$CCRE_ID, 
  rowID = CCRE_lookup$transcriptID, 
  connectivityBasedCollapsing = TRUE, 
  method = "MaxMean")

# Merge with non-CCRE transcripts


# Save
