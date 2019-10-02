## Cluster gonad expression by eigenexon, normalize and then cluster by CCRE

library(magrittr)
library(stringr)
library(dplyr)

library(WGCNA)
library(FESTA)

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
names(CCREs) <- c("transcriptID", "CCRE_ID", "Representative_Node")

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

# Replace each CCRE only with its representative node, based on table used for main experiment
eigenexon_CCRE_representative <-
  merge(
    CCREs, eigenexon_CCRE, by.x = 'transcriptID', by.y = 'row.names') %>% 
  filter(., Representative_Node == TRUE) %>%
  set_rownames(., .$transcriptID) %>%
  select(., grep(pattern = "[0-3]$", x = colnames(.))
  )

# Merge with non-CCRE transcripts
gonad_expr_CCRE <- 
  gonad_expr_eigenexons_rescaled[-which(row.names(gonad_expr_eigenexons_rescaled) %in% CCREs$eigenexonID),] %>%
  rbind(., eigenexon_CCRE_representative)
  
# Save

write.csv(x = gonad_expr_CCRE, file = file.path(newdir, "gonad_expr_CCRE.csv"))