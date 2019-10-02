### Gonad analyses

library(magrittr)
library(stringr)
library(plyr)

## Create output directories
newdir<-file.path(getwd(), "Output/DE_gonads")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/DE_gonads")
dir.create(graphdir)


# Load eigenexon assignments
eigenexons_assignments <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexons_assignments.csv")

## Load CCRE assignments
CCREs <- read.csv(file = "./Output/CCREs/CCRE_assignments.csv")[,-1]
names(CCREs) <- c("eigenexonID", "CCRE_ID", "Representative_Node")


## Import expression for gonads and respective whole-animal stages
# Import exon-level expression
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/2015_05_19_Justin/nasonia_devtesova.exon.r99.score2")
# Select only gonads and respective whole-animal stages (male pupae, female adults)
gonad_expr_raw <- 
  c("EXON|testes.[1-3]|ovaries.[1-3]|adult_female.[1-3]|pupyel_male.[1-3]") %>% 
  grep(pattern = ., names(nasonia_devtesova.exon.r99), value = T) %>%
  nasonia_devtesova.exon.r99[,.]

# remove duplicated exons add gene IDs
gonad_expr_raw<-gonad_expr_raw[which(duplicated(gonad_expr_raw$EXON)==F),]
row.names(gonad_expr_raw) <- gonad_expr_raw$EXON
gonad_expr_raw <- gonad_expr_raw[,-1]
# gonad_expr_raw$geneID<-str_extract(row.names(gonad_expr_raw), "^[^t]*") 
# including the t postscript to disambiguate *b and non b gene names

## Normalize exon expression
# Select same threshold as main experiment
threshold = 1.66
# Remove exons that not expressed in at least 2 out of 3 samples
biorep <- as.factor(str_extract(string = names(gonad_expr_raw), pattern = "^[^.]*"))
t_gonad_expr_raw <- data.frame(biorep, t(gonad_expr_raw)) # remove subselection of exons after debugging
l_gonad_expr_raw <- dlply(.data = t_gonad_expr_raw, .variables = .(biorep))

expressed <- sapply(l_gonad_expr_raw, function(x){
  apply(X = x, MARGIN = 2, FUN = function(y){
    sum(y>threshold)>1
  })
})
expressed <- expressed[-1,]
row.names(expressed) <- row.names(gonad_expr_raw)
# Select only exons expressed in at least one biological condition
pass <- apply(X = expressed[-1,], MARGIN = 1, FUN = function(x){any(x==T)})
pass <- names(pass)[which(pass==T)]
# Subset exon dataset only to exons expressed in at least one biological condition
gonad_expr_raw <- gonad_expr_raw[which(row.names(gonad_expr_raw)%in%pass),]
# How many exons are left?
nrow(gonad_expr_raw)


## Bin eigenexons


## Bin according to CCREs



## Perform DE between ovaries and adult females


## Perform DE between testes and male pupae


## Categorize in: male gonads, female gonads, general gonads


## Intersect with other DE data from development


## Categorize proportion of Development DE from gonad expression