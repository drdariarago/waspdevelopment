### Extract expression values for gonads and associated whole-organism samples
### Perform filtering for low-expressed genes and genes expressed in few samples

library(magrittr)
library(stringr)
library(plyr)
library(WGCNA)

## Create output directories
newdir<-file.path(getwd(), "Output/gonads_filter")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/gonads_filter")
dir.create(graphdir)

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

# Replace expression values below threshold with threshold value (effectively new zero)
gonad_expr_raw <- apply(X = gonad_expr_raw, MARGIN = c(1,2), function(x){
  ifelse(x<threshold, threshold, x)
})
# Subtract threshold value (rescales to new zero)
gonad_expr_raw <- gonad_expr_raw-threshold
# Filter out zero variance exons
summary(goodGenes(datExpr = t(gonad_expr_raw)))
gonad_expr <- gonad_expr_raw[which(goodGenes(datExpr = t(gonad_expr_raw))),]
# Save as csv
write.csv(x = gonad_expr, file = file.path(newdir, "gonad_exon_expr_postfiltering.csv"))