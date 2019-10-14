### Create table of expressed vs sex-biased genes and transcripts at each stage
## Annotate proportion of genes which also have testes/ovary expression

# Initialize output path
newdir<-file.path(getwd(),"Output/GeneExpressionTable")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/GeneExpressionTable")
dir.create(graphdir)

# Load libraries
library(tidyverse)
library(magrittr)

## Import normalized dataset of gene expression per stage
nasdevexpression <- 
  read_csv(file = "./Output/CCREs/nasdevgeneexon.csv") %>%
  column_to_rownames(., "X1")

## Annotate transcript and gene IDs
nasdevexpression$transcriptID <- row.names(nasdevexpression)
nasdevexpression$geneID <- 
  row.names(nasdevexpression) %>%
  str_extract(pattern = "^[^_]*")

# Annotate stage, sex and replicate
tidyexpression <- pivot_longer(
  data = nasdevexpression, 
  cols = tidyselect::matches(match = "*[0-3]$"), 
  names_to = c("Stage", "Sex", "Replicate"), 
  names_pattern = "(.*)_(.*)\\.([0-3])",
  names_ptypes = list(
    Sex = factor(levels = c("female", "male")),
    Stage = factor(levels = c("emb10", "emb18", "lar51", "pupyel", "adult"))
  ),
  values_to = "Expression"
)

## Calculate how many transcripts have expression >0 for 2/3 replicates in at least one sex

tidyexpression_2 <- 
  group_by(.data = tidyexpression, transcriptID, Stage, Sex) %>%
  summarise(., 
    geneID = first(geneID),
    Expressed = sum(Expression > 0) > 2
  ) %>% 
  summarise(., 
    geneID = first(geneID),
    Expressed = max(Expressed)
  )

## Report number of expressed transcripts per stage
table(tidyexpression_2$Expressed, tidyexpression_2$Stage)

## Merge with differential expression from main analysis
# This step helps us calculate proportions and compare with gene-level analyses

nodedata <- 
  # Import dataset
  read_csv(
    file = "./Output/Results_compiler_CCRE/nodedata_full.csv")[,-1] %>%
  select(., 
    nodeID = nodeID,
    transcriptID = eigenexonID,
    geneID = geneID, 
    devsexbias = devsexbias
  ) %>%
  # Decompress devsexbias
  tidyr::extract(
    data = ., col = devsexbias, 
    into = c("emb10", "emb18", "lar51", "pupyel", "adult"), 
    regex = "(.)(.)(.)(.)(.)"
  ) %>%
  # Stack stages in new factor
  pivot_longer(
    ., cols = -ends_with("ID"), 
    names_to = "Stage", 
    names_ptypes = list( 
      Stage = factor(levels =  c("emb10", "emb18", "lar51", "pupyel", "adult"))),
    values_to = "sexbias"
  ) 

# Merge datasets

biasdata <- merge(
  x = nodedata, y = tidyexpression_2, 
  by = c("transcriptID", "Stage", "geneID"), 
  all.x = T, all.y = F
) %>%
  mutate(.data = .,
    F_bias = sexbias == "f",
    M_bias = sexbias == "m"
  ) %>%
  group_by(.data = .,
    nodeID, transcriptID, geneID, Stage)

# Create report table