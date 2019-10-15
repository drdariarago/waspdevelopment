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

transcriptdata <- 
  # Import dataset
  read_csv(
    file = "./Output/Results_compiler_CCRE/transcriptdata_full.csv")[,-1] %>%
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
  x = transcriptdata, y = tidyexpression_2, 
  by = c("transcriptID", "Stage", "geneID"), 
  all.x = T, all.y = F
) %>%
  mutate(.data = .,
    F_bias = sexbias == "f",
    M_bias = sexbias == "m"
  ) %>%
  group_by(.data = .,
    nodeID, transcriptID, geneID, Stage)

# Create report table for transcripts
transcript_counts <- group_by(.data = biasdata, transcriptID, Stage) %>%
  summarise(.data = .,
    Expressed = max(Expressed),
    F_bias = max(F_bias),
    M_bias = max(M_bias)
  ) %>%
  group_by(.data = ., Stage) %>%
  summarise(.data = ., 
    Expressed = sum(Expressed),
    F_bias = sum(F_bias),
    M_bias = sum(M_bias)
  )

transcript_counts
# Calculate percent DE genes of the total genes expressed in that stage
transcript_counts[,4:3]/transcript_counts$Expressed * 100


## Create report table for genes

gene_counts <- group_by(.data = biasdata, geneID, Stage) %>%
  summarise(.data = .,
    Expressed = max(Expressed),
    F_bias = max(F_bias),
    M_bias = max(M_bias)
  ) %>%
  group_by(.data = ., Stage) %>%
  summarise(.data = ., 
    Expressed = sum(Expressed),
    F_bias = sum(F_bias),
    M_bias = sum(M_bias)
  )

gene_counts
# Calculate percent DE genes of the total genes expressed in that stage
gene_counts[,4:3]/gene_counts$Expressed * 100

## Import gonad expression and compare with stage-by-stage gene expression
gonaddata <- 
  read_csv(file = "./Output/gonads_DE/gonads_DE_summary.csv") %>%
  set_colnames(., c("nodeID", "Gonad_vs_Soma", "Ovaries_vs_Testes", "gonad_bias"))

biasdata_gonads <- 
  merge(
    x = gonaddata[,c("nodeID", "gonad_bias")], y = transcriptdata,
    all.x = T, all.y = T, 
    by = "nodeID"
  ) %>%
  filter(., !is.na(sexbias)) %>%
  mutate(.,
    gonad_bias = ifelse(is.na(gonad_bias), 0, gonad_bias)
  )

transcript_counts_gonads <- 
  group_by(.data = biasdata_gonads, transcriptID, Stage) %>%
  summarise(.data = .,
    M_testes = sexbias == "m" & gonad_bias == "Testes",
    M_soma = sexbias == "m" & gonad_bias != "Testes",
    F_ovaries = sexbias == "f" & gonad_bias == "Ovaries",
    F_soma = sexbias == "f" & gonad_bias != "Ovaries"
  ) %>%
  group_by(.data = ., Stage) %>%
  summarise(.data = ., 
    M_testes = sum(M_testes, na.rm = T),
    M_soma = sum(M_soma, na.rm = T),
    F_ovaries = sum(F_ovaries, na.rm = T),
    F_soma = sum(F_soma, na.rm = T)
  )

transcript_counts_gonads


gene_counts_gonads <- 
  group_by(.data = biasdata_gonads, geneID, Stage) %>%
  summarise(.data = .,
    M_testes = max(sexbias == "m" & gonad_bias == "Testes"),
    M_soma = max(sexbias == "m" & gonad_bias != "Testes"),
    F_ovaries = max(sexbias == "f" & gonad_bias == "Ovaries"),
    F_soma = max(sexbias == "f" & gonad_bias != "Ovaries")
  ) %>%
  group_by(.data = ., Stage) %>%
  summarise(.data = ., 
    M_testes = sum(M_testes, na.rm = T),
    M_soma = sum(M_soma, na.rm = T),
    F_ovaries = sum(F_ovaries, na.rm = T),
    F_soma = sum(F_soma, na.rm = T)
  )

gene_counts_gonads
