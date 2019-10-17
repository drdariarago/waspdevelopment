## Cluster-level summary of gonad bias statistics

#library(stringr)
library(magrittr)
library(tidyverse)

## Create output directories
newdir<-file.path(getwd(), "Output/ClusterStats_gonads")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/ClusterStats_gonads")
dir.create(graphdir)

## Import cluster-level data

clusterdata <- read_csv(file = "./Output/Results_compiler_CCRE/clusterdata_full.csv")[,-1]

## Cross-compare clusters with DE, DC and gonad expression
clusterbias <- clusterdata %>%
  select(., 
    clusterID = clusterID,
    devsexbias = cluster_devsexbias,
    nGenes = nGenes,
    DI = DiffIntegrated, 
    testes_enrichment = testes_enrichment,
    ovaries_enrichment = ovaries_enrichment,
    gonads_enrichment = gonads_enrichment
  ) %>%
  # Decompress devsexbias
  tidyr::extract(
    data = ., col = devsexbias, 
    into = c("emb10", "emb18", "lar51", "pupyel", "adult"), 
    regex = "(.)(.)(.)(.)(.)"
  ) %>%
  # Stack stages in new factor
  pivot_longer(
    ., cols = 2:6, 
    names_to = "Stage", 
    names_ptypes = list( 
      Stage = factor(levels =  c("emb10", "emb18", "lar51", "pupyel", "adult"))),
    values_to = "sexbias"
  ) %>%
  group_by(.data = ., clusterID
  ) %>%
  # Recode sex bias as male, female, unbiased or both
  summarise(.data = .,
    nGenes = first(nGenes),
    DI = first(DI),
    testes = first(testes_enrichment),
    ovaries = first(ovaries_enrichment),
    gonads = first(gonads_enrichment),
    sexbias = if(any(sexbias == "m") & all(sexbias != "f")){
      "male"
    } else if(any(sexbias == "f") & all(sexbias != "m")){
      "female"
    } else if(any(sexbias =="f") & any(sexbias == "m")) {
      "both"
    } else {"unbiased"}
  ) %>%
  # convert NAs in testes and ovary bias to 1 (unbiased), then binarize
  mutate(.data = ., 
    sexbias = as.factor(sexbias),
    testes = ifelse(is.na(testes), 1, testes) %>% is_less_than(0.05),
    ovaries = ifelse(is.na(ovaries), 1, ovaries) %>% is_less_than(0.05),
    gonads= ifelse(is.na(gonads), 1, gonads) %>% is_less_than(0.05)
  )
      
clusterbias

## Select clusters with bot DE and DC

clusterbias %>%
  filter(.,
    grepl('adult', DI),
    sexbias != 'unbiased'
  ) %>%
  view()

## Count malebias genes with and without testes enrichment
clusterbias %>%
  group_by(sexbias == "male", testes
    ) %>%
  summarise(.data = .,
    nGenes = sum(nGenes),
    nClusters = n()
  )

## Count femalebiased genes with and without testes enrichment
clusterbias %>%
  group_by(sexbias == "female", ovaries
  ) %>%
  summarise(.data = .,
    nGenes = sum(nGenes),
    nClusters = n()
  )

