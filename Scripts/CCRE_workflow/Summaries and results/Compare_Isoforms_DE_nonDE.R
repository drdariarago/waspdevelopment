### Compare the number of exons and splicing nodes in DE vs non DE genes

## load libraries
library(tidyverse)
library(magrittr)
library(ggplot2)

## Load transcript-level data
transcriptdata <- 
  # Import transcript annotation
  read_csv(
    file = "./Output/Results_compiler_CCRE/transcriptdata_full.csv")[,-1]



## Count exons for DE vs non DE nodes

splicedata <- 
  transcriptdata %>%
  select(.,
    "eigenexonID", "geneID", "devsexbias"
  ) %>%
  group_by(., geneID) %>%
  mutate(.,
    DE_gene = grepl("_con", eigenexonID) & devsexbias != ".....",
    DE_trans = grepl("_fac", eigenexonID) & devsexbias != "....."
  ) %>%
  summarize(.,
    n_transcripts = n(),
    spliced = any(grepl("_fac", eigenexonID)),
    DE_gene = any(DE_gene),
    DE_trans = any(DE_trans)
  )

# Report average transcript number of transcripts for spliced and un-spliced genes..
group_by(splicedata, spliced) %>%
  summarize(
    n_transcripts = sum(n_transcripts),
    p_transcripts = n_transcripts/n()
  )

# Plot distribution of transcript number for genes with and without spicing
splicedata %>%
  # filter(., n_transcripts > 1) %>%
  ggplot(., aes(x = n_transcripts, fill = DE_trans)) +
  geom_histogram(position = 'dodge') +
  ylim(c(0,4000))

splicedata %>%
#  filter(., n_transcripts > 1) %>%
  ggplot(., aes(x = n_transcripts, fill = DE_gene)) +
  geom_histogram(position = 'dodge') +
  ylim(c(0,4000))

# Do any sex-biased genes have only 1 transcript?
filter(splicedata, n_transcripts == 1) %>%
  group_by(., DE_gene) %>%
  summarize(
    n_genes = n()
  )

## 2D proportion plot using only genes with >1 transcript

counts_multiexon <- filter(splicedata, n_transcripts > 1) %>%
  mutate(., 
    DE_trans = factor(
      DE_trans, 
      levels = c("TRUE", "FALSE"), 
      labels = c("Switch", "noSwitch")),
    DE_gene = factor(
      DE_gene, 
      levels = c("TRUE", "FALSE"), 
      labels = c("Biased", "Unbiased"))
    
  ) %>%
  with(., table(switching = DE_trans, bias = DE_gene)) 

fisher.test(counts_multiexon)

mosaicplot(counts_multiexon, color = T)

