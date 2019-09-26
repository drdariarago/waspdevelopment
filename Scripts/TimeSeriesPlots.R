### Make plot for Nasonia cluster DE across development

library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)
library(viridis)

graphdir <- file.path(getwd(),"Graphics/TimeSeriesPlots")
dir.create(graphdir)

# Import cluster data
clusterdata = read.csv(file = "./Output/Results_compiler_CCRE/clusterdata_full.csv", row.names = 1)

# Summarize by expression pattern
# also include: raw gene counts, transcript counts, cluster count
sexbias <- ddply(
  clusterdata, 
  .variables = .(cluster_devsexbias), 
  summarize,
  nClusters = length(nGenes),
  nGenes = sum(nGenes),
  nTranscripts = sum(nTranscripts)
)

# Convert character string to matrix
bias_mat <- t(matrix(unlist(str_extract_all(sexbias$cluster_devsexbias, pattern = ".")), nrow = 5))
# Convert matrix to ternary data (1/0/-1)
bias_mat <- apply(bias_mat, c(1,2), function(x){
  if (x=="m") {
    1
  } else if (x=="f") {
    -1
  } else {0}
})

colnames(bias_mat) <- c("Early Embryo", "Late Embryo", "Larva", "Pupa", "Adult")

# Merge matrix with annotation, and add stage
sexbias <- cbind(bias_mat, sexbias)

# Melt matrix entries as data, to obtain flat file with stage column
sexbias <- melt(data = sexbias, measure.vars = c(1:5))
names(sexbias)[5:6] <- c("Stage", "SexBias")

# Reorder cluster_devsexbias according to number of clusters
sexbias$cluster_devsexbias <- factor(
  x = sexbias$cluster_devsexbias, 
  levels = unique(sexbias$cluster_devsexbias[order(-sexbias$nTranscripts)]))
# Put unbiased clusters last
sexbias$cluster_devsexbias <- factor(
  x = sexbias$cluster_devsexbias,
  levels = c(levels(sexbias$cluster_devsexbias)[-1], levels(sexbias$cluster_devsexbias)[1])
)

# Plot as line-plot with lines weighted by cluster size, annotate number of genes and trancripts
ggplot(data = sexbias, 
  mapping = aes(x = Stage, y = SexBias, col = nTranscripts)) +
  geom_line(group = "cluster_devsexbias", size = 1.5) +
  facet_wrap(~cluster_devsexbias, ncol = 4) +
  geom_text(aes(x = 1, y = 0.8, label = nClusters, col = 0.01)) +
  scale_color_viridis(name = "Number of \nTranscripts", option = "D",
    trans='log10', breaks=c(50,5E2,5E3,2E4), limits = c(50, 20000)) +
  scale_y_continuous(name = "Sex Bias", breaks = c(-1, 0, 1), limits = c(-1.5, 1.5), 
    labels = c("Female", "Unbiased", "Male")) +
  scale_x_discrete(name = "Developmental Stage", 
    labels = c("Early\nEmbryo", "Late\nEmbryo", "Larva", "Pupa", "Adult")) +
  theme(
    # # Remove numeric annotation from y axis
    # axis.ticks.y = element_blank(),
    # axis.text.y = element_blank(),
    # Remove labels from facets
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    # Remove grids
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(),
    # Set y grid (bias)
    panel.grid.major.y = element_line(color = "lightgrey"),
    # Set background
    panel.background = element_rect(fill = "white")
  )

ggsave(filename = file.path(graphdir, "TimeSeries.pdf"), device = "pdf", 
  units = 'mm', width = 297)
ggsave(filename = file.path(graphdir, "TimeSeries.png"), device = "png", 
  units = 'mm', width = 297)
