## Check for emergence of motifs that lead to sex-specific patterns of gene expression from basal unbiased states
## Or rather, check for new motifs/tfs that make groups of paralogs diverge
date()
rm(list=ls())
newdir<-file.path(getwd(), "Output/Paralog_coevolution")
graphdir<-file.path(getwd(), "Graphics/Paralog_coevolution")
dir.create(newdir)
dir.create(graphdir)
library(boot)
library(fdrtool)
library(ggplot2)
library(plyr)
library(tidyr)
library(reshape)
library(stringr)
library(cooccur)
library(igraph)


# Load compiled dataset
clusterdata_compiled<-read.csv(file="./Output/Results_compiler/clusterdata_compiled.csv")[,-1]
# Select only transcription nodes and remove genes without OG assignment (avoids assignments of genes to multiple clusters, we are focusing on TRANSCRIPTION alone)
# Also remove multiple occurrences of the same OG/cluster (using presence/absence only, enables filtering of OGs with only 1 cluster below)
genedata<-droplevels(unique(na.exclude(clusterdata_compiled[which(grepl("con", clusterdata_compiled$eigenexonID)),c("clusterID","odb8_og_id")])))
# 9261 OGs, 100 clusters

# Subselect only OGs with more than one member (restrict to genes with paralogs present in more than one group)
OGcounts<-(table(genedata$odb8_og_id)>1)
genedata<-droplevels(genedata[which(genedata$odb8_og_id%in%names(OGcounts)[which(OGcounts==T)]),])
length(levels(droplevels(genedata$clusterID))) #310 clusters
length(levels(droplevels(genedata$odb8_og_id))) # 1398 gene families

# Plot data as a table of ortholog groups by clusters (presence/absence only)
OG_cluster_data<-as.matrix(table(genedata$odb8_og_id, genedata$clusterID))

# Perform co-occurrence analysis (from ecological data analyses)
# Note: this method assumes equally random probabilities of assignment to each cluster, could be non-conservative if duplicates tend to retain expression patterns (genes will also retain expression patterns, leading to inflated co-occurrences regardless of co-evolution, possibly not relevant since we only consider genes that occur in multiple clusters and have therefore broken free of the initial constraint)

# Remove highly coduplicating transposons
transposons<-c(
  "EOG8CZDZX", "EOG8WWV1R", "EOG8DV81S","EOG83XXMT","EOG8D2962","EOG8CNT6C", "EOG8FJBS3","EOG8CG2R8", "EOG8VDSFF", "EOG893614","EOG8WDGTK", "EOG8HTC74","EOG8Z91BB","EOG8VX4M0", "EOG8B2WDC", "EOG8QNQD6", "EOG80ZTCF", "EOG8935ZZ", "EOG84XN1N", "EOG88GZJM", "EOG8M3CQ4", "EOG8MGVNF", "EOG8QJV2S", "EOG8Q5C4B", "EOG834ZP6", "EOG8DRCV3", "EOG8RV579", "EOG8ZGRV6", "EOG812PMP", "EOG841SSS", "EOG896286", "EOG8B8MV0", "EOG8D2968", "EOG8FXTNG", "EOG8J6V5R", "EOG8PVRGG", "EOG8X9B0N", "EOG8ZCWMQ","EOG8NZX77","EOG8MSGGG")

OG_cluster_data <- OG_cluster_data[-which(row.names(OG_cluster_data)%in%transposons),]

# Computationally intensive (~20 min on laptop without pairs, 18 hrs with pairs)
OG_cluster_probs<-cooccur(mat = OG_cluster_data, spp_names = T, thresh = F, prob = "comb")
proc.time() # User time 64214 seconds
# Save image
save.image(file=file.path(newdir,"cooccurr_image"))
# Plot pairwise co-occurrences vs expected of OGs
plot(OG_cluster_probs)

# p_lt is the probability that the pair is not less frequently co-occurring than expected
# p_gt is the probability that the pair is not more frequently co-occurring than expected
prob.table(mod = OG_cluster_probs)

# Correct for multiple testing
OG_cluster_probs$results$gt_qval<-fdrtool(OG_cluster_probs$results$p_gt, statistic = "pvalue")$qval
# # Version without correction
# OG_cluster_probs$results$gt_qval<-OG_cluster_probs$results$p_gt

# Tabulate by occurrencies
# Filter by qval<0.025, expected ~6 false positives
cooc_freq<-OG_cluster_probs$results[which(OG_cluster_probs$results$gt_qval<0.025),c("sp1_name","sp2_name","gt_qval")]
cooc_freq<-as.data.frame(table(c(as.character(cooc_freq$sp1_name),as.character(cooc_freq$sp2_name))))
cooc_freq<-cooc_freq[order(cooc_freq$Freq, decreasing = T),]
names(cooc_freq)<-c("ODB8_OG", "Frequency")
ggplot(cooc_freq, aes(x=Frequency))+geom_histogram()+theme_bw()+scale_x_log10()+ggtitle(label = "Co-Duplications per\n Orthologous Group")+ylab(label = "Co-Duplicated OGs")+xlab(label = "Number of Co-Duplications")
# Top co-duplicating OGs
cooc_freq[1:20,]
# Total co-duplicating OGs
nrow(cooc_freq) #94
table(cooc_freq$Frequency)

# # Manually remove genes annotated with transposase domains
# transposons<-c(
#   "EOG8CZDZX", # also contains homeodomain-like
#   "EOG8WWV1R", # also contains zinc-finger
#   "EOG8DV81S",
#   "EOG83XXMT",
#   "EOG8D2962",
#   "EOG8CNT6C",
#   "EOG8FJBS3",
#   "EOG8CG2R8", # also contains peptidase
#   "EOG8VDSFF", # also includes domains for cell division
#   "EOG893614",
#   "EOG8WDGTK", # viral capsid
#   "EOG8HTC74",
#   "EOG8Z91BB",
#   "EOG8VX4M0", # also contains peptidase domains
#   "EOG8B2WDC",
#   "EOG8QNQD6", # also contains peptidase
#   "EOG80ZTCF", # DNA replicase with bacterial domains
#   "EOG8935ZZ", # also contains peptidase domains
#   "EOG84XN1N", 
#   "EOG88GZJM", # also contains chromatine restructuring domains
#   "EOG8M3CQ4",
#   "EOG8MGVNF", # also contains zinc-finger domains
#   "EOG8QJV2S",
#   "EOG8Q5C4B",
#   "EOG834ZP6", # also contains zinc-finger domains
#   "EOG8DRCV3", # also contains zinc-finger domains
#   "EOG8RV579", # viral packaging
#   "EOG8ZGRV6", # also contains zinc-finger domains
#   "EOG812PMP",
#   "EOG841SSS", # also contains nuclease domains
#   "EOG896286",
#   "EOG8B8MV0",
#   "EOG8D2968", # also contains homeodomain-like domains
#   "EOG8FXTNG", 
#   "EOG8J6V5R",
#   "EOG8PVRGG", # viral retrotranscriptase
#   "EOG8X9B0N", # transposase derived nuclease
#   "EOG8ZCWMQ", # helitron helicase transposon
#   "EOG8NZX77",
#   "EOG8MSGGG"
#   )

# Extract OG pairs and annotate in data.frame as OG1/OG2/PairQval
pairdata<-OG_cluster_probs$results[which(OG_cluster_probs$results$gt_qval<0.025),c("sp1_name","sp2_name","gt_qval")]
# # Remove pairs that include interactions with transposon OGs
# pairdata<-pairdata[-which(pairdata$sp1_name%in%transposons),]
# pairdata<-pairdata[-which(pairdata$sp2_name%in%transposons),]
pairdata<-droplevels(pairdata)
pairdata$pair_id<-paste(pairdata$sp1_name, pairdata$sp2_name, sep=":")
summary(pairdata)
dim(pairdata)
as.data.frame(table(c(as.character(pairdata$sp1_name), as.character(pairdata$sp2_name))))
hist(table(c(as.character(pairdata$sp1_name), as.character(pairdata$sp2_name))))

### Check via network graph
# Create network matrix 
networkgraph<-graph.edgelist(as.matrix(pairdata[,1:2]),directed = F)
tkplot(networkgraph)
pdf(file = file.path(graphdir, "coduplication_network.pdf"), paper = "a4r")
plot(networkgraph)
dev.off()
as_adjacency_matrix(networkgraph)
# lapply(graph.neighborhood(networkgraph, 1),plot) # plot first order neighbourhood for each point
cliques(networkgraph) # find cliques
# All ternary cliques include EOG805VG7 and EOG8C5G07, each minor node is only connected to hubs

##### Merge with OG to gene data
# OG/ninteractions/gene/cluster (on every individual gene)

clusterdata_2<-na.exclude(clusterdata_compiled[,c("geneID","clusterID","odb8_og_id")])

groups<-apply(pairdata, 1, function(x){
  OG1=x[1]
  OG2=x[2]
  clusters<-intersect(
    clusterdata_2[grep(OG1, clusterdata_2$odb8_og_id),"clusterID"],
    clusterdata_2[grep(OG2, clusterdata_2$odb8_og_id),"clusterID"])
  OGindexes<-union(which(clusterdata_2$odb8_og_id%in%OG1),which(clusterdata_2$odb8_og_id%in%OG2))
  clusterdata_2[intersect(which(clusterdata_2$clusterID%in%clusters), OGindexes),]
})
names(groups)<-paste(pairdata$sp1_name,pairdata$sp2_name, sep = ":")

#### Create list of pairID/cluster/OG1+OG2 (gene IDs)

# Function to bind colums together adding NAs if different lengths
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
# Create lists of data.frames for each pair, containing a list with one entry per cluster, containing a data.frame with genes from family a and b
grouped_genes<-lapply(groups, function(x){
  dlply(.data = x, .variables = .(clusterID), .fun = summarize, 
        # a=levels(droplevels(odb8_og_id))
        geneID_1=cbind.fill(geneID[which(odb8_og_id%in%levels(droplevels(odb8_og_id))[1])],
                          geneID[which(odb8_og_id%in%levels(droplevels(odb8_og_id))[2])])[,1],
        geneID_2=cbind.fill(geneID[which(odb8_og_id%in%levels(droplevels(odb8_og_id))[1])],
                            geneID[which(odb8_og_id%in%levels(droplevels(odb8_og_id))[2])])[,2]
  )
})

# Import GxM matrix
motifs_u5d2 <- read.delim("./Input/motif_target_sets/0.1_223_b1hmotifs_all_u1_mask.gXm.txt")
# Set score threshold (either >0.25 or less than 2000 genes, applied to each factor)
pval_treshold<-0.25
maxgenes_threshold<-2000
group_motifs<-motifs_u5d2
group_motifs[,-1]<- apply(group_motifs[,-1], 2, function(x){
  if (sum(x<pval_treshold, na.rm = T)>maxgenes_threshold) {
    x[which(rank(x, na.last = T, ties.method = "max")>maxgenes_threshold)]<-NA
  }  else {
    x[which(x>pval_treshold)]<-NA
  }
  x  }
)
# Subset to genes of interest
group_motifs<-group_motifs[which(group_motifs$X%in%ldply(groups)[,2]),]
# Convert to affinity scores abs(log10), higher is better
group_motifs[,-1]<-apply(group_motifs[,-1], c(1,2), function(x) abs(log10(x)))
# Rescale using z-transform on affinity scores for each motif
z<-function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}
group_motifs[,-1]<-apply(group_motifs[,-1], 2, z)
# Apply inverse logit (i.e. logistic) transformation, converts from -inf/0/+inf to 0/0.5/1 and robust towards outliers
group_motifs[,-1]<-apply(group_motifs[,-1], 2, inv.logit)

# Get average binding affinities per motif per OG in each pair
OG_affinities<-merge(clusterdata_2,group_motifs, by.x = "geneID", by.y = "X", all.y = T, all.x = F)
OG_affinities<-ddply(.data = OG_affinities, .variables = .(odb8_og_id,clusterID), .fun = numcolwise(mean), na.rm = F)
# Reshape as list of OG by cluster matrix (one per motif)
OG_affinities<-melt(data = OG_affinities, id.vars = c("odb8_og_id","clusterID"), variable_name = "motif")
# Replace NAs with 0
OG_affinities$value<-ifelse(is.na(OG_affinities$value), 0, OG_affinities$value)
# For every motif M in every cluster C calculate the score (Sadj)mc=Smc*(1-max(Sm-c)) 
# This weights the presence of a motif in an OG against the presence of the same motif in the same OG in different clusters (uniqueness score)
OG_affinities$index<-1:nrow(OG_affinities)
OG_adjaffinities<-vector(mode="numeric",length=nrow(OG_affinities))
OG_affinities<-dlply(.data = OG_affinities, .variables =.(motif))
lapply(X = OG_affinities, FUN = function(x) {for (i in 1:nrow(x)){
  OG_adjaffinities[x$index[i]]<<-x$value[i]*(1-max(x$value[-i], na.rm = T))
}})
OG_affinities <- ldply(.data = OG_affinities, .id = "motif")
OG_affinities <- cbind(OG_affinities[,c("motif","clusterID","odb8_og_id","value")],OG_adjaffinities)
rm(OG_adjaffinities)

# Convert to cluster by motif matrices split by OG
OG_affinities<-dlply(.data = OG_affinities[,c("motif","clusterID","odb8_og_id","OG_adjaffinities")], .variables = .(odb8_og_id))
OG_affinities<-lapply(X = OG_affinities, FUN = function(x){
  cast(x, formula = motif ~ clusterID)
})
OG_affinities
# Multiply matrices for matching pairs and plot scores
OG_pairaffinities<-vector(mode = "list", length = nrow(pairdata))
apply(pairdata, 1, function(x){
  OG_pairaffinities <- OG_affinities[[which(names(OG_affinities)==x[1])]]*OG_affinities[[which(names(OG_affinities)==x[2])]]
})
# Extract max for each motif for each OG pair for each cluster
# Filter out NAs and plot distribution


# # Get average binding affinities per motif per OG in each pair, then calculate difference between OGs in each group
# # Function generates NaNs if all genes in a mean have NA motif scores
# d_affinities<-lapply(X = grouped_genes, FUN = function(x){
#   lapply(X = x, FUN = function(y){
#     OG1<-apply(group_motifs[which(group_motifs[,1]%in%y$geneID_1),-1], MARGIN =  2, mean, na.rm = T)# create a vector of average affinities from the genes of OG1
#     OG2<-apply(group_motifs[which(group_motifs[,1]%in%y$geneID_2),-1], MARGIN =  2, mean, na.rm = T)# create a vector of average affinities from the genes of OG2
#     OG1*OG2 # calculate similarity in affinity, weighted by absolute expression
#   })
# })
# 
# d_affinities<-lapply(X = d_affinities, FUN = ldply, .id = "clusterID")
# d_affinities<-ldply(d_affinities, .id = "OG_pair")
# 
# # Replace NaNs and NAs with zero (absence of motif)
# d_affinities[,-c(1:2)]<-apply(X = d_affinities[,-c(1:2)], c(1,2), function(x){ifelse(is.na(x),0,x)})

# # Remove motifs that score zero across all pairs
# present<-ddply(.data = d_affinities, .variables = .(OG_pair), sapply, FUN=function(x){any(x!=0)})
# d_affinities<-droplevels(d_affinities[,which(apply(present, 2, sum)>1)])

# Compare to empirical null distribution across all clusters (fdr or quantiles)
d_affinities_2<-gather(data = d_affinities, key = motif, value = similarity_score, 3:ncol(d_affinities))
# ggplot(d_affinities, aes(x=similarity_score, group=OG_pair))+geom_density()+theme_bw()+facet_wrap( ~ motif)
ggplot(d_affinities_2[which(d_affinities_2$OG_pair=="EOG88D38W:EOG8C5G07"),], aes(y=similarity_score, x=motif))+geom_boxplot()+theme_bw()
ggplot(d_affinities_2, aes(x=similarity_score))+geom_density()+theme_bw()+scale_x_log10()
ggplot(d_affinities_2[which(d_affinities_2$similarity_score!=0),], aes(x=similarity_score))+geom_density()+theme_bw()

# How many motifs per cluster per pair do we observe at different thresholds?
table(ddply(.data = d_affinities_2, .variables = .(OG_pair, clusterID), .fun = summarize, nScoring=sum(similarity_score>0.292))[,3])
quantile(d_affinities_2$similarity_score[which(d_affinities_2$similarity_score>0)]) # using 3rd quartile of non-zero values as cutoff

# Retrieve only motifs observed once per cluster in each pair of co-duplicating paralogs 
d_unique <- d_affinities_2
d_unique$similarity_score <- d_unique$similarity_score>0.29 #2011 motifs overall
d_unique<-spread(data = d_unique, key = clusterID, value = similarity_score)
# retrieve and store first cluster with a non-zero score for that motif (we later select only motifs present in one cluster)
d_unique$clusterID<-names(d_unique[,-c(1:2)])[apply(X = d_unique[,-c(1:2)], MARGIN = 1, FUN = function(x){which(x>0)[1]})]
unique_motifs<-droplevels(d_unique[which(apply(d_unique[,-c(1:2,ncol(d_unique))], 1, function(x){sum(x, na.rm = T)==1})),c(1:2,ncol(d_unique))])

unique_motifs<- unique_motifs[order(unique_motifs$OG_pair,unique_motifs$clusterID),c(1,3,2)]
unique_motifs$clusterID<-as.factor(unique_motifs$clusterID)

# Concatenate factors in each module on a row and save as csv
write.csv(ddply(unique_motifs, .variables = .(OG_pair, clusterID), .fun = summarize, significant_motifs=paste(motif, collapse=", ")), file = file.path(newdir, "unique.motifs.csv"))

# Check cluster-specific motifs
unique_motifs [order(unique_motifs$clusterID, unique_motifs$OG_pair, unique_motifs$motif),]

# Break down by individual OGs (checks for gains of individual novel motifs)
OG_motifs<-unique_motifs
OG_motifs$OG1<- str_extract(string = OG_motifs$OG_pair, "^[[:alnum:]]*")
OG_motifs$OG2<- str_extract(string = OG_motifs$OG_pair, "[[:alnum:]]*$")
OG_m1<- OG_motifs[,-4]
names(OG_m1)[4]<-"odb8_og_id"
OG_m2<- OG_motifs[,-5]
names(OG_m2)[4]<-"odb8_og_id"
OG_motifs<-rbind(OG_m1, OG_m2)
rm(list = c("OG_m1", "OG_m2"))

table(table(OG_motifs$motif, OG_motifs$clusterID))

# Test number of motif by number of paralogs (paralogs as a proxy of cluster age, testing for accumulation of motifs over time)

############### Outdated code snippets
# Plot observed vs expected co-occurrences of OG pairs
# obs.v.exp_2<-function (mod) 
# {
#   ptab <- mod$results
#   ptab$signs <- ifelse(ptab$p_gt >= 0.05, 0, 1) + ifelse(ptab$p_lt >= 
#                                                            0.05, 0, -1)
#   exp_cooccur <- ptab$exp_cooccur
#   obs_cooccur <- ptab$obs_cooccur
#   signs <- ptab$signs
#   p <- ggplot(ptab, aes(x = exp_cooccur, y = obs_cooccur)) + 
#     geom_point(aes(fill = factor(signs, levels = c(-1, 0, 
#                                                    1))), colour = "black", pch = 21, size = 4, position="jitter", alpha="0.3")
#   p <- p + scale_fill_manual(values = c("#FFCC66", "dark gray", 
#                                         "light blue"), name = "", labels = c("negative", "random", 
#                                                                              "positive"), drop = FALSE)
#   p <- p + theme(plot.title = element_text(vjust = 2, size = 20, 
#                                            face = "bold"), legend.text = element_text(size = 18), 
#                  axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
#                  axis.text.x = element_text(hjust = 0, vjust = 1)) + xlab("Expected Co-occurrences") + 
#     ylab("Observed Co-occurrences")
#   p <- p + ggtitle("Observed-Expected Plot") + geom_abline(color = "green")
#   p
# }
# obs.v.exp_2(OG_cluster_probs)
# # Visualize as network
# # Create proximity matrix (og X og) containing effect sizes of paralogs in same cluster
# # Modify standard function (fails because subscripts in the original data file are more than those in the matrix)
# effect.sizes_2<-function (mod, standardized = TRUE, matrix = FALSE) 
# {
#   ptab <- mod$results
#   if ("sp1_name" %in% colnames(ptab)) {
#     sp1 <- as.character(ptab$sp1_name)
#     sp2 <- as.character(ptab$sp2_name)
#   }
#   else {
#     sp1 <- ptab$sp1
#     sp2 <- ptab$sp2
#   }
#   if (standardized == T) {
#     effs <- data.frame(sp1, sp2, effect = as.numeric(((ptab$obs_cooccur - 
#                                                          ptab$exp_cooccur)/mod$sites)))
#   }
#   else {
#     effs <- data.frame(sp1, sp2, effect = as.numeric((ptab$obs_cooccur - 
#                                                         ptab$exp_cooccur)))
#   }
#   if (matrix == F) {
#     effs
#   }
#   else {
#     effs_1 <- effs
#     effs_2 <- effs[, c(2, 1, 3)]
#     colnames(effs_2) <- c("sp1", "sp2", "effect")
#     effs <- rbind(effs_1, effs_2)
#     effs <- recast(data = effs, formula = sp1 ~ sp2, measure.var = "effect", 
#                    id.var = c("sp1", "sp2"))
#     row.names(effs) <- effs$sp1
#     effs$sp1 <- NULL
#     effs <- as.matrix(effs)
#     effs <- effs[, match(row.names(effs), colnames(effs))]
#     as.matrix(effs)
#   }
# }
# # Generate distance matrix
# distmat<-effect.sizes_2(OG_cluster_probs, standardized = T, matrix = T)
# # Bring all values to positives
# distmat<-distmat+abs(min(distmat, na.rm = T))
# # Remove NaNs
# 
# # Convert to distance
# dist<-as.dist(m = distmat)
# 
# # Apply hierarchical clustering
# clust<-hclust(d = dist)
# 
# # Plot as dendrogram
# 
# # What if we have a proximity matrix containing the p-vaues for positive co-occurrence instead?
# # Pros, network representation of significant co-associations, weighted undirected? network
# # Cons? Reflects probabilities, not actual co-occurrencies
# # or effect sizes from function effect.sizes (possibly restricted by significance?)
# Once per every group (lapply)
# Once per every cluster & Once per every motif (ddply)

# # Create list of pairs of genes to match on each column (one pair per group per cluster)
# lapply(groups[1:3], FUN = function(x){
#   ddply(.data = x, .variables = .(clusterID, odb8_og_id), .fun = summarize, geneID=unique(geneID))
# })

# #### Check for supergroups
# 
# ## Check for three-way interactions
# # given the groups a:b and b:c they form a supergroup if there is a significant co-duplication of a:c
# # check only for pairs of pairs with 1 OG in common
# otherpairs<-apply(pairdata, 1, function(x){
#   otherpairs<-union(which(pairdata$sp1_name%in%x[1]),which(pairdata$sp1_name%in%x[2]))
# })
# # Add pair names
# names(otherpairs)<-paste(pairdata$sp1_name,pairdata$sp2_name, sep = ":")
# # # Remove self-matches
# # for (i in 1:length(otherpairs)) function(x){
# #   otherpairs[[x]][-which(otherpairs[[x]]==as.character(x))]
# # }
# 
# patternlist<-as.list(1:length(otherpairs))
# names(patternlist)<-names(otherpairs)
# for (i in 1:length(otherpairs)) {
#   patternlist[i]<-lapply(otherpairs[[i]],function(x){
#     setdiff(union(unlist(pairdata[i,1:2]),unlist(pairdata[x,1:2])),intersect(unlist(pairdata[x,1:2]),unlist(pairdata[i,1:2])))
#   })
# }
# patternlist
# 
# # convert to data.frame
# patternlist<-ldply(patternlist[which(sapply(patternlist, function(x) length(x)>0 ))])
# 
# # Grep list of patterns on validated pairs (remember, list of patterns is a:c, a:b and b:c are already validated)
# patternlist$confirmed<-apply(patternlist, 1, function(x){
#   any(grepl(paste(x[2],x[3], sep=":"), pairdata$pair_id)|grepl(paste(x[3],x[2], sep=":"), pairdata$pair_id))
# })
