#### Search for linkage groups enriched in DE genes or clusters
# Initialize script
library(fdrtool)
library(plyr)
newdir <- file.path(getwd(),"./output/LinkageClusterEnrichment")
dir.create(newdir)

# load datasets
transcriptdata <- read.csv(file = "./Output/Results_compiler_CCRE/transcriptdata.csv")[,-1]

# Reduce to gene by cluster table (annotate each geneID once per cluster)
linkage <- transcriptdata[,c("clusterID","geneID","LinkageCluster","recombrate","Scaffold", "devsexbias")]
linkage <- linkage[which(duplicated(paste(linkage$clusterID,linkage$geneID))==F),]
# Remove genes with no linkage cluster assignment
linkage <- linkage[which(linkage$LinkageCluster!="None"),]
linkage <- droplevels(linkage[order(linkage$clusterID, linkage$geneID),])
ddply(.data = linkage, .variables = .(clusterID), .fun = summarize, allOK = anyDuplicated(geneID))


### Apply fisher's exact test comparing male and female biased genes in cluster vs outside cluster
linkage$sexbias <- factor(ifelse(grepl(pattern = "m", linkage$devsexbias), 
                                 ifelse(grepl(pattern = "f", x = linkage$devsexbias), "mf", "m"), 
                                 ifelse(grepl(pattern = "f", x = linkage$devsexbias), "f", "0")
))

# Version with four tests, m, f, mf and no bias
# tests <- expand.grid(levels(linkage$sexbias),levels(linkage$LinkageCluster))
# linkag2_p.values <- apply(X = tests, MARGIN = 1, function(x){
#   fisher.test(x = linkage$sexbias==x[1],y = linkage$LinkageCluster==x[2], alternative = "greater")$p.value
# })

# Version with two tests, m and f (mf are counted in both)
# Alternative version, using DE annotation from cluster
# tests <- expand.grid(c("m","f"),levels(linkage$LinkageCluster))
# linkag2_p.values <- apply(X = tests, MARGIN = 1, function(x){
#   fisher.test(x = grepl(pattern = x[1], x = linkage$sexbias), 
#               y = linkage$LinkageCluster==x[2], 
#               alternative = "greater")$p.value
# })

# Using DE annotation from individual nodes, collapsed onto genes 
genesummary <- read.csv(file = "./Output/Results_compiler_CCRE/genesummary.csv")[,-1]
LinkageGroups <- unique(transcriptdata[,c("geneID","LinkageCluster")])
linkag2 <- merge(genesummary, LinkageGroups, by = "geneID", all.x = T, all.y = F)
linkag2$inclusiveMale <- linkag2$malebiased|linkag2$bothbiased
linkag2$inclusiveFemale <- linkag2$femalebiased|linkag2$bothbiased
tests <- levels(droplevels(linkag2$LinkageCluster))

linkag2_p.values_m <- sapply(X = tests, function(t){
  fisher.test(x = linkag2$malebiased, 
              y = linkag2$LinkageCluster==t, 
              alternative = "greater")$p.value
})
linkag2_p.values_m <- data.frame(LinkageCluster = tests, pval = linkag2_p.values_m)
linkag2_p.values_m$Sex <- "Male"

linkag2_p.values_f <- sapply(X = tests, function(t){
  fisher.test(x = linkag2$femalebiased, 
              y = linkag2$LinkageCluster==t, 
              alternative = "greater")$p.value
})
linkag2_p.values_f <- data.frame(LinkageCluster = tests, pval = linkag2_p.values_f)
linkag2_p.values_f$Sex <- "Female"

linkag2_p.values <- rbind(linkag2_p.values_f, linkag2_p.values_m)
# Filter via fdr
linkag2_fdr <- fdrtool(x = linkag2_p.values$pval, statistic = "pvalue")$lfdr
linkag2_fdr <- cbind(linkag2_p.values, linkag2_fdr)
names(linkag2_fdr) <- c("LinkageCluster", "pval", "sexbias", "fdr")
write.csv(x = linkag2_fdr, file = file.path(newdir, "Node_linkage_enrichment_fdr.csv"))

# Filter only significant linkage groups
table(linkag2_fdr$sexbias, linkag2_fdr$fdr<.05)
linkag2_names <- linkag2_fdr[which(linkag2_fdr$fdr<0.05),c("sexbias", "LinkageCluster","fdr")]
# modify for paper: add linkage rates, number of male, female and total genes, and medians across all (non-chromosome) linkage clusters
allcounts <-  ddply(.data = linkage, .variables = .(LinkageCluster), .fun = summarize, 
                    recombrate = recombrate[1],
                    maleproportion = sum(grepl(pattern = "m", x = devsexbias))/length(devsexbias),
                    femaleproportion = sum(grepl(pattern = "f", x = devsexbias))/length(devsexbias),
                    totalgenes = length(devsexbias)
)
linkag2_names <- merge(x = linkag2_names, y = allcounts, by = "LinkageCluster", all.x = T, all.y = F)

medians <- c("Median Values", NA, NA,
             median(allcounts$recombrate), 
             median(allcounts$maleproportion), 
             median(allcounts$femaleproportion), 
             median(allcounts$totalgenes))

linkag2_names <- rbind(medians, linkag2_names)
write.csv(x = linkag2_names, file = file.path(newdir, "Node_linkage_enrichment_0.05_names.csv"))

# How many sex-biased genes do all m and f clusters contain, and how many are there in the genome?
sexcounts <-  ddply(.data = linkage, .variables = .(LinkageCluster), .fun = summarize, 
                    recombrate = recombrate[1],
                    malecounts = sum(grepl(pattern = "m", x = devsexbias)),
                    femalecounts = sum(grepl(pattern = "f", x = devsexbias)),
                    totalgenes = length(devsexbias)
)
apply(X = sexcounts[,3:5],MARGIN = 2, sum)
linkag2_names2 <- linkag2_fdr[which(linkag2_fdr$fdr<0.05),c("sexbias", "LinkageCluster","fdr")]
sexcounts <- merge(x = linkag2_names2, y = sexcounts, by = "LinkageCluster", all.x = T, all.y = F)
ddply(.data = sexcounts, .variables = .(sexbias), .fun = summarize, 
      male = sum(malecounts),
      female = sum(femalecounts)
)
ddply(.data = sexcounts, .variables = .(LinkageCluster), .fun = summarize, 
      male = sum(malecounts),
      female = sum(femalecounts),
      total = sum(totalgenes)
)

# Write gene lists for each linkage cluster with confirmed enrichment in sexDE genes
sapply(X = as.character(linkag2_names$LinkageCluster), FUN = function(x){
  write(x = as.character(unique(transcriptdata[which(transcriptdata$LinkageCluster==x), "geneID"])), 
        file = file.path(newdir, paste("LinkageCluster",x,".txt", sep = ""))
  )
})

# Do enriched LinkageCluster show reduced recombination rates compared to others?
# Plot distribution and annotate distance in SD
linkagedata <- na.exclude(unique(transcriptdata[,c("LinkageCluster", "recombrate")]))
linkagedata <- linkagedata[-which(linkagedata$LinkageCluster%in%c(1:5)),]
library(lattice)
densityplot(log10(linkagedata$recombrate))
linkagedata$logNormRecombrate <- log10(linkagedata$recombrate)
linkagedata$logNormRecombrate <- linkagedata$logNormRecombrate-mean(linkagedata$logNormRecombrate)
linkagedata$logNormRecombrate <- linkagedata$logNormRecombrate/sd(linkagedata$logNormRecombrate)
densityplot(linkagedata$logNormRecombrate)

# Extract SD of Linkage Clusters of interest (also quantiles) 
linkagedata[which(linkagedata$LinkageCluster%in%linkag2_names$LinkageCluster),]$logNormRecombrate
fivenum(linkagedata$recombrate)
# merge allcounts with sex bias index
allcounts$sexbias <- linkag2_names[match(x = allcounts$LinkageCluster, table = linkag2_names$LinkageCluster),"sexbias"]

ggplot(data = allcounts, mapping = aes(x = is.na(sexbias), y = recombrate)) + geom_boxplot(notch = T) + scale_y_log10()

wilcox.test(formula = recombrate ~ is.na(sexbias), data = allcounts, conf.int = T)

sexlinkage_rates <- ddply(.data = SexLinkage, .variables = .(LinkageCluster), summarize, 
                          recombrate = recombrate[1],
                          sexbias = sexbias[1]
)

sexlinkage_rates[which(sexlinkage_rates$sexbias!="Baseline"),]
ggplot(data = sexlinkage_rates, mapping = aes(x = sexbias!="Baseline", y = recombrate)) + geom_boxplot(notch = T, varwidth = T) + scale_y_log10() + theme_bw() + geom_point()


# # Apply fisher's exact test for each cluster per linkage cluster vs global (one sided, enrichment)
# # generate list of clusters and linkage clusters
# tests <- expand.grid(levels(droplevels(linkage$clusterID)), levels(droplevels(linkage$LinkageCluster)))
# names(tests) <- c("CorClust","LinkClust")
# 
# # tabulate for matching vs non-matching and apply fisher.test (mapply)
# linkage_p.values <- apply(X = tests, MARGIN = 1, function(x){
#   fisher.test(x = linkage$clusterID==x[1],y = linkage$LinkageCluster==x[2], alternative = "greater")$p.value
# })
# 
# # Filter p-values via fdr
# linkage_fdr <- fdrtool(x = linkage_p.values, statistic = "pvalue")
# linkage_fdr <- cbind(tests, linkage_fdr$lfdr)
# write.csv(x = linkage_fdr, file = file.path(newdir, "Cluster_linkage_enrichment_fdr.csv"))