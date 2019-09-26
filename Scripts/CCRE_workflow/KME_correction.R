## Calculate KMEs and correct for their effect on cluster parameters
## Initialize script
date()
rm(list=ls())
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(WGCNA)
allowWGCNAThreads()

relativescale <- function(x){
  y = x+abs(min(x, na.rm = T))
  y/max(y, na.rm = T)
}


## initialize output path
newdir<-file.path(getwd(), "Output/KME_correction")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/KME_correction")
dir.create(graphdir)

## Load data
load(file = "./Output/WGCNA_CCRE_clustering_biweight/CCRE_MEList.RData")
eigenexons <- MEList$eigengenes
clusterassignments <- read.csv(file = "./Output/WGCNA_CCRE_clustering_biweight/CCRE_clusterassignments_simple.csv")[,-1]
CCRE_nasdevgeneexon <- read.csv(file = "./Output/CCREs/CCRE_nasdevgeneexon.csv")
row.names(CCRE_nasdevgeneexon) <- CCRE_nasdevgeneexon$X
CCRE_nasdevgeneexon <- CCRE_nasdevgeneexon[,-1]
# CCRE_nasdevgeneexon <- CCRE_nasdevgeneexon[sample(1:nrow(CCRE_nasdevgeneexon),1000),]
t_CCRE_nasdevgeneexon <- t(CCRE_nasdevgeneexon)

## Control for avg clustering coef and MAR to eigenexon proximity (or module KME)
# Calculate KMEs
KMEmat <- signedKME(datExpr = t_CCRE_nasdevgeneexon, datME = eigenexons, outputColumnName = "",
                    corFnc = "bicor", corOptions = "use = 'p'")
# Select only KME of each node to its respective module
KME <- mapply(function(x,y){KMEmat[x,y]}, 
              x = match(clusterassignments$transcriptID, row.names(KMEmat)),
              y = match(clusterassignments$clusterID, colnames(KMEmat))
)
KME <- data.frame(nodeID = clusterassignments$transcriptID, 
                  clusterID = clusterassignments$clusterID,
                  KME = KME)

# Check for accuracy
testclusterID <- "yellowgreen"
KMEmat[which(row.names(KMEmat)%in%clusterassignments[which(clusterassignments$clusterID==testclusterID),"transcriptID"]),testclusterID]
KME[which(KME$clusterID==testclusterID),"KME"]
KMEmat[which(row.names(KMEmat)%in%clusterassignments[which(clusterassignments$clusterID==testclusterID),"transcriptID"]),testclusterID]-KME[which(KME$clusterID==testclusterID),"KME"]
match(KMEmat[which(row.names(KMEmat)%in%clusterassignments[which(clusterassignments$clusterID==testclusterID),"transcriptID"]),testclusterID], KME[which(KME$clusterID==testclusterID),"KME"])

# Reorder and save
KME <- KME[order(KME$clusterID, KME$KME, KME$nodeID, decreasing = T),]
write.csv(x = KME, file = file.path(newdir, "KME.csv"))


## Load and merge other network parameters
WithinClusterNetworkConcepts_nodes <- read.csv(file = "./Output/NetworkConcepts_CCRE/CCRE_nodeparams.csv")[,-1]
names(WithinClusterNetworkConcepts_nodes)[2] <-  "nodeID"
WithinClusterBetweenness <- read.csv(file = "./Output/w_betweenness_CCREs/module_emp_betweenness.csv")[,-1]

WithinClusterNetworkConcepts <- merge(WithinClusterNetworkConcepts_nodes, WithinClusterBetweenness, by = c("clusterID","nodeID"))
WithinClusterNetworkConcepts <- merge(WithinClusterNetworkConcepts, KME, by = c("clusterID","nodeID"))

## Rescale Connectivity and Betweenness to max possible (densities, or absolute params)
AbsParams <- ddply(.data = WithinClusterNetworkConcepts, .variables = .(clusterID), .fun = summarise, 
                   nodeID = nodeID,
                   AbsConnectivity = Connectivity/length(Connectivity-1),
                   AbsBetweenness = Betweenness/((length(Betweenness)-1)*(length(Betweenness)-2)/2),
                   MaxBetweenness = rep((length(Betweenness)-1)*(length(Betweenness)-2)/2,length(Betweenness)),
                   MaxConnectivity = rep(length(Connectivity-1),length(Connectivity))
)

## ReScale all network parameters to max observed within cluster (relative parameters)
RelParams <- ddply(.data = WithinClusterNetworkConcepts, .variables = .(clusterID), .fun = summarise, 
                   nodeID = nodeID,
                   RelKME = relativescale(KME),
                   RelConnectivity = ScaledConnectivity,
                   # RelC2  = Connectivity/max(Connectivity), # Debugging line, checks for consistency with pre-rescaled connectivity via WGCNA
                   RelBetweenness = Betweenness/max(Betweenness),
                   RelClusterCoef = ClusterCoef/max(ClusterCoef),
                   RelMAR = MAR/max(MAR)
)

# Merge
allParams <- merge(AbsParams,   RelParams, by = c("clusterID","nodeID"))
allParams <- merge(allParams, WithinClusterNetworkConcepts[,c("clusterID","nodeID","Connectivity","Betweenness","ClusterCoef","MAR","KME")])
allParams <- allParams[,c("clusterID","nodeID","KME", "Connectivity","Betweenness","AbsConnectivity","AbsBetweenness","ClusterCoef","MAR","RelKME","RelConnectivity","RelBetweenness","RelClusterCoef","RelMAR","MaxBetweenness","MaxConnectivity")]

# Add hubness
Hubness <- ddply(.data = allParams, .variables = .(clusterID), .fun = summarize, 
                 nodeID = nodeID, 
                 Hubness = AbsConnectivity*(1-ClusterCoef))
Hubness <- ddply(.data = allParams, .variables = .(clusterID), .fun = summarize, 
                 nodeID = nodeID, 
                 Hubness = AbsConnectivity*(1-ClusterCoef),
                 RelHubness = Hubness/max(Hubness))

allParams <- merge(allParams, Hubness, by = c("nodeID","clusterID"))
# Save as csv
write.csv(x = allParams, file = file.path(newdir, "WithinCluster_NetworkParams_rescaled.csv"))

## Plot to detect CORE
ggplot(data = allParams, mapping = aes(y = AbsConnectivity, x = KME, col = clusterID=="grey"))+geom_smooth()+scale_y_log10()+scale_x_continuous(breaks = c(1:10)*0.1) + geom_vline(xintercept = 0.6)
ggplot(data = allParams, mapping = aes(y = ClusterCoef, x = KME, col = clusterID=="grey"))+geom_smooth()+scale_y_log10()+scale_x_continuous(breaks = c(1:10)*0.1) + geom_vline(xintercept = 0.6)
ggplot(data = allParams, mapping = aes(y = MAR, x = KME, col = clusterID=="grey"))+geom_smooth()+scale_y_log10()+scale_x_continuous(breaks = c(1:10)*0.1) + geom_vline(xintercept = 0.6)
ggplot(data = allParams, mapping = aes(y = AbsBetweenness+0.01, x = KME, col = clusterID=="grey"))+geom_smooth()+scale_y_log10()+scale_x_continuous(breaks = c(1:10)*0.1) + geom_vline(xintercept = 0.6)
ggplot(data = allParams, mapping = aes(y = AbsBetweenness+0.01, x = KME+1, col = clusterID=="grey"))+geom_smooth(se = F)+scale_y_log10()+scale_x_continuous(breaks = c(1:10)*0.1+1) + geom_vline(xintercept = 1.6)

## Add CORE annotation
allParams$CORE <- as.factor(ifelse(allParams$KME>=0.65&allParams$clusterID!="grey","CORE","periphery"))

table(allParams$clusterID, allParams$CORE)
prop.table(table(allParams$clusterID, allParams$CORE), margin = 1)

CORE <- allParams[which(allParams$CORE=="CORE"),-which(names(allParams)=="CORE")]

## Apply GLMs to correct for KME correlation in core nodes
# library(robust)
gammalog1 <- glm(formula = AbsConnectivity ~ KME, family = Gamma(link="log"), data = CORE)
gammainv1 <- glm(formula = AbsConnectivity ~ KME, family = Gamma, data = CORE)
gammalogit1 <- glm(formula = AbsConnectivity+0.1 ~ KME, family = Gamma(link = "logit"), data = CORE)
# gammarobinv1 <- glmrob(formula = AbsConnectivity ~ KME, family = Gamma, data = CORE)
# gammaroblog1 <- glmrob(formula = AbsConnectivity ~ KME, family = Gamma(link = "log"), data = CORE)
gauslog1 <- glm(formula = AbsConnectivity ~ KME, family = gaussian(link = "log"), data = CORE)
gauslogit1 <- glm(formula = AbsConnectivity+0.1 ~ KME, family = gaussian(link = "logit"), data = CORE)
gausid1 <- glm(formula = AbsConnectivity ~ KME, family = gaussian(link = "identity"), data = CORE)
pois1 <- glm(formula = AbsConnectivity ~ KME, family = quasipoisson(link = "log"), data = CORE)
binomlogit1 <- glm(formula = AbsConnectivity ~ KME, family = binomial(link = "logit"), data = CORE)

AIC(gammalog1, gammainv1, gauslog1, gammalogit1, gauslogit1, gausid1, pois1)

# normplots <- melt(data.frame(gammalog1 = resid(gammalog1),
#                              gammaroblog1 = resid(gammaroblog1),
#                              gauslog1 = resid(gauslog1),
#                              gammalogit1 = resid(gammalogit1),
#                              simulres_sd.3 =  rnorm(n = 15836, mean = 0, sd = 0.3),
#                              simulres_sd.1 =  rnorm(n = 15836, mean = 0, sd = 0.065))
# )
# ggplot(data = normplots, aes(x = value, col = variable))+geom_density()+theme_bw()

### Log link function models perform well, with a gamma distribution of residuals
gammalog2 <- glm(formula = AbsConnectivity ~ KME+RelClusterCoef, family = Gamma(link="log"), data = CORE)
gammalog2 <- glm(formula = AbsConnectivity ~ KME+RelClusterCoef+RelMAR, family = Gamma(link="log"), data = CORE)

### Apply gaus log model to the three main node parameters
K_Conn <- glm(formula = AbsConnectivity ~ KME, family = Gamma(link="log"), data = CORE)
qqnorm(resid(K_Conn))
qqline(resid(K_Conn))

K_Clus <- glm(formula = ClusterCoef ~ KME, family = Gamma(link="log"), data = CORE)
qqnorm(resid(K_Clus))
qqline(resid(K_Clus))

K_MAR <- glm(formula = MAR ~ KME, family = Gamma(link="log"), data = CORE)
qqnorm(resid(K_MAR))
qqline(resid(K_MAR))

K_corr <- data.frame(
  clusterID = CORE$clusterID,
  nodeID = CORE$nodeID,
  K_Conn = resid(K_Conn),
  K_Clus = resid(K_Clus),
  K_MAR = resid(K_MAR)
)

pairs(K_corr[,-c(1:2)],
      c("AbsConnectivity", "ClusterCoef", "MAR"))
pairs(log10(K_corr[,-c(1:2)]+5),
      c("AbsConnectivity", "ClusterCoef", "MAR"))

cor.test(log(resid(K_Conn)+10),log(resid(K_Clus)+10))
cor.test(log(resid(K_Conn)+10),log(resid(K_MAR)+10))
cor.test(log(resid(K_Clus)+10),log(resid(K_MAR)+10))

### Test model fit for betweenness

binomlogit <- glm(formula = cbind(as.integer(Betweenness),as.integer(MaxBetweenness-Betweenness)) ~ KME*MAR, 
                  family = binomial(link="logit"), data = CORE)
quasibinom_null <- glm(formula = AbsBetweenness ~ 1, 
                      family = quasibinomial, data = CORE)   
quasibinom_KME <- glm(formula = AbsBetweenness ~ KME, 
                  family = quasibinomial, data = CORE)  
quasibinom_full <- glm(formula = AbsBetweenness ~ (KME+MAR+AbsConnectivity+ClusterCoef)^2, 
                       family = quasibinomial, data = CORE)  
quasibinom_best <- glm(formula = AbsBetweenness ~ KME+MAR+AbsConnectivity+ClusterCoef+KME:MAR+KME:AbsConnectivity+KME:ClusterCoef+MAR:AbsConnectivity, 
                       family = quasibinomial, data = CORE)

# quasibinom_best2 <- glm(formula = AbsBetweenness ~ (KME+MAR+AbsConnectivity)^2+ClusterCoef+ClusterCoef:MAR+ClusterCoef:KME+I(clusterID=="grey"), 
#                        family = quasibinomial, data = CORE)  

par(mfrow = c(2,2))
qqnorm(resid(quasibinom_null))
qqline(resid(quasibinom_null))
qqnorm(resid(quasibinom_KME))
qqline(resid(quasibinom_KME))
qqnorm(resid(quasibinom_full))
qqline(resid(quasibinom_full))
qqnorm(resid(quasibinom_best))
qqline(resid(quasibinom_best))
par(mfrow = c(1,1))

# Add to corrected dataset
K_corr <- cbind(K_corr, K_betw = residuals(quasibinom_KME))
# Rescale from 0 to 1 relative to highest and lowest in cluster CORE
K_corr_rescaled <- ddply(.data = K_corr, .variables = .(clusterID), .fun = summarize,
                         nodeID = nodeID,
                         K_Conn_rel = relativescale(K_Conn),
                         K_Clus_rel = relativescale(K_Clus),
                         K_MAR_rel = relativescale(K_MAR),
                         K_betw_rel = relativescale(K_betw)
)
# 18838 generates NA, clusterID saddlebrown, has only 1 core gene..
## Save as csv
write.csv(x = K_corr_rescaled, file = file.path(newdir, "KME_corrected_params.csv"))

### Compare hub score vs betweenness centrality distributions
ggplot(data = allParams, mapping = aes(x = AbsBetweenness+1E-6, y = Hubness)) + geom_point() + geom_smooth() + scale_x_log10() + facet_wrap(facets = ~clusterID)

ggplot(data = allParams, mapping = aes(x = RelBetweenness+1E-5, y = RelHubness)) + geom_point() + geom_smooth() + scale_x_log10()
cor.test(x = log(allParams$RelBetweenness+1e-5), allParams$RelHubness)

ggplot(data = allParams, mapping = aes(x = RelBetweenness>0, y = RelHubness)) + geom_boxplot() + theme_bw()
ggplot(data = allParams, mapping = aes(col = RelBetweenness>0, x = RelHubness)) + geom_density() + theme_bw()


fulldata <- read.csv(file = "./Output/Results_compiler_CCRE/nodedata_full.csv")[,-1]
pdf(file = file.path(graphdir, "betweenness_vs_hubness.pdf"))
ggplot(data = fulldata, mapping = aes(x = RelBetweenness_wc+1E-5, y = RelHubness_wc)) + geom_point() + geom_smooth() + scale_x_log10() + facet_wrap(~strata) + theme_bw()
ggplot(data = fulldata, mapping = aes(col = RelBetweenness_wc>0, x = RelHubness_wc)) + geom_density() + theme_bw()
dev.off()

# #########
# # Load compiled nodedata file
# nodedata <- read.csv("./Output/Results_compiler_CCRE/nodedata_full.csv")[,-1]
# # Check KME distribution
# ggplot(nodedata, aes(x = KME, col = node_type))+geom_density()+facet_grid(.~Grey)
# 
# # Plot relationship between eigenmodule and connectivity, split by type
# nodedata_graph <-  nodedata
# nodedata_graph$corsplit <- nodedata_graph$KME>0.5
# nodedata_graph$cortype <- paste(nodedata_graph$corsplit, nodedata_graph$node_type, sep = "")
# 
# ggplot(data = droplevels(nodedata_graph), mapping = aes(x = KME, y = ScaledConnectivity_WithinCluster, group = node_type)) + 
#   geom_point(alpha = 0.3)  + geom_smooth(method = "lm", mapping = aes(group = corsplit)) +
#   scale_y_log10() + facet_grid(node_type~Grey)
# 
# ggplot(data = droplevels(nodedata_graph), 
#        mapping = aes(x = KME, y = NormalizedBetweenness_WithinCluster+1, group = node_type)) + 
#   geom_point(alpha = 0.3)  + geom_smooth(method = "lm", mapping = aes(group = corsplit)) +
#   scale_y_log10() + facet_grid(node_type~Grey)
# 
# ## Compare within same graph
# pdf(file = file.path(graphdir, "EigenmoduleVsCentralities.pdf"), paper = "a4r")
# # Scaled Connectivity
# ggplot(data = droplevels(nodedata_graph), mapping = aes(x = KME, y = ScaledConnectivity_WithinCluster, group = node_type, col = node_type)) + 
#   geom_point(alpha = 0.3)  + geom_smooth(method = "lm", mapping = aes(group = cortype)) +
#   scale_y_log10() + facet_grid(.~Grey)
# # Scaled Betweenness
# ggplot(data = droplevels(nodedata_graph), mapping = aes(x = KME, y = ScaledBetweenness_WithinCluster, group = node_type, col = node_type)) + 
#   geom_point(alpha = 0.3)  + geom_smooth(method = "lm", mapping = aes(group = cortype)) +
#   scale_y_log10() + facet_grid(.~Grey)
# # MAR
# ggplot(data = droplevels(nodedata_graph), mapping = aes(x = KME, y = MAR_WithinCluster, group = node_type, col = node_type)) + 
#   geom_point(alpha = 0.3)  + geom_smooth(method = "lm", mapping = aes(group = cortype)) +
#   scale_y_log10() + facet_grid(.~Grey)
# # Clustering Coefficient
# ggplot(data = droplevels(nodedata_graph), mapping = aes(x = KME, y = ClusterCoef_WithinCluster, group = node_type, col = node_type)) + 
#   geom_point(alpha = 0.3)  + geom_smooth(method = "lm", mapping = aes(group = cortype)) +
#   scale_y_log10() + facet_grid(.~Grey)
# dev.off()
# 
# ### Preliminary observations:
# # There are two main groups of KME nodes in each cluster, with a splt at ~0.4 in correlation
# # Nodes below this threshold show independent (possibly decreasing) relationship of centrality measures to KME
# # Nodes above this threshold show strong positive dependence of all centrality measures to KME
# # CCREs have lower mean KME in the first class and higher in the second class
# # CCREs have higher intercept in both cases even after controlling for KME values
# 
# nodedata_graph <- nodedata_graph[,c("clusterID","nodeID", "node_type", "KME", "cortype", "ScaledConnectivity_WithinCluster", "MAR_WithinCluster", "ClusterCoef_WithinCluster", "ScaledBetweenness_WithinCluster", "Grey")]
# nodedata_graph <- melt(nodedata_graph, id.vars = c("clusterID","nodeID", "node_type", "KME", "cortype", "Grey"))
# 
# ggplot(data = nodedata_graph,
#        mapping = aes(x = KME, y = value, group = node_type, col = Grey)) + 
#   geom_point(alpha = 0.1)  + geom_smooth(method = "loess", mapping = aes(group = Grey)) +
#   scale_y_log10() + facet_grid(node_type~variable) + theme_bw() 
# 
# nodedata_graph$highcor <- factor(ifelse(nodedata_graph$KME>0.5, "high", "low"), levels = c("low","high"))
# 
# ggplot(data = nodedata_graph, mapping = aes(x = highcor, y = value+0.001, col = node_type)) + 
#   stat_boxplot(notch = T) + scale_y_log10() + facet_grid(variable~Grey) + theme_bw() 
# 
# ### Additional observations
# # Clusters can be divided bi-modally by KME between "core" and "subclusters"
# # The two groups show the same coarse distribution of parameters between different node types in clusters
# # CCREs are found only in the periphery of Grey
# 
# ## Plot dependencies between variables in cores vs subclusters
# nodedata$highcor <- factor(ifelse(nodedata$KME>0.5, "high", "low"), levels = c("low","high"))
# 
# ggplot(data = droplevels(nodedata[which(nodedata$Grey=="cluster"),]),
#        mapping = aes(x = ScaledConnectivity_WithinCluster-KME, y = ClusterCoef_WithinCluster, col = node_type)) + 
#   geom_point(alpha = 0.1)  + geom_smooth(method = "lm") +
#   theme_bw() + facet_grid(.~highcor)
# 
# ggplot(data = droplevels(nodedata[which(nodedata$Grey=="cluster"),]),
#        mapping = aes(x = ScaledConnectivity_WithinCluster-KME, y = MAR_WithinCluster, col = node_type)) + 
#   geom_point(alpha = 0.1)  + geom_smooth(method = "lm") +
#   theme_bw() + facet_grid(.~highcor)
# 
# ggplot(data = droplevels(nodedata[which(nodedata$Grey=="cluster"),]),
#        mapping = aes(x = ScaledConnectivity_WithinCluster-KME, y = ScaledBetweenness_WithinCluster+1E-6, col = node_type)) + 
#   geom_point(alpha = 0.1)  + geom_smooth(method = "lm") +
#   scale_y_log10() + theme_bw() + facet_grid(.~highcor)

### Outdated code snippets
# ## Exclude grey
# nodedata_graph_nogrey <- droplevels(nodedata_graph[which(nodedata_graph$clusterID!="grey"),])
# pdf(file = file.path(graphdir, "EigenmoduleVsCentralities_noGrey.pdf"), paper = "a4r")
# # Scaled Connectivity
# ggplot(data = droplevels(nodedata_graph_nogrey), mapping = aes(x = KME, y = ScaledConnectivity_WithinCluster, group = node_type, col = node_type)) + 
#   geom_point(alpha = 0.3)  + geom_smooth(method = "lm", mapping = aes(group = cortype)) +
#   scale_y_log10()
# # Scaled Betweenness
# ggplot(data = droplevels(nodedata_graph_nogrey), mapping = aes(x = KME, y = ScaledBetweenness_WithinCluster, group = node_type, col = node_type)) + 
#   geom_point(alpha = 0.3)  + geom_smooth(method = "lm", mapping = aes(group = cortype)) +
#   scale_y_log10()
# # MAR
# ggplot(data = droplevels(nodedata_graph_nogrey), mapping = aes(x = KME, y = MAR_WithinCluster, group = node_type, col = node_type)) + 
#   geom_point(alpha = 0.3)  + geom_smooth(method = "lm", mapping = aes(group = cortype)) +
#   scale_y_log10()
# # Clustering Coefficient
# ggplot(data = droplevels(nodedata_graph_nogrey), mapping = aes(x = KME, y = ClusterCoef_WithinCluster, group = node_type, col = node_type)) + 
#   geom_point(alpha = 0.3)  + geom_smooth(method = "lm", mapping = aes(group = cortype)) +
#   scale_y_log10()
# dev.off()
# 
# 
# ## Only grey
# nodedata_graph_grey <- droplevels(nodedata_graph[which(nodedata_graph$clusterID=="grey"),])
# nodedata_graph_grey <- nodedata_graph_grey[,c("clusterID","nodeID", "node_type", "KME", "ScaledConnectivity_WithinCluster", "MAR_WithinCluster", "ClusterCoef_WithinCluster", "ScaledBetweenness_WithinCluster")]
# nodedata_graph_grey <- melt(nodedata_graph_grey, id.vars = c("clusterID","nodeID", "node_type", "KME"))
# 
# pdf(file = file.path(graphdir, "grey_ClusterparamsVsKME.pdf"), paper = "a4r")
# 
# ggplot(data = droplevels(nodedata_graph_grey), mapping = aes(x = KME, y = value, group = node_type, col = node_type)) + 
#   geom_point(alpha = 0.3)  + geom_smooth(method = "loess", mapping = aes(group = node_type)) +
#   scale_y_log10() + facet_grid(node_type~variable) + theme_bw() 
# 
# dev.off()
# 
# ## Conforms to a random distribution of proportion (maximal variance toward the center, minimal at the extremes)
# 
# nodedata_graph_nogrey <- droplevels(nodedata_graph[which(nodedata_graph$clusterID!="grey"),])
# nodedata_graph_nogrey <- nodedata_graph_nogrey[,c("clusterID","nodeID", "node_type", "KME", "cortype", "ScaledConnectivity_WithinCluster", "MAR_WithinCluster", "ClusterCoef_WithinCluster", "ScaledBetweenness_WithinCluster")]
# nodedata_graph_nogrey <- melt(nodedata_graph_nogrey, id.vars = c("clusterID","nodeID", "node_type", "KME", "cortype"))
# nodedata_graph_nogrey <- droplevels(nodedata_graph_nogrey[-which(is.na(nodedata_graph_nogrey$node_type)),])
# 
# pdf(file = file.path(graphdir, "nogrey_ClusterparamsVsKME.pdf"), width = "")
# 
# ggplot(data = droplevels(nodedata_graph_nogrey), mapping = aes(x = KME, y = value, group = node_type, col = node_type)) + 
#   geom_point(alpha = 0.1)  + geom_smooth(method = "loess", mapping = aes(group = cortype)) +
#   scale_y_log10() + facet_grid(node_type~variable) + theme_bw() 
# 
# dev.off()
# 
# ggplot(data = droplevels(nodedata_graph_nogrey), mapping = aes(x = KME, y = value, group = node_type, col = node_type)) + 
#   geom_smooth(method = "loess", mapping = aes(group = cortype)) +
#   scale_y_log10() + facet_grid(.~variable) + theme_bw()
# 
# nodedata_graph_2 <- droplevels(nodedata_graph[,c("clusterID","nodeID", "node_type", "KME", "cortype", "ScaledConnectivity_WithinCluster", "MAR_WithinCluster", "ClusterCoef_WithinCluster", "ScaledBetweenness_WithinCluster")])
# nodedata_graph_2$grey <- as.factor(ifelse(nodedata_graph_2$clusterID=="grey", "grey","cluster"))
# nodedata_graph_2 <- melt(nodedata_graph_2, id.vars = c("clusterID","nodeID","node_type","KME","cortype","grey"))
# 
# ggplot(data = nodedata_graph_2, mapping = aes(x = KME, y = value, group = node_type, col = node_type, lty = grey)) + 
#   geom_smooth(method = "loess", mapping = aes(group = cortype)) +
#   scale_y_log10() + facet_grid(grey~variable) + theme_bw()