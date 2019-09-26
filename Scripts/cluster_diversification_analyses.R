## Analyze how duplication and AS assist transcriptional diversification of OGs
## We analyze transcriptional diversification as number of clusters per OG

# Clear workspace
rm(list=ls())
# Initialize script
library(reshape2)
library(plyr)
library(stringr)
library(vcd)
library(ggplot2)
library(lattice)
library(pscl)
library(MuMIn)
options(na.action="na.fail")
newdir<-file.path(getwd(),"Output/OG_diversification_analyses")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/OG_diversification_analyses")
dir.create(graphdir)
source(file = "./Scripts/multiplot.R")

# Load data
clusterdata<-read.csv(file = "./Output/Results_compiler//clusterdata_compiled.csv")[,-1]

# Count number of clusters per OG (assesses diversification between paralogs)

cluster_paralog_count<-ddply(clusterdata, .(odb8_og_id), summarize, 
                             nClusters = length(unique(clusterID)), 
                             nParalogs = length(unique(geneID)),
                             nIsoforms = length(unique(eigenexonID)),
                             .progress = "text")

# Remove clusters with no OG
cluster_paralog_count<-na.exclude(cluster_paralog_count)
# # Remove 1 from all numeric entries (normalizes model to baseline status of 1 cluster, 1 gene, 1 isoform)
# transcript beyond the first transcript per gene
cluster_paralog_count$nIsoforms<-cluster_paralog_count$nIsoforms-cluster_paralog_count$nParalogs
# Clusters above the first
cluster_paralog_count$nClusters<-cluster_paralog_count$nClusters-1
# Additional gene copies
cluster_paralog_count$nParalogs<-cluster_paralog_count$nParalogs-1

# Limit to subset of clusters with parameters within 0.05 and 0.95th quantiles
apply(cluster_paralog_count[,-1], 2, quantile, probs=c(0.01, 0.99))
cluster_paralog_count<-cluster_paralog_count[which(apply(ddply(cluster_paralog_count, .(odb8_og_id), summarize, nClusters<=10, nParalogs<=8, nIsoforms<=20)[,-1],1,function(x)all(x==T))),]

# Plot distribution
ggplot(cluster_paralog_count, aes(x=(nParalogs+nIsoforms)/nClusters))+geom_density()+theme_bw()
ggplot(cluster_paralog_count, aes(x=(nParalogs+nIsoforms)/nClusters))+geom_histogram()+theme_bw()+scale_x_log10()

# Analyze diversification as GLM (poisson events?)
nClust_model<-glm(nClusters~-1+log10(nParalogs+1)*log10(nIsoforms+1), data = cluster_paralog_count, family = poisson)
nClust_model_logiso<-glm(nClusters~-1+nParalogs*log10(nIsoforms+1), data = cluster_paralog_count, family = poisson)
nClust_model_logpara<-glm(nClusters~-1+log10(nParalogs+1)*nIsoforms, data = cluster_paralog_count, family = poisson)
AICc(nClust_model,nClust_model_logiso, nClust_model_logpara)

summary(nClust_model_logiso)
plot(nClust_model_logiso)

model.sel(dredge(nClust_model_logiso))

## Subselect only OGs present in more than 1 module
cluster_paralog_count_2<-cluster_paralog_count[which(cluster_paralog_count$nClusters>0),]

# Analyze diversification as GLM (poisson events?)
nClust_model_2<-glm(nClusters~-1+log10(nParalogs+1)*log10(nIsoforms+1), data = cluster_paralog_count_2, family = poisson)
nClust_model_logiso_2<-glm(nClusters~-1+nParalogs*log10(nIsoforms+1), data = cluster_paralog_count_2, family = poisson)
nClust_model_logpara_2<-glm(nClusters~-1+log10(nParalogs+1)*nIsoforms, data = cluster_paralog_count_2, family = poisson)
AICc(nClust_model_2,nClust_model_logiso_2, nClust_model_logpara_2)

summary(nClust_model_logiso_2)
plot(nClust_model_logiso_2)

# distribution is now regular, suggesting two processes that lead to initial differentiation of OG and then increased divergence potential
# Trying with hurdle based model

nClust_model_hurdle<-hurdle(nClusters~-1+log10(nParalogs+1)*log10(nIsoforms+1), data = cluster_paralog_count, dist = "poisson", zero.dist = "binomial")
nClust_model_hurdle_logiso<-hurdle(nClusters~-1+nParalogs*log10(nIsoforms+1), data = cluster_paralog_count, dist = "poisson", zero.dist = "binomial")
nClust_model_hurdle_logpara<-hurdle(nClusters~-1+log10(nParalogs+1)*nIsoforms, data = cluster_paralog_count, dist = "poisson", zero.dist = "binomial")
nClust_model_hurdle_nolog<-hurdle(nClusters~-1+poly(nParalogs, degree = 5)*poly(nIsoforms, degree = 5), data = cluster_paralog_count, dist = "poisson", zero.dist = "binomial")

AICc(nClust_model_hurdle, nClust_model_hurdle_logiso, nClust_model_hurdle_logpara, nClust_model_hurdle_nolog)

pdf(file.path(graphdir, "qqplots.pdf"))
par(mfrow = c(1,3))
qqnorm(residuals(nClust_model_hurdle), ylim = c(-5,5), xlim = c(-5,5))
abline(0,1)
qqnorm(residuals(nClust_model_hurdle_logiso), ylim = c(-5,5), xlim = c(-5,5))
abline(0,1)
qqnorm(residuals(nClust_model_hurdle_logpara), ylim = c(-5,5), xlim = c(-5,5))
abline(0,1)
par(mfrow = c(1,1))
dev.off()

AICc(nClust_model_logiso, nClust_model_hurdle_logiso, nClust_model_hurdle_logpara)

summary(nClust_model_hurdle_logpara)

hurdledata_logpara<-as.data.frame(as.matrix(expand.grid(0:8,0:20)))
names(hurdledata_logpara)<-c("nParalogs","nIsoforms")
hurdledata_logpara$nClusters<-predict(object = nClust_model_hurdle_logpara,
                              newdata = hurdledata,
                              type="response"
)

hurdledata_logiso<-as.data.frame(as.matrix(expand.grid(0:8,0:20)))
names(hurdledata_logiso)<-c("nParalogs","nIsoforms")
hurdledata_logiso$nClusters<-predict(object = nClust_model_hurdle_logiso,
                              newdata = hurdledata,
                              type="response"
)

nClusters_obs<-ggplot(data = cluster_paralog_count, aes(x=nParalogs, y=nIsoforms, fill=nClusters))+geom_raster()+theme_bw()+
  scale_fill_distiller(type="div",  palette = 9, breaks=c(1,(1:5)*2), limits=c(0,10))+
  scale_x_continuous(breaks=c((0:4)*2), limits=c(0,8))+
  scale_y_discrete(breaks=c((0:10)*2), limits=c(0:20))

nClusters_logpara<-ggplot(data = hurdledata_logpara, aes(x=nParalogs, y=nIsoforms, fill=nClusters))+geom_raster()+theme_bw()+
  scale_fill_distiller(type="div",  palette = 9, breaks=c(1,(1:5)*2), limits=c(0,10))+
  scale_x_continuous(breaks=c((0:4)*2), limits=c(0,8))+
  scale_y_discrete(breaks=c((0:10)*2), limits=c(0:20))

nClusters_logiso<-ggplot(data = hurdledata_logiso, aes(x=nParalogs, y=nIsoforms, fill=nClusters))+geom_raster()+theme_bw()+
  scale_fill_distiller(type="div",  palette = 9, breaks=c(1,(1:5)*2), limits=c(0,10))+
  scale_x_continuous(breaks=c((0:4)*2), limits=c(0,8))+
  scale_y_discrete(breaks=c((0:10)*2), limits=c(0:20))

pdf(file.path(graphdir, "model_projection_hurdle.pdf"), width = 8, height = 4)
multiplot(nClusters_obs, nClusters_logiso, nClusters_logpara, cols = 1)
dev.off()

poisdata<-as.data.frame(as.matrix(expand.grid(0:8,0:20)))
names(poisdata)<-c("nParalogs","nIsoforms")
poisdata$nClusters<-predict(object = nClust_model_logiso,
                            newdata = poisdata,
                            type="response"
)

nClusters_exp_pois<-ggplot(data = poisdata, aes(x=nParalogs, y=nIsoforms, fill=nClusters))+geom_raster()+theme_bw()+
  scale_fill_distiller(type="div",  palette = 9, breaks=c(1,(1:5)*2), limits=c(0,10))+
  scale_x_continuous(breaks=c((0:4)*2), limits=c(0,8))+
  scale_y_discrete(breaks=c((0:10)*2), limits=c(0:20))

pdf(file.path(graphdir, "model_projection_pois.pdf"), width = 8, height = 4)
multiplot(nClusters_exp_pois, nClusters_obs, cols = 2)
dev.off()

### Outdated code snippets
# OG EOG6G79G8 and EOG60GB79 are strong outliers with huge number of paralogs and clusters, maps to retrotranscriptases, endo/exonucleases-phosphatases therefore most likely a family of retrotransposons
# Excluding them from further analyses, no longer necessary after data trim
# cluster_paralog_count<-cluster_paralog_count[-which(cluster_paralog_count$ODB6_OG_ID%in%c("EOG6G79G8","EOG60GB79")),]
