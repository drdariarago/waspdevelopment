# test appropriate distribution for glm based analyses
modulardensities <- read.csv("./Output/sex_modulewise_modularity/sex_modulewise_modulardensities.csv")
library(robustbase)
library(lme4)

color<-levels(droplevels(modulardensities$clusterID))

bestmodel<-list(color)

param="dOut"
for(f in 1:length(color)){
  glmgausslogit<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),],family=gaussian(link="logit"))
  #   glmgausslog<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),],family=gaussian(link="log"))
  #   glmgaussidentity<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),], family=gaussian)
  glmGamma<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),], family=Gamma)
#   glmPoislogit<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),], family=poisson(link="logit"))
  
#   glmGammalog<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),], family=Gamma(link="log"))
  #   glmQuasibinomial<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),], subset = get(param)>0, family=quasibinomial)
  AICtable<-AIC(glmGamma, glmgausslogit, glmPoislogit)
  #   bestmodel[f]<-AICtable
  bestmodel[[f]]<-row.names(AICtable)[which(AICtable$AIC<=(min(AICtable$AIC)+2))]
}
table(unlist(bestmodel))/length(bestmodel)
table(unlist(lapply(X = bestmodel, FUN=paste0, collapse="")))/length(bestmodel)
# prop.table(table(unlist(bestmodel)))

color<-levels(droplevels(modulardensities$clusterID))
bestmodel<-list()
param="dWithin"
for(f in 1:length(color)){
  glmgausslogit<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),],family=gaussian(link="logit"))
#   glmgausslog<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),],family=gaussian(link="log"))
  #   glmgaussidentity<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),], family=gaussian)
  glmGamma<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),], family=Gamma)
#   glmGammalog<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),], family=Gamma(link="log"))
#   glmQuasibinomial<-glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==color[f]),], subset = get(param)>0, family=quasibinomial)
  AICtable<-AIC(glmGamma, glmgausslogit)
  bestmodel[[f]]<-as.vector(row.names(AICtable)[which(AICtable$AIC<=(min(AICtable$AIC)+2))])
}
table(unlist(bestmodel))/length(bestmodel)
table(unlist(lapply(X = bestmodel, FUN=paste0, collapse="")))/length(bestmodel)
# prop.table(table(unlist(bestmodel)))

# color<-levels(modulardensities$clusterID)
# param="dTotal"
# 
# diagnostics<-list()
# glmlist<-lapply(color, function(x){glm(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==x),], subset = get(param)>0, family=gaussian(link="logit"))})
# reslist<-llply(.data = glmlist, .fun = function(x){cbind(residuals(x), fitted(x))})
# names(reslist)<-color
# resdata<-ldply(reslist)
# names(resdata)<-c("clusterID", "residuals", "fitted")
# ggplot(data = resdata, aes(x=residuals, y=fitted, group=clusterID))+theme_bw()+geom_point()+geom_smooth()+facet_wrap(~clusterID)
# ggplot(data = resdata, aes(sample=z(residuals)))+theme_bw()+facet_wrap(~clusterID)+stat_qq()+geom_abline(aes(a=0,b=1))

color<-levels(droplevels(modulardensities$clusterID))
param="dWithin"

diagnostics<-list()
glmlist<-lapply(color, function(x){glmrob(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==x),],subset = get(param)>0, family=Gamma)})
reslist<-llply(.data = glmlist, .fun = function(x){cbind(residuals(x), fitted(x))})
names(reslist)<-color
resdata<-ldply(reslist)
names(resdata)<-c("clusterID", "residuals", "fitted")
ggplot(data = resdata, aes(x=residuals, y=fitted, group=clusterID))+theme_bw()+geom_point()+geom_smooth()+facet_wrap(~clusterID)
qqmath(~resdata$residuals|resdata$clusterID)

param="dOut"

diagnostics<-list()
glmlist<-lapply(color, function(x){glmrob(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==x),],subset = get(param)>0, family=Gamma)})
reslist<-llply(.data = glmlist, .fun = function(x){cbind(residuals(x), fitted(x))})
names(reslist)<-color
resdata<-ldply(reslist)
names(resdata)<-c("clusterID", "residuals", "fitted")
ggplot(data = resdata, aes(x=residuals, y=fitted, group=clusterID))+theme_bw()+geom_point()+geom_smooth()+facet_wrap(~clusterID)
qqmath(~resdata$residuals|resdata$clusterID)

param="dTotal"

diagnostics<-list()
glmlist<-lapply(color, function(x){glmrob(formula = get(param)~stage*sexbias, data = modulardensities[which(modulardensities$clusterID==x),],subset = get(param)>0, family=Gamma)})
reslist<-llply(.data = glmlist, .fun = function(x){cbind(residuals(x), fitted(x))})
names(reslist)<-color
resdata<-ldply(reslist)
names(resdata)<-c("clusterID", "residuals", "fitted")
ggplot(data = resdata, aes(x=residuals, y=fitted, group=clusterID))+theme_bw()+geom_point()+geom_smooth()+facet_wrap(~clusterID)
qqmath(~resdata$residuals|resdata$clusterID)


## test for mean/variance relationships
meanvartest<-ddply(.data = modulardensities, .variables = .(clusterID, stage, sexbias), .fun = summarize, withinmean=mean(dWithin, na.rm = T), withinvariance=var(dWithin, na.rm = T), outmean=mean(dOut, na.rm=T), outvariance=var(dOut,na.rm = T))
ggplot(data = meanvartest, aes(x=withinmean, y=withinvariance))+geom_point()+scale_y_log10()+geom_smooth()+scale_x_log10()
ggplot(data = meanvartest, aes(x=outmean, y=outvariance))+geom_point()+scale_y_log10()+geom_smooth()+scale_x_log10()
