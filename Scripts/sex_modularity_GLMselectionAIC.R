# initialize script
date()
rm(list=ls())
newdir<-file.path(getwd(), "Output/sex_modularity_GLMselectionAIC")
graphicsdir<-file.path(getwd(), "Graphics/sex_modularity_GLMselectionAIC")
dir.create(newdir)
dir.create(graphicsdir)
library(MuMIn)
library(plyr)
options(na.action = "na.fail")   #  prevent fitting models to different datasets

# load permutations
modulardensities <- read.csv("./Output/sex_modulewise_modularity/sex_modulewise_modulardensities.csv")

# hotfix: five permutations have a negative or greater than zero dWithin, replacing with NAs
modulardensities$dWithin[which(modulardensities$dWithin<0)]<-NA
modulardensities$dWithin[which(modulardensities$dWithin>1)]<-NA
modulardensities$dOut[which(modulardensities$dOut<0)]<-NA
modulardensities$dOut[which(modulardensities$dOut>1)]<-NA
# and remove rows with NAs
modulardensities<-na.exclude(modulardensities)

# split by module into list of data.frames
modulardensities<-dlply(.data = modulardensities[,-1], .variables = .(clusterID))


## Begin version specific for dWithin
# generate list of one glm per data.frame

modulenames<-names(modulardensities)
dWithinGLMs<-as.list(modulenames)
names(dWithinGLMs)<-modulenames

for(i in modulenames){
  dWithinGLMs[[i]]<-do.call(glm, args = list(
    formula = dWithin~stage-1+stage:sexbias+dMain+dOut,
    family = gaussian(link="logit"),
    data = modulardensities[[i]])
  )
}

# average every model in list
dWithin.avg<-lapply(dWithinGLMs, function(x){
  dredge(x)
  y<-dredge(x)
  model.avg(y)
})

# save effect sizes
dWithin.avg2<-lapply(X = dWithin.avg, FUN = function(x){summary(x)$avg.model})
factorName<-unlist(lapply(dWithin.avg2, row.names))
dWithin.avg2<-cbind(factorName,ldply(dWithin.avg2))
names(dWithin.avg2)[2]<-"clusterID"

# calculate whether effects are significant (if 95% intervals overlap zero not significant)
dWithin.avg2$significant<-sign(dWithin.avg2$"Upper CI"*dWithin.avg2$"Lower CI")==1

# Save results as csv
write.csv(dWithin.avg2, file = file.path(newdir, "dWithin_avg.csv"))


## Begin version specific for dOut
# generate list of one glm per data.frame

modulenames<-names(modulardensities)
dOutGLMs<-as.list(modulenames)
names(dOutGLMs)<-modulenames

for(i in modulenames){
  dOutGLMs[[i]]<-do.call(glm, args = list(
    formula = dOut~stage-1+stage:sexbias+dMain+dWithin,
    family = gaussian(link="logit"),
    data = modulardensities[[i]])
  )
}

# average every model in list
dOut.avg<-lapply(dOutGLMs, function(x){
  y<-dredge(x)
  model.avg(y)
})

# save effect sizes
dOut.avg2<-lapply(X = dOut.avg, FUN = function(x){summary(x)$avg.model})
factorName<-unlist(lapply(dOut.avg2, row.names))
dOut.avg2<-cbind(factorName,ldply(dOut.avg2))
names(dOut.avg2)[2]<-"clusterID"

# calculate whether effects are significant (if 95% intervals overlap zero not significant)
dOut.avg2$significant<-sign(dOut.avg2$"Upper CI"*dOut.avg2$"Lower CI")==1

# Save results as csv
write.csv(dOut.avg2, file = file.path(newdir, "dOut_avg.csv"))

#####


# run model selection on each element of the list
# save remaining factors and effect sizes
