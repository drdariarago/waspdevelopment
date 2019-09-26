# Find genes and linkage clusters enriched for sexDE genes


##### binomial model
# load data and select an Fdr threshold
source('.//Scripts/BaseStats.R', echo=F)
sexdedata<-as.data.frame(apply(fit2$Fdr[,6:10],1, function(x){any(x<0.01)}))
names(sexdedata)<-c("sexDE")
# add chromosome and linkage group information
OGS2.0<- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")
sex_data<-merge(sexdedata,OGS2.0[,c("geneID","LinkageCluster","recombrate","cM","Kb","ODB6_OG_ID","adult_female_meth_status","quality7","isoforms","IntMatch","quality2","ratio")], by.x="row.names", by.y="geneID", all.x=T)
sex_data$LinkageCluster<-as.factor(sex_data$LinkageCluster)
sex_data$Chromosome<-as.factor(substr(sex_data$LinkageCluster, 1, 1))
sex_data$KbCM<-sex_data$Kb/sex_data$cM
# # remove genes with no linkage cluster or with chromosome alone
# sex_data<-droplevels(sex_data[-which(sex_data$LinkageCluster%in%c(1:5,"None")|is.na(sex_data$LinkageCluster)),])
# add number of genes on each cluster
sex_data<-merge(sex_data, as.data.frame(table(sex_data$LinkageCluster)), by.x="LinkageCluster", by.y="Var1", all.x=T)
names(sex_data)[ncol(sex_data)]<-"Ngenes"
# # remove clusters with less than 10 genes
# sex_data<-droplevels(sex_data[-which(sex_data$Ngenes<5),])

# Select only models with strong expression support
sex_data<-droplevels(sex_data[which(sex_data$quality2=="Express:Strong"),])
# # Select only models with orthology/paralogy assignments (genes with no ortholog/paralog have too few isoforms)
# mosaicplot(quality7~sexDE+(isoforms>1), data=sex_data)
# sex_data<-droplevels(sex_data[which(sex_data$quality7!="None"),])

# How much of OGS2 is present in the final model?
table(OGS2.0$geneID%in%sex_data$Row.names)
prop.table(table(OGS2.0$geneID%in%sex_data$Row.names))

# Fit maximal model
library(lme4)
binomfull<-glmer(sexDE~quality7+adult_female_meth_status+log10(Ngenes)+log10(recombrate+0.000001)+log10(IntMatch+1)+log10(ratio)+(1|LinkageCluster)+(1|ODB6_OG_ID), family="binomial", data=sex_data)

# previous AICc based data-exploration supports retaining OG. Support for Linkage cluster is low if ngenes and recombrate are present

# Calculate bootstrap based parameter estimates
# NOTE: bootMer cannot handle NAs, so remove those before fitting main model (data=na.omit(dataset))

# Parameter sampling function: sample all fixed effect estimates, variance estimates for random effects and sigma (returns error for sigma, so put it last to avoid loop skipping)
mySumm <- function(.) {
  c(beta=fixef(.),sig01=unlist(VarCorr(.)),sigma=sigma(.))
}

boo01<-bootMer(binomfull, mySumm, nsim=5, verbose=T, .progress="txt")

boo01

# Create list to store confidence intervals with different methods
bootCIs<-list(normal=as.data.frame(matrix(nrow=length(boo01$t0), ncol=3)),
              basic=as.data.frame(matrix(nrow=length(boo01$t0), ncol=5))
)

# Rename parameters
bootCIs<-lapply(bootCIs, function(x){
  row.names(x)<-names(boo01$t0)
  x
  })
# Rename variables
names(bootCIs$normal)<-c("alpha","lowerCI","upperCI")
names(bootCIs$basic)<-c("alpha","lowerEND","upperEND","lowerCI","upperCI")

# Calculate confidence intervals for each factor with normal and basic method (when error skips line, must skip before)
for (i in 1:length(boo01$t0)){
  a<-boot.ci(boo01, index=i, conf=0.95, type=c("norm","basic"))
  bootCIs[["normal"]][i,]<-a[["normal"]]
  bootCIs[["basic"]][i,]<-a[["basic"]]
}
# and plot resampled parameter densities (title not printing correctly)
pdf("./Graphics/Sex_bootparameters.pdf", paper="a4", onefile=T)
for (i in 1:length(boo01$t0)){
  print(plot(boo01, index=i),
        title(main=names(boo01$t0[i])))
}
dev.off()

# Save results to workspace
save("./Output/Sex_bootstrapGLMM")

# for (v in names(bootCIs)){
#   row.names(bootCIs[[v]])[i]<-row.names(a[[v]])
#   bootCIs[[v]][i,]<-a[[v]][,c(1,ncol(bootCIs[[v]])-1,ncol(bootCIs[[v]]))]
# }


# binomnull<-glmer(sexDE~1+(1|ODB6_OG_ID)+(1|LinkageCluster), family="binomial", data=sex_data)
# 


# binomfullglm<-glm(sexDE~0+quality7+adult_female_meth_status+log10(Ngenes)+log10(recombrate)+log10(ratio)+log10(IntMatch+1)+I(isoforms>1), family="binomial", data=sex_data)
# plot(binomfullglm)
# 
# library(MuMIn)
# AICc(binomfull, binomfull2, binomfull3, binomfullglm, binomnull)
# 
# #binomfull<-glmer(sexDE~LinkageCluster-1+(1|LinkageCluster), family="binomial", data=sex_data)
# 
# plot(binomfull)
# summary(binomfull2)
# 
# # model simplification
# library(MuMIn)
# binomset<-dredge(binomfull2)
# model.sel(binomset)
# summary(model.avg(binomset))
# 
# binomnull<-glm(sexDE~1, family="binomial", data=sex_data)
# AICc(binomfull, binomnull)
# 
# # get approximate 95% CIs
# se<-sqrt(diag(vcov(binomfull)))
# tab<-data.frame(Est=binomfull$coefficients, LL=binomfull$coefficients-1.96*se, UL=binomfull$coefficients+1.96*se)
# tab$sign<-apply(tab, 1, function(x){
#   sign(x[2])+sign(x[3])
# })
# table(tab$sign)
# prop.table(table(tab$sign))
# 
# # Only 2 clusters are enriched in sex DE genes, but several are depleted
# # Show genes on the 2 enriched clusters
# sex_data[which(sex_data$LinkageCluster%in%c("3.018","4.096")),]
# 
# # Validation of significant linkage groups by bootstrapping
# library(boot)
# # Using the resampling with replacement of datapoints (genes)
# model.boot<-function(sex_data,i){
#   sub.data<-sex_data[i,]
#   model<-glm(sexDE~LinkageCluster, family="binomial", data=sub.data)
#   coef(model)
# }
# glm.boot<-boot(sex_data, model.boot, R=10)
# boot.ci(glm.boot)
# # using the shuffle-residuals method (Crawley, page 525)
# yhat<-fitted(binomfull)
# resids<-resid(binomfull)
# res.data<-data.frame(resids, sex_data$LinkageCluster)
# names(res.data)<-c("resids", "LinkageCluster")
# bf<-function(res.data,i){
#   y<-yhat+res.data[i,1]
#   nd<-data.frame(y,res.data$LinkageCluster)
#   model<-glm(y~res.data.LinkageCluster, data=nd)
#   coef(model)}
# perms<-boot(res.data, bf, R=10, sim="permutation")

## Outdated trials
# # load data
# nasoniadevgene<- read.csv("./Output/nasoniadevgene.csv", header=T)
# row.names(nasoniadevgene)<-nasoniadevgene[,1]
# nasoniadevgene<-nasoniadevgene[,-1]
# # reformat for GLM analyses
# tnasoniadevgene<-as.data.frame(t(nasoniadevgene))
# # extract stage and sex
# tnasoniadevgene$stage<-factor(as.character(regmatches(row.names(tnasoniadevgene), gregexpr("^[[:alnum:]]+",row.names(tnasoniadevgene)))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
# tnasoniadevgene$male<-factor(grepl("[^[:alnum:]]male|testes", row.names(tnasoniadevgene)))
# # convert into flat format
# library(reshape)
# flatnasoniadevgene<-melt(tnasoniadevgene, variable_name="geneID")
# write.csv(flatnasoniadevgene, file="./Output/nasoniadevgene_flat.csv")
# 
# # add chromosome and linkage group information
# OGS2.0<- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")
# sex_data<-merge(flatnasoniadevgene,OGS2.0[,c("geneID","LinkageCluster","recombrate","cM","Kb")], by="geneID", all.x=T)
# sex_data$LinkageCluster<-as.factor(sex_data$LinkageCluster)
# sex_data$Chromosome<-as.factor(substr(sex_data$LinkageCluster, 1, 1))
# sex_data$KbCM<-sex_data$Kb/sex_data$cM
# 
# # # continuous model (too many levels in matrix)
# # library(lme4)
# # linkage_full<-lmer(value~male*Chromosome/LinkageCluster+(1|stage), data=sex_data)
# # summary(linkage_full)
# # library(MuMIn)
# # linkage_set<-dredge(linkage_full)
# # 
# # linkage_avg<-model.avg(linkage_set)
# # linkage_sel<-model.sel(linkage_set)
# # 
# # # write results on table 
# # write.csv(linkage_set, file="./Output/linkage_sel.csv")
# 
# 
# ### limma-style trial (test individual LGs for sex*stage DE)
# 
# # data import and formatting
# sex_data$replicate<-
#   cast(melt(head(sex_data[,1:5])),formula=LinkageCluster~geneID+stage+male)
# cast(melt(linkage~geneID+sample))
# 
# model.matrix(~0+LinkageCluster*recombrate, data=???)
# 
# # write results on table and workspace
# save.image(file="./Output/linkage_analyses.csv")
