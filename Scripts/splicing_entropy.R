## Comparison of information carried by genes rather than transcripts
## Idea 1, calculate ratio of shannon entropy explained by genes rather than transcripts in each sample, then model as glm

# Initialize script
date()
rm(list=ls())
library(entropy)
library(plyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(MuMIn)
options(na.action = "na.fail")
source('./Scripts/multiplot.R', echo=F)
# initialize output path
newdir<-file.path(getwd(), "Output/splicing_entropy")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/splicing_entropy")
dir.create(graphdir)
# Load and reshape dataset

eigenexon_evalues<-read.csv(file = "./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexon_evalues.csv")
eigenexon_evalues$splicing<-grepl(pattern = "fac", eigenexon_evalues$eigenexonID)
# Replace NAs with zeroes
eigenexon_evalues_2<-apply(eigenexon_evalues[,-1], c(1,2), function(x)(ifelse(is.na(x),0,x)))
# Binarize expression of genes and isoforms
eigenexon_evalues_2<-apply(eigenexon_evalues_2, c(1,2), function(x)(ifelse(x>0,1,0)))

# And save gene names with splicing status as a data.frame
row.names(eigenexon_evalues_2)<-eigenexon_evalues$eigenexonID
eigenexon_evalues_2<-as.data.frame(eigenexon_evalues_2)

# Store number of genes and isoforms for normalization
nGenes<-table(eigenexon_evalues_2$splicing)[1]
nIsoforms<-table(eigenexon_evalues_2$splicing)[2]

# # Calculate information content for every gene/isoform based on overall frequency
# inf_coef<-apply(eigenexon_evalues_2[,-grep("splicing",names(eigenexon_evalues_2))], 1, function(x){
#   y<-sum(x>0)/30
#   ifelse(y==0, 0, log2(1/y))
#   })

# Calculate information coefficient based on sample type (gene expressed if present in more than 1 sample of the same type)
inf_coef<-as.data.frame(t(eigenexon_evalues_2[,-grep("splicing",names(eigenexon_evalues_2))]))
inf_coef$sampletype<-str_extract(string = row.names(inf_coef), pattern = "[[:alnum:]]*_[[:alnum:]]*")
# Use two out of three binary expression principle
inf_coef<-daply(inf_coef, .variables = .(sampletype), colwise(function(x){sum(x)>1}), .progress = "text")

# # Using empirical estimator
# inf_coef<-apply(inf_coef, 2, function(x){
#   y<-sum(x>0)/10
#   ifelse(y==0, 0, log2(1/y))
# })
# 
# # Using Miller correction
# inf_coef<-apply(inf_coef, 2, function(x){
#   y<-sum(x>0)/10
#   ifelse(y==0, 0, log2(1/y)+1/(2*10))
# })

## Using estimators corrected for small sample size (package entropy)
inf_coef<-ifelse(inf_coef==0,0,1)
# inf_freqs<-apply(inf_coef,2,freqs, method="CS")
inf_freq<-apply(inf_coef, 2, function(x){sum(x>0)})
# We double the dataset to avoid numerical problems with single entries, then divide by two to normalize the scores
inf_coef<-rbind(inf_coef, inf_coef)
inf_entr<-apply(X = inf_coef, MARGIN = 2, FUN = entropy, method="CS", unit = "log2")
inf_entr<-inf_entr/2
# inf_coef<-inf_coef/30
# inf_coef<-ifelse(is.na(inf_coef),0,inf_coef)

# Add information coefficient to dataset 
eigenexon_evalues_2$inf_coef<-inf_entr/inf_freq

# Save self-information for every eigenexon
self_information<-data.frame(eigenexonID=row.names(eigenexon_evalues_2), inf_coef=eigenexon_evalues_2$inf_coef)
write.csv(self_information, file=file.path(newdir, "eigenexon_self-information.csv"))

# Remove nodes with zero expression across all samples 
eigenexon_evalues_2<-na.exclude(eigenexon_evalues_2)

# Plot information coefficients of genes vs isoforms
ggplot(data = eigenexon_evalues_2, aes(x=inf_coef, col=as.factor(splicing)))+geom_density()+theme_bw()
# ggplot(data = eigenexon_evalues_2[which(eigenexon_evalues_2$inf_coef>0),], aes(x=inf_coef, col=as.factor(splicing)))+geom_density()+theme_bw()

## Boxplot of information coefficients of transcripts above zero
pdf(file = file.path(graphdir,"Inf_coefs_boxplot.pdf"), paper = "a4")
ggplot(data = eigenexon_evalues_2[which(eigenexon_evalues_2$inf_coef>0),], aes(y=inf_coef, x=as.factor(splicing)))+geom_boxplot(notch = T, varwidth = T)+theme_bw()+stat_summary(fun.y="mean", geom="point", shape=23, size=3)+xlab(label = "Splicing category")+ylab(label = "Self Information \n (bits per occurrence)")+ggtitle(label = "Distribution of information coefficients across informative nodes\n")+scale_x_discrete(labels=c("Expression","Splicing"))
dev.off()

pdf(file = file.path(graphdir,"Inf_coefs_hist.pdf"), height = 10, width = 15)
ggplot(data = eigenexon_evalues_2, aes(x=inf_coef, fill=as.factor(splicing)))+geom_histogram(position="dodge")+theme_bw()+xlab(label = "Self Information \n (bits per occurrence)")+ggtitle(label = "Distribution of self-information \n in constitutive and splicing nodes")+scale_fill_discrete(labels = c("Constitutive","Alternative"), name = "Node Type")+theme(legend.position="bottom")
dev.off()

ggplot(data = eigenexon_evalues_2, aes(x=inf_coef>0.45, fill=as.factor(splicing)))+geom_histogram(position="dodge")+theme_bw()+xlab(label = "Self Information \n (bits per occurrence)")+ggtitle(label = "Distribution of information coefficients across informative nodes\n")+scale_x_discrete(labels = c("Constitutive","Restricted"))

# Plot tabulations of self-information vs splicing
library(vcd)
eigenexon_evalues_2_vcd<-eigenexon_evalues_2
eigenexon_evalues_2_vcd$splicing<-as.factor(ifelse(eigenexon_evalues_2_vcd$splicing==0,"Gene","Splicing"))
eigenexon_evalues_2_vcd$restricted<-as.factor(ifelse(eigenexon_evalues_2_vcd$inf_coef>0.45,"Restricted","Aspecific"))
eigenexon_evalues_2_vcd$inf_coef<-as.factor(as.numeric(as.factor(eigenexon_evalues_2_vcd$inf_coef)))

pdf(file = file.path(graphdir, "cotab_splicing_restriction_exp.pdf"), width = 7, height = 7)
cotabplot(~restricted+splicing, data=na.exclude(eigenexon_evalues_2_vcd), shade=T, type="expected", 
          labeling_args=list(
            rot_labels = c(left = 0),
            rot_varnames = c(left = 0))
          )
dev.off()

pdf(file = file.path(graphdir, "cotab_splicing_restriction_obs.pdf"), width = 7, height = 7)
cotabplot(~restricted+splicing, data=na.exclude(eigenexon_evalues_2_vcd), shade=T, 
          labeling_args=list(
            rot_labels = c(left = 0),
            rot_varnames = c(left = 0))
)
dev.off()

pdf(file = file.path(graphdir, "cotab_splicing_self_inf_exp.pdf"), width = 7, height = 7)
cotabplot(~inf_coef+splicing, data=eigenexon_evalues_2_vcd, shade=T, type="expected", 
          labeling_args=list(
            rot_labels = c(left = 0),
            rot_varnames = c(left = 0))
)
dev.off()

pdf(file = file.path(graphdir, "cotab_splicing_self_inf_obs.pdf"), width = 7, height = 7)
cotabplot(~inf_coef+splicing, data=eigenexon_evalues_2_vcd, shade=T, 
          labeling_args=list(
            rot_labels = c(left = 0),
            rot_varnames = c(left = 0))
)
dev.off()

# # How many isoforms and genes have zero information potential?
# ddply(eigenexon_evalues_2, .variables = .(splicing), .fun = summarize, zeroes=sum(inf_coef==0), nonzeroes=sum(inf_coef!=0), prop.zeroes=sum(inf_coef==0)/length(inf_coef))

# Calculate gene and splicing information in each sample
inf_sample<-rbind(
  apply(eigenexon_evalues_2[eigenexon_evalues_2$splicing=="0",-grep("splicing",names(eigenexon_evalues_2))],2,function(x){sum(x*inf_coef[which(eigenexon_evalues_2$splicing=="0")])}),
  apply(eigenexon_evalues_2[eigenexon_evalues_2$splicing=="1",-grep("splicing",names(eigenexon_evalues_2))],2,function(x){sum(x*inf_coef[which(eigenexon_evalues_2$splicing=="1")])})
)
# Plot proportions of entropy explained by genes vs isoforms
cbind(inf_sample[1,]/apply(inf_sample,2,sum),inf_sample[2,]/apply(inf_sample,2,sum))
apply(cbind(inf_sample[1,]/apply(inf_sample,2,sum),inf_sample[2,]/apply(inf_sample,2,sum))[-31,],2,mean)
# what are the mean differences from the raw expected proportions?
apply(cbind(inf_sample[1,]/apply(inf_sample,2,sum),inf_sample[2,]/apply(inf_sample,2,sum))[-31,]/cbind(inf_sample[1,]/apply(inf_sample,2,sum),inf_sample[2,]/apply(inf_sample,2,sum))[31,],2,mean)
# What are the mean information density of genes and isoforms?
ddply(eigenexon_evalues_2, .(splicing), summarize, mean=mean(inf_coef), median=median(inf_coef))

# Convert to data.frame for further analyses
eigenexon_entropy<-as.data.frame(inf_sample)[,-31]
row.names(eigenexon_entropy)<-c("genes","isoforms")
eigenexon_entropy$splicing<-c("genes","isoforms")

# Reshape data for plotting and glms
eigenexon_entropy<-melt(eigenexon_entropy, id.vars = "splicing", value.name = "entropy")
eigenexon_entropy$stage<-as.factor(str_extract(string = eigenexon_entropy$variable, pattern = "[[:alnum:]]*"))
eigenexon_entropy$stage<-factor(eigenexon_entropy$stage, levels = c("emb10","emb18","lar51","pupyel","adult"))
eigenexon_entropy$sex<-ifelse(grepl(pattern = "female", x = eigenexon_entropy$variable),"female","male")
eigenexon_entropy$orderedstage<-as.ordered(eigenexon_entropy$stage)
eigenexon_entropy$splicing<-as.factor(eigenexon_entropy$splicing)
eigenexon_entropy$group<-as.factor(paste(eigenexon_entropy$splicing,eigenexon_entropy$sex,eigenexon_entropy$stage))
eigenexon_entropy$group2<-as.factor(paste(eigenexon_entropy$splicing,eigenexon_entropy$sex))

# Plot entropy of genes and alternative isoforms over development  in each sex
ggplot(data = eigenexon_entropy, aes(y=entropy, x=stage, col=sex, pch=splicing))+geom_point()+theme_bw()+geom_smooth(aes(group=group2), se=F)

ggplot(data = eigenexon_entropy, aes(y=entropy, x=stage, col=sex, group=group, linetype=splicing))+geom_boxplot(varwidth = T)+theme_bw()+geom_smooth(aes(group=group2), se=F)+ylab(label = "Entropy (bits)")+xlab(label = "Stage")+ggtitle(label = "Information encoded by transcriptional and gene nodes in each stage\n")

pdf(file = file.path(graphdir, "information_treatments.pdf"), paper = "a4")
ggplot(data = eigenexon_entropy, aes(y=entropy, x=stage, col=sex, group=group, pch=splicing))+geom_point()+theme_bw()+geom_smooth(aes(group=group2), se=F)+ylab(label = "Entropy (bits)")+xlab(label = "Stage")+ggtitle(label = "Information encoded by expression and splicing nodes in each stage\n")+scale_shape_discrete(name="Node type", labels=c("Expression","Splicing"))+scale_color_discrete(name="Sex", label=c("Female","Male"))
dev.off()

# Test for significance of entropy changes

dataplot<-ggplot(data = eigenexon_entropy, aes(y=entropy, x=stage, col=sex, group=group, pch=splicing))+geom_point()+theme_bw()+geom_smooth(aes(group=group2), se=F)+ylab(label = "Entropy (bits)")+xlab(label = "Stage")+ggtitle(label = "Information encoded by expression and splicing nodes \n in each stage from raw data")+scale_shape_discrete(name="Node type", labels=c("Expression","Splicing"))+scale_color_discrete(name="Sex", label=c("Female","Male"))+ylim(7500,22000)+theme(legend.position="bottom", panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())#+facet_grid(splicing~.)

dataplot

# Check for differences between splicing and gene other than the intercept
model4<-glm(formula = entropy~poly(as.numeric(orderedstage),degree = 4)/sex*splicing, data = eigenexon_entropy)
model3<-glm(formula = entropy~poly(as.numeric(orderedstage),degree = 3)/sex*splicing, data = eigenexon_entropy)
model2<-glm(formula = entropy~poly(as.numeric(orderedstage),degree = 2)/sex*splicing, data = eigenexon_entropy)
model1<-glm(formula = entropy~poly(as.numeric(orderedstage),degree = 1)/sex*splicing, data = eigenexon_entropy)
model0<-glm(formula = entropy~1, data = eigenexon_entropy)

mainmodellist<-ldply(lapply(list(model1, model2, model3, model4, model0), FUN = dredge))
mainmodellist<-mainmodellist[order(mainmodellist$AICc),]
head(mainmodellist)

dredge(model4)
finalmodel<-get.models(dredge(model4), subset = delta==0)
finalmodel<-finalmodel[[1]]
summary(finalmodel)

plot(finalmodel)

predictions2<-data.frame(orderedstage = rep(1:5,4),
                         splicing = c(rep("genes",10),rep("isoforms",10)),
                         sex=rep(c(rep("male",5),rep("female",5)),2)
)

pred_values<-predict(finalmodel, se.fit = T, 
                     newdata = predictions2,
                     type="response")
pred_values<-as.data.frame(pred_values)

pred_values$se95<-pred_values$fit+1.96*pred_values$se.fit
pred_values$se05<-pred_values$fit-1.96*pred_values$se.fit

predictions2<-cbind(predictions2, pred_values)

predplot<-ggplot(data = predictions2, aes(x=orderedstage, y=fit, ymin=se05, ymax=se95, col=sex, pch=splicing))+geom_pointrange()+geom_line()+theme_bw()+ylab(label = "Entropy (bits)")+xlab(label = "Stage")+ggtitle(label = "Information encoded by expression and splicing nodes \n in each stage from best fitting models")+scale_shape_discrete(name="Node type", labels=c("Expression","Splicing"))+scale_color_discrete(name="Sex", label=c("Female","Male"))+theme(legend.position="bottom", panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())+scale_x_discrete(labels=levels(eigenexon_entropy$orderedstage))+ylim(7500,22000)

pdf(file=file.path(graphdir, "entropy_TS_joint_models.pdf"), height = 10, width = 11)
multiplot(dataplot, predplot, cols = 2)
dev.off()


#### Outdated code snippets

# # Add offset (proportion of gene to splicing nodes)
# propGenes<-(table(eigenexon_evalues_2$splicing))[1]
# propSplicing<-(table(eigenexon_evalues_2$splicing))[2]
# eigenexon_entropy$props<-ifelse(eigenexon_entropy$splicing=="genes",propGenes,propSplicing)
# 
# # Run model using stage as a factor
# fullmodel<-lm(formula = entropy~(stage-1)/sex*splicing, data = eigenexon_entropy, offset = props)
# plot(fullmodel)
# summary(fullmodel)
# model.sel(dredge(fullmodel))
# summary(model.avg(dredge(fullmodel)))
# # Run model using stage as an ordered factor
# 
# fullmodel<-glm(formula = entropy~orderedstage/sex*splicing, data = eigenexon_entropy, family=gaussian(link="log"))
# plot(fullmodel)
# summary(fullmodel)
# model.sel(dredge(fullmodel))
# summary(model.avg(dredge(fullmodel)))
# 
# # Select between different polynomials of the stage variable
# model4NS<-glm(formula = entropy~poly(as.numeric(orderedstage),degree = 4)/splicing, data = eigenexon_entropy)
# model4<-glm(formula = entropy~poly(as.numeric(orderedstage),degree = 4)/sex*splicing, data = eigenexon_entropy)
# model4NI<-glm(formula = entropy~splicing-1/poly(as.numeric(orderedstage),degree = 4)/sex, data = eigenexon_entropy)
# model3<-glm(formula = entropy~poly(as.numeric(orderedstage),degree = 3)/sex*splicing, data = eigenexon_entropy)
# model2<-glm(formula = entropy~poly(as.numeric(orderedstage),degree = 2)/sex*splicing, data = eigenexon_entropy)
# model1<-glm(formula = entropy~poly(as.numeric(orderedstage),degree = 1)/sex*splicing, data = eigenexon_entropy)
# 
# AICc(model1, model2, model3, model4, model4NS)
# 
# summary(model4)
# 
# # Now fit an intercept and polynomial for each splicing category
# nest4<-glm(formula = entropy~-1+splicing*poly(as.numeric(orderedstage),degree = 4)/sex, data = eigenexon_entropy)
# nest3<-glm(formula = entropy~-1+splicing*poly(as.numeric(orderedstage),degree = 3)/sex, data = eigenexon_entropy)
# nest2<-glm(formula = entropy~-1+splicing*poly(as.numeric(orderedstage),degree = 2)/sex, data = eigenexon_entropy)
# nest1<-glm(formula = entropy~-1+splicing*poly(as.numeric(orderedstage),degree = 1)/sex, data = eigenexon_entropy)
# 
# AICc(nest4S1, nest4S, nest4, nest3, nest2, nest1)
# 
# Nest4<-glm(formula = entropy~-1+splicing/(as.numeric(orderedstage)+I(as.numeric(orderedstage)^2)+I(as.numeric(orderedstage)^3)+I(as.numeric(orderedstage)^4)+(as.numeric(orderedstage):sex+I(as.numeric(orderedstage)^2):sex+I(as.numeric(orderedstage)^3):sex+I(as.numeric(orderedstage)^4):sex)), data = eigenexon_entropy)
# dNest<-dredge(Nest4)
# 
# 
# modgenes<-glm(formula = entropy ~ -1+(as.numeric(orderedstage)+I(as.numeric(orderedstage)^2)+I(as.numeric(orderedstage)^3)+I(as.numeric(orderedstage)^4)+(as.numeric(orderedstage):sex+I(as.numeric(orderedstage)^2):sex+I(as.numeric(orderedstage)^3):sex+I(as.numeric(orderedstage)^4):sex)), data = eigenexon_entropy, subset = splicing=="genes")
# 
# modisof<-glm(formula = entropy ~ -1+(as.numeric(orderedstage)+I(as.numeric(orderedstage)^2)+I(as.numeric(orderedstage)^3)+I(as.numeric(orderedstage)^4)+(as.numeric(orderedstage):sex+I(as.numeric(orderedstage)^2):sex+I(as.numeric(orderedstage)^3):sex+I(as.numeric(orderedstage)^4):sex)), data = eigenexon_entropy, subset = splicing!="genes")
# summary(model.avg(dredge(modisof, subset = )))


#### Fit isoforms and genes separately
# # Validate polynomial complexity and interactions with sex in genes
# modgenes<-list(
#   modgenes4S1=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 4)+as.numeric(orderedstage):sex, data = eigenexon_entropy, subset = splicing=="genes"),
#   modgenes4S2=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 4)+poly(as.numeric(orderedstage), degree = 2):sex, data = eigenexon_entropy, subset = splicing=="genes"),
#   modgenes4S3=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 4)+poly(as.numeric(orderedstage), degree = 3):sex, data = eigenexon_entropy, subset = splicing=="genes"),
#   modgenes4=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 4)/sex, data = eigenexon_entropy, subset = splicing=="genes"),
#   modgenes3=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 3)/sex, data = eigenexon_entropy, subset = splicing=="genes"),
#   modgenes2=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 2)/sex, data = eigenexon_entropy, subset = splicing=="genes"),
#   modgenes1=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 1)/sex, data = eigenexon_entropy, subset = splicing=="genes"),
#   modgenesNS4=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 4), data = eigenexon_entropy, subset = splicing=="genes"),
#   modgenesNS3=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 3), data = eigenexon_entropy, subset = splicing=="genes"),
#   modgenesNS2=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 2), data = eigenexon_entropy, subset = splicing=="genes"),
#   modgenesNS1=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 1), data = eigenexon_entropy, subset = splicing=="genes")
# )
# 
# summary(model.avg(modgenes))
# model.sel(modgenes)
# 
# Validate polynomial complexity and interactions with sex in isoforms
# modisof<-list(
#   modisof4S1=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 4)+as.numeric(orderedstage):sex, data = eigenexon_entropy, subset = splicing!="genes"),
#   modisof4S2=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 4)+poly(as.numeric(orderedstage), degree = 2):sex, data = eigenexon_entropy, subset = splicing!="genes"),
#   modisof4S3=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 4)+poly(as.numeric(orderedstage), degree = 3):sex, data = eigenexon_entropy, subset = splicing!="genes"),
#   modisof4=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 4)/sex, data = eigenexon_entropy, subset = splicing!="genes"),
#   modisof3=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 3)/sex, data = eigenexon_entropy, subset = splicing!="genes"),
#   modisof2=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 2)/sex, data = eigenexon_entropy, subset = splicing!="genes"),
#   modisof1=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 1)/sex, data = eigenexon_entropy, subset = splicing!="genes"),
#   modisofNS4=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 4), data = eigenexon_entropy, subset = splicing!="genes"),
#   modisofNS3=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 3), data = eigenexon_entropy, subset = splicing!="genes"),
#   modisofNS2=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 2), data = eigenexon_entropy, subset = splicing!="genes"),
#   modisofNS1=glm(formula = entropy~poly(as.numeric(orderedstage), degree = 1), data = eigenexon_entropy, subset = splicing!="genes")
# )
# 
# summary(model.avg(modisof))
# model.sel(modisof)
# 
# # Plot model predictions
# predictionsG<-matrix(ncol=2, 
#                      predict(modgenes$modgenesNS4, 
#                              newdata = data.frame(orderedstage=rep(1:5,2), 
#                                                   sex=rep(c(rep("male",5),rep("female",5)),1)),
#                              type="response")
# )
# predictionsI<-matrix(ncol=2, 
#                      predict(modisof$modisofNS4, 
#                              newdata = data.frame(orderedstage=rep(1:5,2), 
#                                                   sex=rep(c(rep("male",5),rep("female",5)),1)),
#                              type="response")
# )
# 
# predictions<-cbind(predictionsI, predictionsG)
# 
# colnames(predictions)<-c("male_isoforms", "female_isoforms", "male_genes", "female_genes")
# predictions<-as.data.frame(t(predictions))
# colnames(predictions)<-levels(eigenexon_entropy$orderedstage)
# predictions$splicing<-c(rep("isoforms",2),rep("genes",2))
# predictions$sex<-rep(c("male","female"),2)
# predictions$grouper<-row.names(predictions)
# predictions<-melt(predictions)
# 
# predplot<-ggplot(data = predictions, aes(x=variable, y=value, col=sex, pch=splicing, group=grouper))+geom_point()+geom_line()+theme_bw()+ylab(label = "Entropy (bits)")+xlab(label = "Stage")+ggtitle(label = "Information encoded by expression and splicing nodes \n in each stage from best fitting models")+scale_shape_discrete(name="Node type", labels=c("Expression","Splicing"))+scale_color_discrete(name="Sex", label=c("Female","Male"))+ylim(8000,20000)+theme(legend.position="bottom")#+facet_grid(splicing~.)
# 
# pdf(file=file.path(graphdir, "entropy_TS_models.pdf"), height = 10, width = 11)
# multiplot(dataplot, predplot, cols = 2)
# dev.off()
