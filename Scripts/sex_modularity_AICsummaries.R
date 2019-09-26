# Initialize script
date()
rm(list=ls())
newdir<-file.path(getwd(), "Output/sex_modularity_AICsummaries")
graphicsdir<-file.path(getwd(), "Graphics/sex_modularity_AICsummaries")
dir.create(newdir)
dir.create(graphicsdir)
library(reshape)
library(ggplot2)
library(MuMIn)
options(na.action="na.fail")

######## NOTE ######
# Sexbias is expressed as number of female samples removed (0-3)
# Negative terms indicate that the parmeter decreases with removal of female samples (is female-specific)
# Positive terms indicate that the parameter increases with removal of female samples (is male-specific)
####################

### Version for dOut (increase in pleiotropy coefficient at the removal of female samples)
# Load data
dOut.avg<-read.csv(file = "./Output/sex_modularity_GLMselectionAIC/dOut_avg.csv")

# Select only rows with non-constitutive factors
table(dOut.avg$factorName, dOut.avg$significant)
dOut.avg<-droplevels(dOut.avg[c(grep(pattern = "sexbias",x = dOut.avg$factorName)),])

# reshape into cluster~stage filling positions with significance (and direction)
dOut.avg$factorName<-factor(dOut.avg$factorName, levels = c("sexbias:stageemb10","sexbias:stageemb18","sexbias:stagelar51","sexbias:stagepupyel","sexbias:stageadult"))
dOut.avg$direction<-ifelse(dOut.avg$significant==T, dOut.avg$Estimate, NA)
dOut_matrix<-recast(data = dOut.avg, formula = clusterID~factorName, measure.var = "direction")

# Order by cluster profile similarity
dOut_ternary<-as.matrix(sign(dOut_matrix[,-1]))
dOut_ternary<-apply(dOut_ternary, c(1,2),function(x){ifelse(is.na(x),0,x)})
row.names(dOut_ternary)<-dOut_matrix[,1]
colnames(dOut_ternary)<-colnames(dOut_matrix[,-1])
dOut_similarity<-princomp(dOut_ternary)$scores[,1]
dOut.avg<-merge(dOut.avg, rank(dOut_similarity, ties.method = "random"), by.x = "clusterID", by.y = 0)
names(dOut.avg)[length(dOut.avg)]<-"similarity"

# How many clusters have sex-bias of dOut in at least one stage?
table(apply(dOut_matrix,1,function(x){any(is.na(x)==F)}))

# How many cluster have multi-stage sex bias?
table(apply(dOut_matrix,1,function(x){sum(is.na(x)==F)}))
prop.table(table(apply(dOut_matrix,1,function(x){sum(is.na(x)==F)})))
plot(table(apply(dOut_matrix,1,function(x){sum(is.na(x)==F)})))

# Which stages have the most convergence based on module parcellation?
princomp_dOut<-princomp(dist(t(dOut_matrix)))
plot(princomp_dOut)
biplot(princomp_dOut)

# Small multiple of module fluctuation in sex bias over time
ggplot(dOut.avg, aes(x=factorName, y=log10(abs(Estimate))*significant*sign(Estimate), group=clusterID))+geom_line()+facet_wrap(~similarity)+ggtitle(label = "Pleiotropy Coefficient")+xlab(label = "stage")+ylab(label = "Sex-Specific Parcellation")+theme_bw()+scale_x_discrete(breaks=NULL)+ theme(strip.background = element_blank(),strip.text.x = element_blank())
ggplot(dOut.avg, aes(x=factorName, y=Estimate*significant, group=clusterID))+geom_line()+facet_wrap(~similarity)+ggtitle(label = "dOut")

# How many clusters are sexually partitioned during embryonal development (y) and adulthood (x)?
table(apply(dOut_matrix[,-c(1,6)],1,function(x){any(is.na(x)==F)}),is.na(dOut_matrix[,6])==F)
prop.table(table(apply(dOut_matrix[,-c(1,6)],1,function(x){any(is.na(x)==F)}),is.na(dOut_matrix[,6])==F))
wilcox.test(table(apply(dOut_matrix[,-c(1,6)],1,function(x){any(is.na(x)==F)}),is.na(dOut_matrix[,6])==F))

# How many clusters display unidirectional parcellation over development (male only or female only)?
directions<-apply(dOut_ternary,1,function(x){
  y<-ifelse(x==0,NA,x)
  length(levels(as.factor(y)))
  })
table(directions)
prop.table(table(directions))

# Histogram of parcellation over developmental stages
ternsummary<-melt(apply(dOut_ternary, 2, table))
names(ternsummary)<-c("bias","stage","count")
ternsummary$stage<-factor(x = ternsummary$stage, levels=c("sexbias:stageemb10","sexbias:stageemb18","sexbias:stagelar51","sexbias:stagepupyel","sexbias:stageadult"))
ggplot(ternsummary, aes(x=stage, y=count, fill=as.factor(bias)))+geom_bar(stat="identity", position="dodge")+theme_bw()

# subset only significant rows

### Version for dWithin (increase in integration at the removal of female samples)
# Load data
dWithin.avg<-read.csv(file = "./Output/sex_modularity_GLMselectionAIC/dWithin_avg.csv")

# Select only rows with non-constitutive factors
table(dWithin.avg$factorName, dWithin.avg$significant)
dWithin.avg<-droplevels(dWithin.avg[c(grep(pattern = "sexbias",x = dWithin.avg$factorName)),])

# reshape into cluster~stage filling positions with significance (and direction)
dWithin.avg$factorName<-factor(dWithin.avg$factorName, levels = c("sexbias:stageemb10","sexbias:stageemb18","sexbias:stagelar51","sexbias:stagepupyel","sexbias:stageadult"))
dWithin.avg$direction<-ifelse(dWithin.avg$significant==T, dWithin.avg$Estimate, NA)
dWithin_matrix<-recast(data = dWithin.avg, formula = clusterID~factorName, measure.var = "direction")

# Order by cluster profile similarity
dWithin_ternary<-as.matrix(sign(dWithin_matrix[,-1]))
dWithin_ternary<-apply(dWithin_ternary, c(1,2),function(x){ifelse(is.na(x),0,x)})
row.names(dWithin_ternary)<-dWithin_matrix[,1]
colnames(dWithin_ternary)<-colnames(dWithin_matrix[,-1])
dWithin_similarity<-princomp(dWithin_ternary)$scores[,1]
dWithin.avg<-merge(dWithin.avg, rank(dWithin_similarity, ties.method = "random"), by.x = "clusterID", by.y = 0)
names(dWithin.avg)[length(dWithin.avg)]<-"similarity"

# How many cluster have multi-stage sex bias?
table(apply(dWithin_matrix,1,function(x){sum(is.na(x)==F)}))
prop.table(table(apply(dWithin_matrix,1,function(x){sum(is.na(x)==F)})))
plot(table(apply(dWithin_matrix,1,function(x){sum(is.na(x)==F)})))

# Which stages have the most convergence based on module parcellation?
princomp_dWithin<-princomp(dist(t(dWithin_matrix)))
plot(princomp_dWithin)
biplot(princomp_dWithin)

# Small multiple of module fluctuation in sex bias over time
ggplot(dWithin.avg, aes(x=factorName, y=log10(abs(Estimate))*significant*sign(Estimate), group=clusterID))+geom_line()+facet_wrap(~similarity)
ggplot(dWithin.avg, aes(x=factorName, y=Estimate*significant, group=clusterID))+geom_line()+facet_wrap(~similarity)+ggtitle(label = "dWithin")

# How many clusters are sexually partitioned during embryonal development (y) and adulthood (x)?
table(apply(dWithin_matrix[,-c(1,6)],1,function(x){any(is.na(x)==F)}),is.na(dWithin_matrix[,6])==F)
prop.table(table(apply(dWithin_matrix[,-c(1,6)],1,function(x){any(is.na(x)==F)}),is.na(dWithin_matrix[,6])==F))
wilcox.test(table(apply(dWithin_matrix[,-c(1,6)],1,function(x){any(is.na(x)==F)}),is.na(dWithin_matrix[,6])==F))

# How many clusters display unidirectional parcellation over development (male only or female only)?
dWithin_directions<-apply(dWithin_ternary,1,function(x){
  y<-ifelse(x==0,NA,x)
  length(levels(as.factor(y)))
})
table(dWithin_directions)
prop.table(table(dWithin_directions))

# Histogram of parcellation over developmental stages
dWithin_ternsummary<-melt(apply(dWithin_ternary, 2, table))
names(dWithin_ternsummary)<-c("bias","stage","count")
dWithin_ternsummary$stage<-factor(x = dWithin_ternsummary$stage, levels=c("sexbias:stageemb10","sexbias:stageemb18","sexbias:stagelar51","sexbias:stagepupyel","sexbias:stageadult"))
ggplot(dWithin_ternsummary, aes(x=stage, y=count, fill=as.factor(bias)))+geom_bar(stat="identity", position="dodge")+theme_bw()

# subset only significant rows
cor(dWithin.avg$Estimate, dOut.avg$Estimate, use = "pairwise", method="pearson")
cor(dWithin.avg$Estimate, dOut.avg$Estimate, use = "pairwise", method="spearman")
cor(dWithin.avg$Estimate, dOut.avg$Estimate, use = "pairwise", method="kendall")

cor(dWithin.avg$Estimate*dWithin.avg$significant, dOut.avg$Estimate*dOut.avg$significant, use = "pairwise", method="pearson")
cor(dWithin.avg$Estimate*dWithin.avg$significant, dOut.avg$Estimate*dOut.avg$significant, use = "pairwise", method="spearman")
cor(dWithin.avg$Estimate*dWithin.avg$significant, dOut.avg$Estimate*dOut.avg$significant, use = "pairwise", method="kendall")

# Compare directions in dOut/dWithin
prop.table(table(dOut_ternary, dWithin_ternary))
# Compare probability of having significant direction in dOut/dWithin
prop.table(table(dOut_ternary!=0, dWithin_ternary!=0))
# Is it statistically significant? (no)
wilcox.test(table(dOut_ternary!=0, dWithin_ternary!=0))

# What is the correlation between increase of pleiotropy and integration once we remove the zeroes?
cor(c(apply(dOut_ternary, c(1,2),function(x){ifelse(x==0,NA,x)})), c(apply(dWithin_ternary, c(1,2),function(x){ifelse(x==0,NA,x)})), use = "pairwise", method="spearman")
# Are the proportions statistically significant? (no)
wilcox.test(table(c(apply(dOut_ternary, c(1,2),function(x){ifelse(x==0,NA,x)})), c(apply(dWithin_ternary, c(1,2),function(x){ifelse(x==0,NA,x)}))))


# How many expression profiles are possible with ternary status (1/0/-1, converted to a/0/b for convenience)?
3^5
# How many are observed (dOut)?
dOut_ternaryprofiles<-apply(dOut_ternary,c(1,2), function(x){
  if (x=="1"){"a"} else if (x=="-1"){"b"} else {x}
})
dOut_ternaryprofiles<-apply(dOut_ternaryprofiles,1,function(x){paste(x,collapse="")})
table(dOut_ternaryprofiles)
length(table(dOut_ternaryprofiles))
# 76 total profiles, 31% of the expected ones
length(table(dOut_ternaryprofiles))/(3^5)
# How many are observed (dWithin)?
dWithin_ternaryprofiles<-apply(dWithin_ternary,c(1,2), function(x){
  if (x=="1"){"a"} else if (x=="-1"){"b"} else {x}
})
dWithin_ternaryprofiles<-apply(dWithin_ternaryprofiles,1,function(x){paste(x,collapse="")})
table(dWithin_ternaryprofiles)
length(table(dWithin_ternaryprofiles))
# 72 total profiles, 30% of the expected ones
length(table(dWithin_ternaryprofiles))/(3^5)

# are there any profiles contained in only one of the two?

# Which profiles are missing from either?
observed_profiles<-union(levels(as.factor(dWithin_ternaryprofiles)),levels(as.factor(dOut_ternaryprofiles)))
expected_profiles<-apply(expand.grid(c("a","b","0"),c("a","b","0"),c("a","b","0"),c("a","b","0"),c("a","b","0")),1,function(x){paste(x, collapse="")})
setdiff(expected_profiles,observed_profiles)
# Could this be because the number of zeroes is inflated (i.e. there is a detection limit)?
# plot distribution of zeroes in observed vs expected patterns
obs_0<-unlist(lapply(strsplit(observed_profiles, split="*"), function(x){sum(x==0)}))
exp_0<-unlist(lapply(strsplit(expected_profiles, split="*"), function(x){sum(x==0)}))
# how many a type unbalances are present in all possible patterns?
exp_a<-unlist(lapply(strsplit(expected_profiles, split="*"), function(x){sum(x=="a")}))
# are there any diverging directions in this pattern?
unbalanced<-unlist(lapply(strsplit(expected_profiles, split="*"), function(x){any(x=="a")&any(x=="b")}))
# how unbalanced are the significant changes in direction?
unbalance<-unlist(lapply(strsplit(expected_profiles, split="*"), function(x){sum(x=="a")-sum(x=="b")}))

# analyze probability of observing profile from number of zeroes (and a) values
profile_probability<-data.frame(profile=expected_profiles, observed=expected_profiles%in%observed_profiles, zeroes=exp_0, a=exp_a, unbalanced=unbalanced, unbalance=unbalance)

profile_probabilityGLM<-glm(formula = observed~zeroes*unbalance*unbalanced, data = profile_probability, family = binomial)
model.sel(dredge(profile_probabilityGLM))
profile_probabilityGLM<-glm(formula = observed~zeroes, data = profile_probability, family = binomial)
# Underrepresented patterns (residuals below 1%)
expected_profiles[which(residuals(profile_probabilityGLM)<quantile(x = residuals(profile_probabilityGLM), probs = 0.01))]
# Overrepresented patterns (residuals above 99%)
expected_profiles[which(residuals(profile_probabilityGLM)>quantile(x = residuals(profile_probabilityGLM), probs = 0.99))]

# analyze number of profiles observed as a zero inflated model (presence and then number of obs)
library(pscl)

# Perform hurdle and negative binomial regression on module frequency for dWithin data only
profile_probability_dWithin<-data.frame(profile=expected_profiles, observed=expected_profiles%in%observed_profiles, zeroes=exp_0, a=exp_a, unbalanced=unbalanced, unbalance=unbalance)
profile_frequency_dWithin<-merge(profile_probability_dWithin, as.data.frame(table(dWithin_ternaryprofiles)), by.x="profile", by.y=1, all.x=T)
profile_frequency_dWithin$Freq<-ifelse(is.na(profile_frequency_dWithin$Freq), 0, profile_frequency_dWithin$Freq)

profile_probability_HURD_dWithin<-hurdle(formula = Freq~zeroes+unbalanced+abs(unbalance)|zeroes+unbalanced+abs(unbalance), data=profile_frequency_dWithin,dist = "negbin", zero.dist = "poisson")
profile_probability_NB_dWithin<-glm.nb(formula = Freq~zeroes+unbalanced+abs(unbalance), data=profile_frequency_dWithin)

model.sel(dredge(profile_probability_NB_dWithin))
model.sel(dredge(profile_probability_HURD_dWithin))

profile_probability_dWithin<-glm.nb(formula = Freq~zeroes, data=profile_frequency_dWithin)
profile_frequency_dWithin$residuals<-residuals(profile_probability_dWithin)

# Perform hurdle and negative binomial regression on module frequency for dOut data only
profile_probability_dOut<-data.frame(profile=expected_profiles, observed=expected_profiles%in%observed_profiles, zeroes=exp_0, a=exp_a, unbalanced=unbalanced, unbalance=unbalance)
profile_frequency_dOut<-merge(profile_probability_dOut, as.data.frame(table(dOut_ternaryprofiles)), by.x="profile", by.y=1, all.x=T)
profile_frequency_dOut$Freq<-ifelse(is.na(profile_frequency_dOut$Freq), 0, profile_frequency_dOut$Freq)

profile_probability_HURD_dOut<-hurdle(formula = Freq~zeroes+unbalanced+abs(unbalance)|zeroes+unbalanced+abs(unbalance), data=profile_frequency_dOut,dist = "negbin", zero.dist = "pois")
profile_probability_NB_dOut<-glm.nb(formula = Freq~zeroes+unbalanced+abs(unbalance), data=profile_frequency_dOut)

model.sel(dredge(profile_probability_NB_dOut))
model.sel(dredge(profile_probability_HURD_dOut))

profile_probability_dOut<-glm.nb(formula = Freq~zeroes+abs(unbalance), data=profile_frequency_dOut)
# Proportion of deviance explained

# Plot dWithin patterns by relative frequency and residual frequency

# Create lookup table for patterns
profiles<-apply(ldply(strsplit(expected_profiles, split = "*")),c(1,2), function(x){ifelse(x=="a",1,ifelse(x=="b",-1,0))})
profiles<-data.frame(expected_profiles, profiles)
names(profiles)<-c("profileID", "emb10", "emb18", "lar51", "pupyel", "adult")
profiles<-recast(profiles, formula = profileID+...~.)
names(profiles)<-c("profileID","stage","status")

# Merge dWithin frequency data with lookup table
dWithinFreqplot<-merge(profile_frequency_dWithin[,c("profile","Freq","residuals")], profiles, by.x="profile", by.y="profileID", all.x=T, all.y=F)
ggplot(data = dWithinFreqplot, aes(x=stage, y=status, group=profile, col=as.factor(sign(residuals))))+geom_line()+facet_grid(residuals~profile)+theme_bw()

##### Outdated code snippets
# ## Hurdle analysis on cumulative dWithin dOut data
# profile_frequency<-merge(profile_probability, as.data.frame(table(cbind(dOut_ternaryprofiles,dWithin_ternaryprofiles))), by.x="profile", by.y=1, all.x=T)
# profile_frequency$Freq<-ifelse(is.na(profile_frequency$Freq), 0, profile_frequency$Freq)
# profile_probabilityHURD<-hurdle(formula = Freq~zeroes+unbalanced|zeroes+unbalanced, data = profile_frequency)
# model.sel(dredge(profile_probabilityHURD))
# summary(model.avg(dredge(profile_probabilityHURD)))
# profile_probabilityHURD<-zeroinfl(formula = Freq~zeroes, data = profile_frequency)
