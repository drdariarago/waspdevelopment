## Find predictors of stage DE

# load data and select an Fdr threshold
load(file="./Output/BaseStats_out.RData")
stagededata<-as.data.frame(apply(fit2$Fdr[,1:5],1, function(x){any(x<0.00001)}))
names(stagededata)<-c("stageDE")
# add OGS2 information
OGS2.0<- read.csv("./Input/NVIT_OGS2_goodannotcomplete.csv")

stage_data<-merge(stagededata,OGS2.0[,c("geneID","LinkageCluster","recombrate","cM","Kb","ODB6_OG_ID","adult_female_meth_status","quality7","isoforms","IntMatch","quality2","ratio")], by.x="row.names", by.y="geneID", all.x=T)
stage_data$LinkageCluster<-as.factor(stage_data$LinkageCluster)
stage_data$Chromosome<-as.factor(substr(stage_data$LinkageCluster, 1, 1))
stage_data$KbCM<-stage_data$Kb/stage_data$cM
names(stage_data)<-c("geneID","stageDE","LinkageCluster","recombrate","cM","Kb","ODB6_OG_ID","Methylated", "Paralog", "Isoforms","IntMatch","quality2", "AA_evol_ratio","Chromosome","KbcM")

# Which factors have the most missing values?
apply(stage_data, 2, function(x){sum(is.na(x))})

# add number of genes on each cluster
stage_data<-merge(stage_data, as.data.frame(table(stage_data$LinkageCluster)), by.x="LinkageCluster", by.y="Var1", all.x=T)
names(stage_data)[ncol(stage_data)]<-"Ngenes"

# Select only models with strong expression support (important when considering isoforms)
stage_data<-droplevels(stage_data[which(stage_data$quality2=="Express:Strong"),])

# Select only models with orthology/paralogy assignments (genes with no ortholog/paralog have too few isoforms)
# mosaicplot(Paralog~stageDE+(Isoforms>1), data=stage_data)
stage_data<-droplevels(stage_data[which(stage_data$Paralog!="None"),])

# How much of OGS2 is present in the final model?
table(OGS2.0$geneID%in%stage_data$geneID)
prop.table(table(OGS2.0$geneID%in%stage_data$geneID))
# After exclusion of entries with NAs?
stage_data<-droplevels(na.omit(stage_data))
table(OGS2.0$geneID%in%stage_data$geneID)
prop.table(table(OGS2.0$geneID%in%stage_data$geneID))


# Fit maximal model
# NOTE: Orthologous group as a random intercept alone is more supported than the rest of the model, must investigate
library(lme4)
binomfull<-glmer(stageDE~Paralog+Methylated+log10(Ngenes)+log10(recombrate+0.000001)+log10(IntMatch+1)+log10(AA_evol_ratio)+(1|LinkageCluster), family="binomial", data=stage_data)

stage_dredge<-dredge(binomfull)
summary(model.avg(stage_dredge))
stage_dredge

binomfull2<-glm(stageDE~ODB6_OG_ID, family="binomial", data=stage_data)

# previous AICc based data-exploration supports retaining OG. Support for Linkage cluster is low if ngenes and recombrate are present

# Calculate bootstrap based parameter estimates
# NOTE: bootMer cannot handle NAs, so remove those before fitting main model (data=na.omit(dataset))

# Parameter sampling function: sample all fixed effect estimates, variance estimates for random effects and sigma (returns error for sigma, so put it last to avoid loop skipping)
mySumm <- function(.) {
  c(beta=fixef(.),sig01=unlist(VarCorr(.)),sigma=sigma(.))
}
# Resample model (specify number of replicates and verbosity)
boo01<-bootMer(binomfull, mySumm, nsim=5, verbose=T, .progress="txt")
# View model
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
pdf(paste(c("./Graphics/",format(Sys.time(), "%y_%m_%d"),"_StageGenes_GLMM_Boot_estimates.pdf"), collapse=""), paper="a4", onefile=T)
for (i in 1:length(boo01$t0)){
  print(plot(boo01, index=i),
        title(main=names(boo01$t0[i])))
}
dev.off()

# Save results to workspace
save.image(paste(c("./Output/",format(Sys.time(), "%y_%m_%d"),"_StageGenes_GLMM_Boot_results.RData"), collapse=""))


## Caveats and model criticism
# 1) Maximal model does not estimate individual linkage groups (will do so if they are a significant variance component in the final model)
# 2) Model relies on OGS2.0 annotation of isoforms, increasing false positive in case of dev/sex specific splicing and restricts the dataset to strong expression support models only (unequal representation in isoforms between different classes of evidence)
# 3) Methylation is based on adult females and is also highly correlated with mean gene expression values and expression support 
# 4) Main sources of NAs result from AA evolution, Methylation and Expression support
# 5) No test for overdispersion/absolute model fit implemented in this version, will run it on minimized model