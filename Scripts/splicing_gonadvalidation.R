## Cross validation of splicing from gonads

## Initialize script
rm(list=ls())
library(plyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(vcd)
source('./Scripts/multiplot.R', echo=F)
# initialize output path
newdir<-file.path(getwd(), "Output/splicing_gonadvalidation")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/splicing_gonadvalidation")
dir.create(graphdir)

## Load expression data
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/2015_05_19_Justin/nasonia_devtesova.exon.r99.score2")
nasoniaE<-nasonia_devtesova.exon.r99[,c(1,grep("(adult|ovaries|testes).*[1-3]$", names(nasonia_devtesova.exon.r99)))]
rownames(nasoniaE)<-nasoniaE$EXON
nasoniaE<-as.matrix(nasoniaE[,-1])
nasoniaE<-as.data.frame(apply(nasoniaE, c(1,2), function(x){ifelse(x>0,x,0)}))
rm(nasonia_devtesova.exon.r99)

## Load eigenexon assignments
eigenexons_assignments <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexons_assignments.csv")
eigenexons_assignments <- eigenexons_assignments[,c("exonID","eigenexonID","constitutive")]

## Merge exons with respective eigenexon ID
eigenexonsE<-merge(eigenexons_assignments, nasoniaE, by.x="exonID", by.y="row.names", all.x=F, all.y=F)
eigenexonsE_con<-eigenexonsE[which(eigenexonsE$constitutive=="constitutive"),]
eigenexonsE_fac<-eigenexonsE[which(eigenexonsE$constitutive=="facultative"),]

## Average expression in adults and gonads by eigenexons
eigenexonsE_con<-ddply(eigenexonsE_con, .(eigenexonID), numcolwise(median), .progress = "text")
eigenexonsE_fac<-ddply(eigenexonsE_fac, .(eigenexonID), numcolwise(median), .progress = "text")

## Divide spliced nodes by constitutive ones
eigenexonsE_fac2<-apply(eigenexonsE_fac, 1, function(x){
  as.numeric(x[-1])/eigenexonsE_con[grep(str_extract(x["eigenexonID"], "^[[:alnum:]]*."), eigenexonsE_con$eigenexonID),-1]
})
eigenexonsE_fac2<-ldply(eigenexonsE_fac2)
eigenexonsE_fac2<-cbind(eigenexonsE_fac[,1], eigenexonsE_fac2)
names(eigenexonsE_fac2)[1]<-"eigenexonID"

# set infinity scores to NA, set scores greater than 1 to 1
# reasons: higher than 1 result from noise on small-scale measurements, therefore possibly complete signal plus noise
# infinity results in constitutive exons ranking below cutoff in samples where nc exons are above cutoff, therefore possibly spurious signal

eigenexonsE_fac2[,sapply(eigenexonsE_fac2, is.numeric)]<-apply(eigenexonsE_fac2[,sapply(eigenexonsE_fac2, is.numeric)],c(1,2),function(x){as.numeric(ifelse(x==Inf, 0, x))})
eigenexonsE_fac2[,sapply(eigenexonsE_fac2, is.numeric)]<-apply(eigenexonsE_fac2[,sapply(eigenexonsE_fac2, is.numeric)],c(1,2),function(x){as.numeric(ifelse(x>1, 1, x))})

## stack data-frames
eigenexonsE<-rbind(eigenexonsE_con, eigenexonsE_fac2)
eigenexonsE<-eigenexonsE[order(eigenexonsE$eigenexonID),]

## Save as csv
write.csv(eigenexonsE, file=file.path(newdir, "gonadadult_eigenexonsE.csv"))

## Reshape as ID+stage+sex*expr
eigenexonsE_2<-melt(eigenexonsE)
eigenexonsE_2$stage<-ifelse(grepl("male",eigenexonsE_2$variable), "adult", "gonads")
eigenexonsE_2$sex<-ifelse(grepl("female|ovaries",eigenexonsE_2$variable), "female", "male")
eigenexonsE_2$con<-ifelse(grepl("con", eigenexonsE_2$eigenexonID), "constitutive", "facultative")

## Average the triplicates
eigenexonsE_2<-ddply(eigenexonsE_2, .variables = .(eigenexonID, stage, sex, con), summarize, E=median(value), .progress="text")
eigenexonsE_2$stage<-as.factor(eigenexonsE_2$stage)
eigenexonsE_2$sex<-as.factor(eigenexonsE_2$sex)
eigenexonsE_2$con<-as.factor(eigenexonsE_2$con)

## Reshape to enable direct adult/gonad comparisons
eigenexonsE_ratios<-recast(data = eigenexonsE_2, formula = eigenexonID+sex+con~stage)
eigenexonsE_ratios$meanE<-(eigenexonsE_ratios$adult+eigenexonsE_ratios$gonads)/2

## Plot Adult/Gonad expression, split by sex and splicing
ggplot(eigenexonsE_ratios, aes(x=adult, y=gonads))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+geom_point(alpha=0.1)+theme_bw()+geom_smooth(method="lm")

## Compare distributions of gonadal bias
ggplot(eigenexonsE_ratios, aes(x=adult/gonads))+facet_grid(sex~con, scales = "fixed")+geom_density()+theme_bw()+scale_x_log10()

## Compare gonadal bias to average expression
ggplot(na.exclude(eigenexonsE_ratios), aes(x=meanE, y=gonads/adult))+facet_grid(sex~con, scales = "free_x")+geom_density2d()+geom_point(alpha=0.1)+theme_bw()+scale_y_log10()

## Calculate correlations between adults and gonads in m/f con/spl
adult_gonad_correlations<-ddply(eigenexonsE_ratios, .variables = .(sex, con), summarize, correlation_pearson=cor(adult, gonads, use = "pairwise",  method="pearson"), correlation_spearman=cor(adult, gonads, use = "pairwise",  method="spearman"), correlation_kendall=cor(adult, gonads, use = "pairwise",  method="kendall"), .progress = "text")
adult_gonad_correlations

write.csv(adult_gonad_correlations, file=file.path(newdir, "adult_gonad_correlations.csv"))

## Load developmental self-information scores
self.information<- read.csv("./Output/splicing_entropy/eigenexon_self-information.csv")

## Attach self-information to dataset
eigenexonsE_ratios_2<-merge(eigenexonsE_ratios, self.information[,-1], by="eigenexonID", all.x=F, all.y=F)
eigenexonsE_ratios_2$inf_coef_discr<-cut(x = eigenexonsE_ratios_2$inf_coef, breaks = 10)
eigenexonsE_ratios_2$group<-paste(eigenexonsE_ratios_2$sex, eigenexonsE_ratios_2$con, eigenexonsE_ratios_2$inf_coef)

## Plot correlations vs self-information (hypothesis, more specific nodes will correlate better between adults and gonads if the adult signal is gonadal, otherwise vice versa)
ggplot(na.exclude(eigenexonsE_ratios_2), aes(x=adult, y=gonads, alpha=inf_coef))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+geom_point()+theme_bw()

## Same, excluding nodes that are zero in either sample
eigenexonsE_ratios_2$adult1<-ifelse(eigenexonsE_ratios_2$adult>0,1,0)
eigenexonsE_ratios_2$gonad1<-ifelse(eigenexonsE_ratios_2$gonads>0,1,0)

inf_coef_cont1<-ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(x=adult, y=gonads, alpha=inf_coef))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+geom_point()+theme_bw()

ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(x=adult, y=gonads, alpha=0.5, col=inf_coef>0.45))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+geom_point()+theme_bw()

inf_coef_cont2<-ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(x=adult, y=gonads, alpha=inf_coef, col=inf_coef>0.45))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+theme_bw()

pdf(file = file.path(graphdir, "inf_coef_contours2.pdf"), height = 20, width = 20)
multiplot(inf_coef_cont1, inf_coef_cont2)
dev.off()

pdf(file = file.path(graphdir, "inf_coef_scatterplot.pdf"), height = 10, width = 10)
inf_coef_cont1
dev.off()

pdf(file = file.path(graphdir, "inf_coef_contours.pdf"), height = 10, width = 10)
inf_coef_cont2
dev.off()

# Clip to remove nonexpressed in either stage

genbias<-ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2&eigenexonsE_ratios_2$con=="constitutive"),]), aes(x=adult, y=gonads, alpha=inf_coef, col=inf_coef>0.55))+facet_wrap(~sex+con, scales = "free", nrow = 2)+geom_density2d()+theme_bw()+xlim(c(2,6.5))+ylim(c(2,6.5))+geom_smooth(method="lm")+scale_color_brewer(type = "qual", name="Specificity", labels=c("low","high"))+theme(legend.position="bottom")
splbias<-ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2&eigenexonsE_ratios_2$con=="facultative"),]), aes(x=adult, y=gonads, alpha=inf_coef, col=inf_coef>0.55))+facet_wrap(~sex+con, scales = "free", nrow = 2)+geom_density2d()+theme_bw()+xlim(c(0.2,1))+ylim(c(0.2,1))+geom_smooth(method="lm")+scale_color_brewer(type = "qual", name="Specificity", labels=c("low","high"))+theme(legend.position="bottom")

pdf(file=file.path(graphdir, "inf_coef_clipped_contours.pdf"), height = 10, width = 10)
multiplot(genbias, splbias, cols = 2)
dev.off()

# Plot gonadbias vs information coef
ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(y=abs(log10(gonads/adult)), x=inf_coef, col=sex, lty=con))+theme_bw()+scale_x_log10()+geom_boxplot(aes(group=group))+geom_smooth(method="lm")

# Split between gonad and adult biased
eigenexonsE_ratios_2$gonadbias<-(eigenexonsE_ratios_2$gonads/eigenexonsE_ratios_2$adult)<0
ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(y=log10(gonads/adult), x=inf_coef, col=sex, lty=con))+theme_bw()+scale_x_log10()+geom_boxplot(aes(group=group))+geom_smooth(method="lm")+facet_grid(.~gonadbias)

# Tabulate gonad/adult specific nodes
spectable<-table(eigenexonsE_ratios_2$adult1, eigenexonsE_ratios_2$gonad1, eigenexonsE_ratios_2$sex, eigenexonsE_ratios_2$con, dnn=c("adult", "gonad", "sex", "con"))
prop.table(table(eigenexonsE_ratios_2$adult1, eigenexonsE_ratios_2$gonad1, eigenexonsE_ratios_2$sex, eigenexonsE_ratios_2$con, dnn=c("adult", "gonad", "sex", "con")), margin = c(3,4))

# Plot tabulations
pdf(file = file.path(graphdir, "cotabplot_gonadspecific.pdf"))
cotabplot(x = spectable, formula = ~adult+con+gonad+sex, shade = T, panel = cotab_coindep, n=5000)
dev.off()

# Plot differences in self-information across different tabulations
eigenexonsE_ratios_2$express1<-as.factor(paste(eigenexonsE_ratios_2$adult1, eigenexonsE_ratios_2$gonad1))

ggplot(na.exclude(eigenexonsE_ratios_2), aes(x=express1, y=inf_coef))+geom_boxplot(varwidth = T, notch = T)+facet_grid(sex~con)+theme_bw()+scale_x_discrete(labels=c("Gonad \n specific", "Adult \n specific", "Aspecific"), name="")+scale_y_continuous(name="self_information")

eigenexonsE_ratios_2$group1<-paste(eigenexonsE_ratios_2$sex, eigenexonsE_ratios_2$con, eigenexonsE_ratios_2$express1)
ggplot(na.exclude(eigenexonsE_ratios_2), aes(x=group1, y=inf_coef, group=group1, col=sex, lty=con))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+scale_x_discrete(labels=c("Gonad \n specific", "Adult \n specific", "Aspecific"))

pdf(file=file.path(graphdir, "self_inf_gonadspecific.pdf"), width = 15, height = 10)
ggplot(na.exclude(eigenexonsE_ratios_2), aes(x=express1, y=inf_coef))+geom_boxplot(varwidth = T, notch = T)+theme_bw()+scale_x_discrete(labels=c("Gonad \n specific", "Adult \n specific", "Aspecific"), name="")+scale_y_continuous(name="self_information")+facet_wrap(~con+sex, nrow=1)
dev.off()

cotabplot_data<-droplevels(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$express1%in%c("1 1", "0 1", "1 0")),])
cotabplot_data$express1<-mapvalues(cotabplot_data$express1, from = c("1 1", "0 1", "1 0"), to = c("Aspecific","Gonad\nSpecific", "Adult\nSpecific"))
cotabplot_data$inf_coef<-mapvalues(as.factor(cotabplot_data$inf_coef), from = levels(as.factor(cotabplot_data$inf_coef)), to = c(1:10))

pdf(file=file.path(graphdir, "cotabplot_express_selfinf.pdf"), width = 12, height = 12)
cotabplot(~inf_coef+express1|sex+con, data=na.exclude(cotabplot_data), split_vertical=T, shade = T, gp= shading_max, 
          labeling_args = list( rot_labels = c(top=0, right=0),
                                offset_labels = c(right = 0.7),
                                tl_labels=c(T,F),
                                varnames = c(T,F),
                                set_varnames = list(inf_coef="Self-Information Rank"))
)
dev.off()

pdf(file=file.path(graphdir, "cotabplot_express_selfinf_expected.pdf"), width = 12, height = 12)
cotabplot(~inf_coef+express1|sex+con, type="expected", data=na.exclude(cotabplot_data), split_vertical=T, shade = T, gp= shading_max, 
          labeling_args = list( rot_labels = c(top=0, right=0),
                                offset_labels = c(right = 0.7),
                                tl_labels=c(T,F),
                                varnames = c(T,F),
                                set_varnames = list(inf_coef="Self-Information Rank"))
)
dev.off()

## Plot density distributions
ggplot(na.exclude(eigenexonsE_ratios_2), aes(x=adult/gonads, col=inf_coef_discr))+facet_wrap(~sex+con, scales = "fixed")+theme_bw()+geom_density()+scale_x_log10()

## Save workspace
save.image(file = file.path(newdir, "splicing_gonadvalidation.R"))

# # Analyze tabulations as lm
# library(lme4)
# eigenexonsE_ratios_2$randomint<-as.factor(eigenexonsE_ratios_2$inf_coef>0.75)
# lm_infcoef<-glmer(droplevels(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$express1%in%c("0 1","1 0", "1 1")),]), 
#                formula = inf_coef>0.75~-1+sex*con*express1+(1|randomint), family=binomial)
# plot(lm_infcoef)