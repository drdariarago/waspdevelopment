## Cross validation of splicing from gonads

## Initialize script
rm(list=ls())
library(plyr)
library(stringr)
library(reshape2)
library(ggplot2)
# initialize output path
newdir<-file.path(getwd(), "Output/splicing_gonadvalidation")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/splicing_gonadvalidation")
dir.create(graphdir)

## Load expression data
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
nasoniaE<-nasonia_devtesova.exon.r99[,c(1,grep("(adult|ovaries|testes).*[1-3]$", names(nasonia_devtesova.exon.r99)))]
rownames(nasoniaE)<-nasoniaE$EXON
nasoniaE<-as.matrix(nasoniaE[,-1])
nasoniaE<-as.data.frame(apply(nasoniaE, c(1,2), function(x){ifelse(x>0,x,0)}))
rm(nasonia_devtesova.exon.r99)

## Load eigenexon assignments
eigenexons_assignments <- read.csv("./Output/exon_to_transcript_clustering_rankingonly/eigenexons_assignments.csv")
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

## Load developmental self-information scores
self.information<- read.csv("./Output/splicing_entropy/eigenexon_self-information.csv")

## Attach self-information to dataset
eigenexonsE_ratios_2<-merge(eigenexonsE_ratios, self.information[,-1], by="eigenexonID", all.x=F, all.y=F)
eigenexonsE_ratios_2$inf_coef_discr<-cut(x = eigenexonsE_ratios_2$inf_coef, breaks = 10)

## Plot correlations vs self-information (hypothesis, more specific nodes will correlate better between adults and gonads if the adult signal is gonadal, otherwise vice versa)
ggplot(na.exclude(eigenexonsE_ratios_2), aes(x=adult, y=gonads, alpha=inf_coef))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+geom_point()+theme_bw()

## Same, excluding nodes that are zero in either sample
eigenexonsE_ratios_2$adult1<-ifelse(eigenexonsE_ratios_2$adult>0,1,0)
eigenexonsE_ratios_2$gonad1<-ifelse(eigenexonsE_ratios_2$gonads>0,1,0)

ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(x=adult, y=gonads, alpha=inf_coef))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+geom_point()+theme_bw()

ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(x=adult, y=gonads, alpha=0.5, col=inf_coef>0.7))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+geom_point()+theme_bw()

ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(x=adult, y=gonads, alpha=inf_coef, col=inf_coef<0.7))+facet_wrap(~sex+con, scales = "free")+geom_density2d()+theme_bw()

# Plot gonadbias vs information coef
eigenexonsE_ratios_2$group<-paste(eigenexonsE_ratios_2$sex, eigenexonsE_ratios_2$con, eigenexonsE_ratios_2$inf_coef)
ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(y=abs(log10(gonads/adult)), x=inf_coef, col=sex, lty=con))+theme_bw()+scale_x_log10()+geom_boxplot(aes(group=group))+geom_smooth(method="lm")

# Split between gonad and adult biased
eigenexonsE_ratios_2$gonadbias<-(eigenexonsE_ratios_2$gonads/eigenexonsE_ratios_2$adult)<1
eigenexonsE_ratios_2$group<-paste(eigenexonsE_ratios_2$sex, eigenexonsE_ratios_2$con, eigenexonsE_ratios_2$inf_coef, eigenexonsE_ratios_2$gonadbias)

ggplot(na.exclude(eigenexonsE_ratios_2[which(eigenexonsE_ratios_2$adult1+eigenexonsE_ratios_2$gonad1==2),]), aes(y=abs(log10(gonads/adult)), x=inf_coef, col=sex, lty=con))+theme_bw()+scale_x_log10()+geom_boxplot(aes(group=group))+geom_smooth(aes(size=gonadbias),method="lm")+scale_size_discrete(range = c(1,1.5))

# Tabulate gonad/adult specific nodes
spectable<-table(eigenexonsE_ratios_2$adult1, eigenexonsE_ratios_2$gonad1, eigenexonsE_ratios_2$sex, eigenexonsE_ratios_2$con, dnn=c("adult", "gonad", "sex", "con"))
prop.table(table(eigenexonsE_ratios_2$adult1, eigenexonsE_ratios_2$gonad1, eigenexonsE_ratios_2$sex, eigenexonsE_ratios_2$con, dnn=c("adult", "gonad", "sex", "con")), margin = c(3,4))

# Plot tabulations
library(vcd)
cotabplot(x = spectable, formula = ~adult+con+gonad+sex, shade = T, panel = cotab_coindep, n=5000)

## Plot density distributions
ggplot(na.exclude(eigenexonsE_ratios_2), aes(x=adult/gonads, col=inf_coef_discr))+facet_wrap(~sex+con, scales = "fixed")+theme_bw()+geom_density()+scale_x_log10()
