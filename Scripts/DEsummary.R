# Load data from BaseStats
load(file="./Output/BaseStats_fit2.RData")

# load packages
library(ggplot2)
library(reshape)
library(stringr)

# How many genes significant for any factor?
sum(apply(fit2$fdr, 1, function(x){any(x<0.00001)}))
# 20674, of which how many estimated false positives?
sum(apply(fit2$fdr, 1, function(x){any(x<0.00001)}))*0.00001
# 0.2

# how many genes are significant for each factor?
FactorSummary<-as.data.frame(apply(fit2$fdr, 2, function(x){sum(x<0.00001)}))
FactorSummary<-cbind(row.names(FactorSummary), FactorSummary, as.data.frame(apply(fit2$fdr, 2, function(x){sum(x<0.0001)})), as.data.frame(apply(fit2$fdr, 2, function(x){sum(x<0.001)})))
names(FactorSummary)<-c("Factor","q-value cutoff E-5","q-value cutoff E-4","q-value cutoff E-3")
# table summary
FactorSummary
# and plot (sex vs stage DE genes)
ggplot(melt(fit2$fdr), aes(x=value, fill=grepl("male",X2)))+geom_histogram(position="dodge")+scale_x_log10()+theme_bw()
# plot all factors
ggplot(melt(fit2$fdr), aes(x=value, fill=X2))+geom_histogram(position="dodge")+scale_x_log10()+theme_bw()+scale_fill_brewer(type="div", palette=2)

# tabulate stage-specific and sex-specific genes
stspec<-apply(fit2$fdr, 1, function(x){any(x[1:5]<0.00001)})
sxspec<-apply(fit2$fdr, 1, function(x){any(x[6:10]<0.00001)})
addmargins(table(stspec, sxspec),margin=c(1,2))
prop.table(table(stspec, sxspec))
prop.table(table(stspec, sxspec), margin=1)
# is there any interaction betwen the two?
fisher.test(table(stspec, sxspec))
# yes, stage-specific genes hage 10 times more sex-specific genes than the others. p-value~0
prop.table(table(stspec, sxspec), margin=1)
# and only 34 of all sex-specific genes (0.01%) are not DE in any stage as well
addmargins(table(stspec, sxspec), margin=1)

# How many genes are DE in any number of stages?
table(apply(fit2$fdr[,1:5],1,function(x) sum(x<0.00001)))

# How many genes are sex DE in any number of stages?
table(apply(fit2$fdr[,6:10],1,function(x) sum(x<0.00001)))

# tabulate single-stage specific expression for sex
stspec1<-apply(fit2$fdr, 1, function(x){sum(x[1:5]<0.00001)==1})
sxspec1<-apply(fit2$fdr, 1, function(x){sum(x[6:10]<0.00001)==1})
addmargins(table(stspec1, sxspec1),margin=c(1,2))
prop.table(table(stspec1, sxspec1))
prop.table(table(stspec1, sxspec1), margin=1)

# How many stage specific genes are also single stage specific?
addmargins(table(stspec, stspec1),margin=c(2))
prop.table(table(stspec, stspec1),margin=c(2))
# How many sex specific genes are also single-stage specific?
addmargins(table(sxspec, sxspec1),margin=c(2))
prop.table(table(sxspec, sxspec1),margin=c(2))


### Comparison with Ranz 2003 (Drosophila, estimated FDR~0.03%)
sxspec<-apply(Fdrfitall, 1, function(x){any(x[6:10]<0.03)})
prop.table(table(sxspec)) # 45% in Nasonia, ~70% in Drosophila
sum(sxspec==T)*0.03 # 312 expected false positives in Nasonia CI=44-46%

### Comparison with Graveley 2011 (0.1% FDR, Drosophila adults)
sxspec<-Fdrfitall[,10]<0.1
prop.table(table(sxspec)) # 42% in Nasonia, ~29% in Drosophila adults, increases to 73% when considering all stages
sum(sxspec==T)*0.1 # 975 expected false positives in Nasonia 
(sum(sxspec==T)-975)/length(sxspec) # Adults CI=38-46%, Development CI=65-80%
# Also report more male biased (1829) than female biased (945), checking for convergence
sxspec<-as.factor((fit2$Fdr[,10]<0.01)*sign(fit2$coefficients[,10]))
table(sxspec) # sex-specific genes are balanced, males predominate with 456 genes extra (less than the expected number of false positives)
sxspec<-matrix(c(2647,2296,945,1829), ncol=2, dimnames=list(Sex=c("Female", "Male"), Species=c("Nasonia","Drosophila")))
fisher.test(sxspec) # significantly different (p-value~0)
# Increase stringency of FDR to decrease noise (<1 false positive expected)
sxspec<-as.factor((fit2$Fdr[,10]<0.0001)*sign(fit2$coefficients[,10]))
table(sxspec) # sex-specific genes are less balanced, males predominate with 583 genes extra (~1 false positive expected in the whole set)
sxspec<-matrix(c(1102,1685,945,1829), ncol=2, dimnames=list(Sex=c("Female", "Male"), Species=c("Nasonia","Drosophila")))
fisher.test(sxspec) # significance in differences decreases but still highly significant (P-value<0.0001)

# Graveley mentions a role for transcriptional complexity of testes, validating this hypothesis with gonad data
sxspec<-as.factor((gonadfit2$Fdr[,4]<0.01)*sign(gonadfit2$coefficients[,4]))
table(sxspec) # sex-specific genes are balanced, males predominate with 98 genes extra (203 expected false positives)
sxspec<-matrix(c(2860,2958,945,1829), ncol=2, dimnames=list(Sex=c("Female", "Male"), Species=c("Nasonia","Drosophila")))
fisher.test(sxspec) # significantly different (p-value~0)
# Increase stringecy to <1 False positive
sxspec<-as.factor((gonadfit2$Fdr[,4]<0.0001)*sign(gonadfit2$coefficients[,4]))
table(sxspec) # Males predominate with 297 genes extra (~1 false positive expected in the set)
sxspec<-matrix(c(1380,1677,945,1829), ncol=2, dimnames=list(Sex=c("Female", "Male"), Species=c("Nasonia","Drosophila")))
fisher.test(sxspec) # significantly different (p-value~0)

### checking wiht the exon dataset would provide an estimate of splicing differentiation vs gene expression differentiation
# Importing Splicing data
load(file="./Output/Splicing_Basestats_fit2_gene.RData")

# How many genes significant for any factor?
sum(apply(fit2_splicing$fdr_splicing, 1, function(x){any(x<0.00001)}))
# 14110, of which how many estimated false positives?
sum(apply(fit2_splicing$fdr_splicing, 1, function(x){any(x<0.00001)}))*0.00001
# 0.1

# how many genes are significant for each factor?
FactorSummary_Splicing<-as.data.frame(apply(fit2_splicing$fdr_splicing, 2, function(x){sum(x<0.00001)}))
FactorSummary_Splicing<-cbind(row.names(FactorSummary_Splicing), FactorSummary_Splicing, as.data.frame(apply(fit2_splicing$fdr_splicing, 2, function(x){sum(x<0.0001)})), as.data.frame(apply(fit2_splicing$fdr_splicing, 2, function(x){sum(x<0.001)})))
names(FactorSummary_Splicing)<-c("Factor","q-value cutoff E-5","q-value cutoff E-4","q-value cutoff E-3")
# table summary
FactorSummary_Splicing
# and plot (sex vs stage DE genes)
ggplot(melt(fit2_splicing$fdr_splicing), aes(x=value, fill=grepl("male",X2)))+geom_histogram(position="dodge")+scale_x_log10()+theme_bw()
# plot all factors
ggplot(melt(fit2_splicing$fdr_splicing), aes(x=value, fill=X2))+geom_histogram(position="dodge")+scale_x_log10()+theme_bw()+scale_fill_brewer(type="div", palette=2)

# tabulate stage-specific and sex-specific genes
stspec_splicing<-apply(fit2_splicing$fdr_splicing, 1, function(x){any(x[1:5]<0.00001)})
sxspec_splicing<-apply(fit2_splicing$fdr_splicing, 1, function(x){any(x[6:10]<0.00001)})
addmargins(table(stspec_splicing, sxspec_splicing),margin=c(1,2))
prop.table(table(stspec_splicing, sxspec_splicing))
prop.table(table(stspec_splicing, sxspec_splicing), margin=1)
# is there any interaction betwen the two?
fisher.test(table(stspec_splicing, sxspec_splicing))
# yes, stage-specific genes have 7 times more sex-specific genes than the others. p-value~0
prop.table(table(stspec_splicing, sxspec_splicing), margin=1)
# and only 423 of all sex-specific genes (0.10%) are not DE in any stage as well
addmargins(table(stspec_splicing, sxspec_splicing), margin=1)

# How many genes are DS in any number of stages?
table(apply(fit2_splicing$fdr_splicing[,1:5],1,function(x) sum(x<0.00001)))

# How many genes are sex DS in any number of stages?
table(apply(fit2_splicing$fdr_splicing[,6:10],1,function(x) sum(x<0.00001)))

## Plot significantly DE vs DS genes at fdr 10-5
FactorSummary2<-as.data.frame(rbind((fit2$fdr<0.00001)[which(row.names(fit2$fdr)%in%row.names(fit2_splicing$fdr_splicing)),], fit2_splicing$fdr_splicing<0.00001))
FactorSummary2$Set<-as.factor(c(rep("DE", nrow((fit2$fdr<0.00001)[which(row.names(fit2$fdr)%in%row.names(fit2_splicing$fdr_splicing)),])), rep("DS", nrow(fit2_splicing$fdr_splicing<0.00001))))
FactorSummary2$geneID<-row.names(FactorSummary2)
FactorSummary2<-melt(FactorSummary2, id.vars=c("geneID","Set"))
names(FactorSummary2)
FactorSummary2$Sex<-as.factor(ifelse(grepl("male", FactorSummary2$variable), "Sex","Stage"))
FactorSummary2$Stage<-factor(substr(str_match(levels((FactorSummary2$variable)), "^stage[[:alnum:]]*")[,1], 6,20), levels=c("emb10","emb18","lar51","pupyel","adult"))

ggplot(data=FactorSummary2[which(FactorSummary2$value==T),], aes(x=Stage, fill=Set))+geom_histogram(position="dodge")+facet_grid(.~Sex)+scale_y_log10()+theme_bw()+ylab(label="Number of significant genes with a fdr of 10E-5")
# Stage-specific differences are more pronounced in DE genes, while Sex-specific differences are most evident in alternatively spliced genes
# DE genes are always more than DS genes in stages, DS genes are always more than DE genes between sexes

## Same for single-stage genes (are genes that drive differences between sexes DE or DS)?
## Are genes that drive sex differences between stages DE or DS?
FactorSummaryDEstage<-apply(fit2$fdr[,1:5], 1, function(x){sum(x<0.00001)==1})
FactorSummaryDEsex<-na.exclude(apply(fit2$fdr[,6:10], 1, function(x){ifelse(sum(x<0.00001)==1, which(x<0.00001), NA)}))
FactorSummaryDSstage<-na.exclude(apply(fit2_splicing$fdr_splicing[,1:5], 1, function(x){ifelse(sum(x<0.00001)==1, which(x<0.00001), NA)}))
FactorSummaryDSsex<-na.exclude(apply(fit2_splicing$fdr_splicing[,6:10], 1, function(x){ifelse(sum(x<0.00001)==1, which(x<0.00001), NA)}))

FactorSummary3<-dlply(FactorSummary2, .(Set, Sex))
FactorSummary3$DE.Sex<-FactorSummary3$DE.Sex[which(FactorSummary3$DE.Sex$geneID%in%names(FactorSummaryDEsex)),]
FactorSummary3$DS.Sex<-FactorSummary3$DS.Sex[which(FactorSummary3$DS.Sex$geneID%in%names(FactorSummaryDSsex)),]
FactorSummary3$DE.Stage<-FactorSummary3$DE.Stage[which(FactorSummary3$DE.Stage$geneID%in%names(FactorSummaryDEstage)),]
FactorSummary3$DS.Stage<-FactorSummary3$DS.Stage[which(FactorSummary3$DS.Stage$geneID%in%names(FactorSummaryDSstage)),]
FactorSummary3<-ldply(FactorSummary3)

FactorSummary3<-FactorSummary3[which(FactorSummary3$value==T),]

ggplot(FactorSummary3, aes(x=Stage, fill=Set))+geom_histogram(position="dodge")+facet_grid(.~Sex)+theme_bw()

## Some problems in picking the genes that are sex/stage specifically expressed/DE. Perhaps a multivariate approach would be more fruitul, i.e. PCA, then check the rank of PCs that separate sex and stage respectively, or clustering with comparison of sex/stage separation

prin_gene<-princomp(as.matrix(fit2$fdr))
plot(prin_gene)
summary(prin_gene)
prin_gene$loadings
biplot(prin_gene, display=c("sites"))

prin_splicing<-princomp(as.matrix(fit2_splicing$fdr))
plot(prin_splicing)
summary(prin_splicing)
prin_splicing$loadings

## PCA on gene values is unable to resolve sex-differences between sexes. PCA on isoforms does that and increases the ability to separate stages (must test statistically), check the Vegan library for PCA measures, some factors are devoid of loadings:why? Now addressing the issue from raw data using the adonis script. Also check cca function (still vegan)

## now on comlete splicing set
prin_splicing2<-princomp(as.matrix(fit2_splicing$qr))
plot(prin_splicing2)
summary(prin_splicing2)
prin_splicing2$loadings