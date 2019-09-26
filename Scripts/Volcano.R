## Volcano plots of gene, exon and alternative splicing regulation in sex and development
## Load libraries
library(ggplot2)
library(plyr)
library(reshape)
library(limma)
library(stringr)

setwd(dir="C:/Users/Aldfreo/Documents/My Box Files/Data/Nasonia_development")

## Load data (gene based estimates)
load(file="./Output/BaseStats_fit2.RData")
load(file="./Output/Exon_Basestats_fit2_gene.RData")
load(file="./Output/Splicing_Basestats_fit2_fdrtool.RData")

# Create volcano dataset with genes annotated by fdr for presence of differential expression, presence of >1 spliced exon and >1 DE exon
Volcano_gene<-merge(melt(fit2$fdr), melt(fit2$coefficients), by=c("X1","X2"))
names(Volcano_gene)<-c("geneID","Factor", "fdr", "coef")
Volcano_gene$Set<-"Gene"
## Names of Factor variable do not match. Correct first, then merge
Volcano_splicing<-merge(melt(fit2_splicing$fdr_splicing), melt(fit2_splicing$coefficients_splicing), by=c("X1","X2"))
names(Volcano_splicing)<-c("geneID","Factor", "fdr", "coef")
Volcano_splicing$Set<-"Splicing"
Volcano_exon<-merge(melt(fit2_exon$fdr_exon), melt(fit2_exon$coefficients_exon), by=c("X1","X2"))
names(Volcano_exon)<-c("geneID","Factor", "fdr", "coef")
Volcano_exon$Set<-"Exon"
# Stack into shared data.frame
Volcano<-rbind(Volcano_gene, Volcano_exon, Volcano_splicing)
Volcano$Sex<-as.factor(ifelse(grepl("male", Volcano$Factor), "Sex", "Stage"))
Volcano$Stage<-factor(substr(str_match((Volcano$Factor), "^stage[[:alnum:]]*")[,1], 6,20), levels=c("emb10","emb18","lar51","pupyel","adult"))
Volcano$N<-1

# Plot as dotplot
ggplot(data=Volcano, aes(x=coef, y=fdr))+geom_point(alpha=0.1)+facet_grid(Stage~Sex+Set)+scale_y_log10()+theme_bw()
# Plot with hexagonal binning
ggplot(data=Volcano, aes(x=coef, y=fdr, z=N))+stat_summary_hex(fun=function(z){log10(sum(z))})+facet_grid(Stage~Sex+Set)+scale_y_log10()+theme_bw()

# Plot GENE DE as dotplot, for poster
# H 8.72 x W 7.43, inches
tiff(filename="./Graphics/20140408_Gene_volcano_poster.tiff", width=7.43, height=8.72, units="in", res=500, compression="lzw")
ggplot(data=subset(Volcano, Volcano$Set=="Gene"), aes(x=coef, y=fdr, col=as.ordered((fdr<0.00001)*sign(coef))))+geom_point(alpha=0.5, size=0.7)+facet_grid(Stage~Sex)+scale_y_log10()+theme_bw()+theme(legend.position="none")+xlab(label="log fold change")+ylab(label="local false discovery rate")+ggtitle(label="Differentially Expressed Genes in Sex and Development\n")+scale_color_brewer(type="qual", palette=1)
dev.off()

# Plot GENE DS as dotplot, for poster
# H 8.72 x W 7.43, inches
tiff(filename="./Graphics/20140408_Splicing_volcano_poster.tiff", width=7.43, height=8.72, units="in", res=500, compression="lzw")
ggplot(data=subset(Volcano, Volcano$Set=="Splicing"), aes(x=coef, y=fdr, col=as.ordered((fdr<0.00001)*sign(coef))))+geom_point(alpha=0.5, size=0.7)+facet_grid(Stage~Sex)+scale_y_log10()+theme_bw()+theme(legend.position="none")+xlab(label="log fold change")+ylab(label="local false discovery rate")+ggtitle(label="Differentially Spliced Genes in Sex and Development\n")+scale_color_brewer(type="qual", palette=1)
dev.off()

## Some genes in the gene dataset are absent in the splicing results: which ones?
fit2$coefficients[((row.names(fit2$coefficients)%in%row.names(fit2_splicing$coefficients_splicing))==F),]
## Same for exon results
fit2$coefficients[((row.names(fit2$coefficients)%in%row.names(fit2_exon$coefficients_exon))==F),]
## They are the same
row.names(fit2$coefficients[((row.names(fit2$coefficients)%in%row.names(fit2_splicing$coefficients_splicing))==F),])%in%row.names(fit2$coefficients[((row.names(fit2$coefficients)%in%row.names(fit2_exon$coefficients_exon))==F),])
## Store missing gene names
missinggenes<-row.names(fit2$coefficients[((row.names(fit2$coefficients)%in%row.names(fit2_exon$coefficients_exon))==F),])
write(missinggenes, file="./Output/Genes_missing_from_Exon.txt")
# How many differentially spliced genes are also DE and how many are not?
table(fit2$fdr[which(row.names(fit2$fdr)%in%row.names(fit2_splicing$fdr_splicing)),]<0.00001, fit2_splicing$fdr_splicing<0.00001, dnn=c("DE","DS"))

# 
# par(mfrow=c(2,1))
# # Plot correlation between gene and exon data
# plot(fit2$coefficients[-which(row.names(fit2$coefficients)%in%missinggenes),], fit2_exon$coefficients_exon)
# # Plot exon coefficients vs splicing coefficients (expect no correlation)
# plot(fit2_exon$coefficients_exon, fit2_splicing$coefficients_splicing)
# # Plot exon coefficients vs exon-random values based on gene coefficients (expect retaining correlation)
# plot(c(fit2_exon$coefficients_exon), c(fit2_exon$coefficients_exon)-rnorm(mean=mean(fit2$coefficients), sd=sd(fit2$coefficients), n=length(fit2_exon$coefficients_exon)))

## Re-create them as individual melted data.frames (format gene-name,stage,value (named as dataset) (factor can be obtained by gregexpr on stage))
## Them melt by geneID and stage, calculate factor from gregexpr
genecor<-melt(fit2$coefficients[-which(row.names(fit2$coefficients)%in%missinggenes),])
names(genecor)<-c("geneID","factor","Gene")
exoncor<-melt(fit2_exon$coefficients_exon)
names(exoncor)<-c("geneID","factor","Exon")
splicingcor<-melt(fit2_splicing$coefficients_splicing)
names(splicingcor)<-c("geneID","factor","Splicing")
corplot<-merge(genecor, exoncor, by=c("geneID","factor"))
corplot<-merge(corplot,splicingcor, by=c("geneID","factor"))
corplot$Sex<-grepl("male",corplot$factor)
## Add male intercepts to males
corplot[which(corplot$Sex==T),3:5]<-corplot[which(corplot$Sex==F),3:5]+corplot[which(corplot$Sex==T),3:5]

## Plot estimates set vs set
## Gene vs exon estimates are correlated
ggplot(corplot, aes(x=Gene, y=Exon))+geom_point()+facet_grid(Sex~.)+geom_smooth()
## Gene vs Splicing estimates lose correlation in stage, retain it in sex (sex genewide estimates are poor predictor of individual exon DE)
ggplot(corplot, aes(x=Gene, y=Splicing, col=Exon))+geom_point()+facet_grid(Sex~.)+geom_smooth()
## Exon vs Splicing estimates lose correlation in stage (most DE exon in stage can be predicted by genewide DE), retain it in sex (genewide DE poor predictor of most DE exon)
ggplot(corplot, aes(x=Exon, y=Splicing, col=Gene))+geom_point()+facet_grid(Sex~.)+geom_smooth()

## Same for fdrs
genefdr<-melt(fit2$fdr[-which(row.names(fit2$fdr)%in%missinggenes),])
names(genefdr)<-c("geneID","factor","Gene")
exonfdr<-melt(fit2_exon$fdr_exon)
names(exonfdr)<-c("geneID","factor","Exon")
splicingfdr<-melt(fit2_splicing$fdr_splicing)
names(splicingfdr)<-c("geneID","factor","Splicing")
fdrcor<-merge(genefdr, exonfdr, by=c("geneID","factor"))
fdrcor<-merge(fdrcor,splicingfdr, by=c("geneID","factor"))
fdrcor$Sex<-grepl("male",fdrcor$factor)


# Now plot correlation in different sets
## Genes and Exons are highly correlated (more in stage than in sex)
ggplot(fdrcor, aes(x=Gene, y=Exon, col=log10(Splicing)))+geom_point()+facet_grid(Sex~.)+scale_x_log10()+scale_y_log10()+geom_smooth()
## Genes and splicing lose correlation in stage, retain it in sex (meaning: sex gene-wide intercepts are poor predictor of exon-expression)
ggplot(fdrcor, aes(x=Gene, y=Splicing, col=log10(Splicing)))+geom_point()+facet_grid(Sex~.)+scale_x_log10()+scale_y_log10()+geom_smooth()
## Exon and splicing are not correlated for stages, are still correlated in sex (stage is led by genewide effects with significant intercepts, sex is led by differential splicing so that the behaviour of the most DE exon is predicted by the most differentially spliced exon)
ggplot(fdrcor, aes(x=Exon, y=Splicing, col=log10(Gene)))+geom_point()+facet_grid(Sex~.)+scale_x_log10()+scale_y_log10()+geom_smooth()


# Select only genes with significant fdr!!!
# prediction: correlation between gene and splicing is high with low fdr, correlation between exon and splicing is high with high fdr
fdrcorplot<-corplot
tempfdrcor<-fdrcor
names(tempfdrcor)<-c("geneID","factor","Genefdr","Exonfdr","Splicingfdr")
fdrcorplot<-merge(fdrcorplot, tempfdrcor, by=c("geneID","factor"))
rm(tempfdrcor)

ggplot(fdrcorplot, aes(x=Gene, y=Exon, col=log10(Splicingfdr)))+geom_point()+geom_smooth()+facet_grid(Sex~.)
## Splicing fdr identifies genes with high exon-gene differences, but also other with no difference between most DE exon and gene (are those genes whose expression is skewed by the DE exon?)
ggplot(fdrcorplot, aes(y=Exon-Gene, x=log10(Splicingfdr)))+geom_point()+geom_smooth()+facet_grid(Sex~.)
ggplot(fdrcorplot, aes(y=log10(Exonfdr/Genefdr), x=Splicingfdr))+geom_point()+geom_smooth()+facet_grid(Sex~.)
ggplot(fdrcorplot, aes(y=log10(Exonfdr/Genefdr), x=Splicing))+geom_point()+geom_smooth()+facet_grid(Sex~.)

## Genes with low splicing-gene difference have low fdr for exon DE
ggplot(fdrcorplot, aes(y=Splicing-Gene, col=Exon, x=log10(Splicingfdr)))+geom_point()+geom_smooth()+facet_grid(Sex~.)

ggplot(fdrcorplot, aes(x=log10(Splicingfdr), y=log10(abs((Gene-Splicing)/(Exon-Splicing))), col=Splicing))+geom_point()+geom_smooth()+facet_grid(Sex~.)

## Now for all significant spliced exons
Volcano_splicing2<-merge(melt(fit2_splicing$fdr), melt(fit2_splicing$coefficients), by=c("X1","X2"))
names(Volcano_splicing2)<-c("geneID","Factor", "fdr", "coef")
Volcano_splicing2$Set<-"Splicing2"
Volcano_splicing2$Sex<-as.factor(ifelse(grepl("male", Volcano_splicing2$Factor), "Sex", "Stage"))
Volcano_splicing2$Stage<-factor(substr(str_match(levels((Volcano_splicing2$Factor)), "^stage[[:alnum:]]*")[,1], 6,20), levels=c("emb10","emb18","lar51","pupyel","adult"))

ggplot(data=Volcano_splicing2, aes(x=(2^coef), y=log10(fdr)))+geom_point()+facet_grid(Stage~Sex)
ggplot(data=Volcano_splicing2, aes(x=coef, y=log10(fdr)))+geom_point()+facet_grid(Stage~Sex)