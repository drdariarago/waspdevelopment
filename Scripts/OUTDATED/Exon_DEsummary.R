# Summary of Exon-based analyses (exon-specific DE rather than gene-specific, based on Exon_BaseStats)
library(reshape)
# Load results
load(file="./Output/Exon_Basestats.RData")

# How many exons are significant for any factor?
sum(apply(Fdrfitall_exon, 1, function(x){any(x<0.00001)}))
# 117061, of which how many estimated false positives?
sum(apply(Fdrfitall_exon, 1, function(x){any(x<0.00001)}))*0.00001
# 1

# Volcano plot (Fdr based, estimates probability of false positives within the whole set of positives)
library(stringr)
library(ggplot2)

Volcano_exon<-merge(melt(fit2_exon$fdr), melt(fit2_exon$coefficients), by=c("X1","X2"))
names(Volcano_exon)<-c("exonID","Factor", "Fdr", "coef")
Volcano_exon$Sex<-ifelse(grepl("male", Volcano_exon$Factor), "Sex", "Stage")
Volcano_exon$Stage<-factor(substr(str_match(levels((Volcano_exon$Factor)), "^stage[[:alnum:]]*")[,1], 6,20), levels=c("emb10","emb18","lar51","pupyel","adult"))


ggplot(Volcano_exon, aes(x=coef, y=Fdr, col=Stage))+geom_point(alpha=0.2)+theme_bw()+geom_hline(aes(yintercept=0.00001))+geom_hline(aes(yintercept=0.0001))+geom_hline(aes(yintercept=0.001))+scale_y_log10()+facet_grid(.~Sex)+scale_color_brewer(type="qual")+guides(colour = guide_legend(override.aes = list(alpha = 1)))


## Statistics on splicing-based gene data
# Dataset for plotting volcano plot (geneID, parameter estimated Sex/Stage, coefficients and Fdr/fdr)
Volcano_gene_exon<-merge(melt(fit2$Fdr_splicing), melt(fit2$coefficients_splicing), by=c("X1","X2"), all=T)
names(Volcano_gene_exon)<-c("geneID","Factor", "Fdr", "coef")
Volcano_gene_exon$Sex<-as.factor(ifelse(grepl("male", Volcano_gene_exon$Factor), "Sex", "Stage"))
Volcano_gene_exon$Stage<-factor(substr(str_match(levels((Volcano_gene_exon$Factor)), "^stage[[:alnum:]]*")[,1], 6,20), levels=c("emb10","emb18","lar51","pupyel","adult"))
Volcano_gene_exon<-merge(Volcano_gene_exon, melt(fit2$fdr_splicing), by.x=c("geneID","Factor"),by.y=c("X1","X2"))
names(Volcano_gene_exon)[grep("value", names(Volcano_gene_exon))]<-"fdr"
Volcano_gene_exon$z<-1

# Fdr based volcano plot
ggplot(Volcano_gene_exon, aes(x=coef, y=Fdr, col=Stage))+geom_point(alpha=0.2)+theme_bw()+geom_hline(aes(yintercept=0.00001))+geom_hline(aes(yintercept=0.0001))+geom_hline(aes(yintercept=0.001))+scale_y_log10()+facet_grid(.~Sex)+scale_color_brewer(type="qual")+guides(col=guide_legend(override.aes= list(alpha = 1)))
# Colorless version with facets for stages
ggplot(Volcano_gene_exon, aes(x=coef, y=Fdr))+geom_point(alpha=0.2)+theme_bw()+geom_hline(aes(yintercept=0.00001))+geom_hline(aes(yintercept=0.0001))+geom_hline(aes(yintercept=0.001))+scale_y_log10()+facet_grid(Stage~Sex)
# fdr based volcano plot (estimates the probability of individual features being false positives)
ggplot(Volcano_gene_exon, aes(x=coef, y=fdr, col=Stage))+geom_point(alpha=0.2)+theme_bw()+geom_hline(aes(yintercept=0.00001))+geom_hline(aes(yintercept=0.0001))+geom_hline(aes(yintercept=0.001))+scale_y_log10()+facet_grid(.~Sex)+scale_color_brewer(type="qual")+guides(col=guide_legend(override.aes= list(alpha = 1)))
# Colorless version with facets for stages
ggplot(Volcano_gene_exon, aes(x=coef, y=fdr))+geom_point(alpha=0.2)+theme_bw()+geom_hline(aes(yintercept=0.00001))+geom_hline(aes(yintercept=0.0001))+geom_hline(aes(yintercept=0.001))+scale_y_log10()+facet_grid(Stage~Sex)
# Hexagonal binning version
ggplot(Volcano_gene_exon, aes(x=coef, y=fdr, z=z))+stat_summary_hex(fun=function(z){log10(sum(z))}, bins=80)+theme_bw()+geom_hline(aes(yintercept=0.00001))+geom_hline(aes(yintercept=0.0001))+geom_hline(aes(yintercept=0.001))+scale_y_log10()+facet_grid(Stage~Sex)
