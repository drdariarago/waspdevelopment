## DRAFT separation in exons vs splicing

# load data
nasdevexon<-read.csv(file="./Output/Nasonia_Dev_Exon.csv")
nasdevsplicing<-read.csv(file="./Output/Nasonia_Dev_Splicing.csv")
nasdevgene<-read.csv(file="./Output/nasoniadevgene.csv")
## Final file will use all available data (including gonads) Must solve problem of splicing in gonads. either include gonads in the complete model (and bias the mean expression levels) or import the intercepts from separate adult-gonad male-female dataset (and have comparability problems)

# reformat into named, transposed matrices
row.names(nasdevexon)<-nasdevexon$X
nasdevexon<-as.matrix(nasdevexon[,-1])
row.names(nasdevsplicing)<-nasdevsplicing$X
nasdevsplicing<-as.matrix(nasdevsplicing[,-1])
row.names(nasdevgene)<-nasdevgene$X
nasdevgene<-as.matrix(nasdevgene[,-1])


## method 1: MDS (2d) between stages or sexes, then visually inspect

# Nonmetric MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

library(MASS)

distexon <- dist(t(nasdevexon)) # euclidean distances between the rows
MDSexon <- isoMDS(distexon, k=4) # k is the number of dim
MDSexon # view results


distsplicing <- dist(t(nasdevsplicing)) # euclidean distances between the rows
MDSsplicing <- isoMDS(distsplicing, k=4) # k is the number of dim
MDSsplicing # view results

distgene <- dist(t(nasdevgene)) # euclidean distances between the rows
MDSgene <- isoMDS(distgene, k=4) # k is the number of dim
MDSgene # view results

save.image(file="./Output/MDS_ExonVsSplicing")

## plot results
Set<-as.factor(c(rep("Exon",nrow(MDSexon$points)),rep("Splicing",nrow(MDSsplicing$points)), rep("Gene",nrow(MDSgene$points))))
MDSNasonia<-as.data.frame(rbind(MDSexon$points, MDSsplicing$points, MDSgene$points), row.names=paste(row.names(rbind(MDSexon$points, MDSsplicing$points, MDSgene$points)), Set, sep="_"))
MDSNasonia$Set<-factor(Set, c("Gene","Exon","Splicing"))
colnames(MDSNasonia)<-c("x", "y","z","w", "Set")
MDSNasonia$Male<-grepl("female", row.names(MDSNasonia))==F
MDSNasonia$Stage<-factor(as.character(regmatches(row.names(MDSNasonia), gregexpr("^[[:alnum:]]+",row.names(MDSNasonia)))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))

## flip y coord of splicing
MDSNasonia[which(MDSNasonia$Set=="Splicing"),c("y")]<-MDSNasonia[which(MDSNasonia$Set=="Splicing"),c("y")]*(-1)

library(ggplot2)
ggplot(MDSNasonia, aes(x=x, y=y, label=row.names(MDSNasonia), col=Male))+geom_text()+facet_grid(Set~.)
ggplot(MDSNasonia, aes(x=x, y=y, col=Male))+geom_point()+facet_grid(Set~Stage)+theme_bw()
ggplot(MDSNasonia, aes(x=x+abs(min(x)), y=y+abs(min(y)), col=Male))+geom_point()+facet_grid(Set~Stage)+theme_bw()+scale_x_log10()+scale_y_log10()

ggplot(MDSNasonia, aes(x=x, y=y, col=Male, shape=Set))+geom_point()+facet_grid(.~Stage)+theme_bw()


## Separate stages by MDS, then compare sexes
ggplot(MDSNasonia, aes(x=x, y=y, col=Stage, shape=Male))+geom_point(size=2)+scale_color_brewer(type="qual", palette=1)+facet_grid(.~Set)+scale_shape_manual(values=c(8,19))+geom_rug()
# ## rotate splicing by 90 degrees
# MDSNasonia2<-MDSNasonia
# MDSNasonia2[which(MDSNasonia2$Set=="Splicing"),"x"]<-MDSNasonia[which(MDSNasonia$Set=="Splicing"),"y"]*(-1)
# MDSNasonia2[which(MDSNasonia2$Set=="Splicing"),"y"]<-MDSNasonia[which(MDSNasonia$Set=="Splicing"),"x"]
# ## re-plot
# ggplot(MDSNasonia2, aes(x=x, y=y, col=Stage, shape=Male))+geom_point(size=10)+scale_color_brewer(type="qual", palette=2)+facet_grid(.~Set)+theme_bw()+scale_shape_manual(values=c("\u2641","\u2642"))
## version for poster
tiff(filename="./Graphics/MDS_poster.tiff", compression="lzw", res=500, height=6.77, width=11.74, units="in")
ggplot(MDSNasonia[which(MDSNasonia$Set%in%c("Exon","Splicing")),], aes(x=x, y=y, col=Stage, shape=Male))+geom_point(size=10)+scale_color_brewer(type="qual", palette=2)+facet_grid(.~Set)+theme_bw()+scale_shape_manual(values=c("\u2641","\u2642"))
dev.off()
# Version for Jack's presentation at Cornell
tiff(filename="./Graphics/MDS_Cornell.tiff", compression="lzw", res=500, height=9, width=16, units="in")
ggplot(MDSNasonia, aes(x=x, y=y, col=Stage, shape=Male))+geom_point(size=10)+scale_color_brewer(type="qual", palette=2)+facet_grid(.~Set)+theme_bw()+scale_shape_manual(values=c("\u2641","\u2642"))+xlab(label="PC1")+ylab(label="PC2")
dev.off()

MDSNasonia3<-MDSNasonia
MDSNasonia3[which(MDSNasonia3$Set=="Splicing"),"x"]<-MDSNasonia[which(MDSNasonia$Set=="Splicing"),"x"]*-1
ggplot(MDSNasonia3, aes(x=x, y=y, col=Stage, shape=Male))+geom_point(size=10)+scale_color_brewer(type="qual", palette=2)+facet_grid(.~Set)+theme_bw()+scale_shape_manual(values=c("\u2641","\u2642"))


## Create summary dataset for confidence intervals
library(plyr)
centroids<-ddply(MDSNasonia2, .(Stage, Male, Set), summarize, y=median(y),x=median(x),ymin=min(y), ymax=max(y), xmin=min(x), xmax=max(x))
## re-plot with confidence interval overlay
ggplot(MDSNasonia2, aes(x=x, y=y, col=Stage, shape=Male))+geom_point(size=2)+scale_color_brewer(type="qual", palette=1)+geom_errorbar(data=centroids, aes(x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax, col="red", size=10))+geom_errorbarh(data=centroids, aes(x=x, y=y, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax, col="red", size=10))+facet_grid(.~Set)+scale_shape_manual(values=c(8,19))

## Separation of adult males-females is parallel to separation between emb-10 and others (ovarian expression). Splicing separation of sexes is perpendicular to adult M-F separation. Adult splicing differences are likely not due to ovarian expression. Males are closer splicing-wise to emb-10 on the emb10-adult axis


## Statistical validation: what is the average (euclidean) distance in the final space between stages and sexes within that stage?

centroids<-ddply(MDSNasonia, .(Stage, Male, Set), summarize, medx=median(x), medy=mean(y))
centroiddist<-ddply(centroids, .(Stage, Male, Set), summarize, distances=dist(x))

## method 2: (semi)-supervised clustering for stages or sexes, then compare performance (semi-supervised: keep stages or sexes together, check separation)
## KODAMA software from http://www.kodama-project.com/, article 10.1073/pnas.1220873111

install.packages(pkgs=c("knnflex","e1071","plsgenomics","sfsmisc"))
source('./Scripts/kodama.r')

### What if the stage separation is on some extra dimension?

MDSexon4 <- isoMDS(distexon, k=4) # k is the number of dim
MDSexon4 # view results

MDSexon4<-as.data.frame(MDSexon4$points)
colnames(MDSexon4)<-c("x", "y", "z","w")
MDSexon4$Male<-grepl("female", row.names(MDSexon4))==F
MDSexon4$Stage<-factor(as.character(regmatches(row.names(MDSexon4), gregexpr("^[[:alnum:]]+",row.names(MDSexon4)))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))

library(rgl)
with(MDSexon4, plot3d(x, y, z, pch=as.numeric(Male)+1, col=as.numeric(Stage), size=3))


MDSsplicing4 <- isoMDS(distsplicing, k=4) # k is the number of dim
MDSsplicing4 # view results

MDSsplicing4<-as.data.frame(MDSsplicing4$points)
colnames(MDSsplicing4)<-c("x", "y", "z","w")
MDSsplicing4$Male<-grepl("female", row.names(MDSsplicing4))==F
MDSsplicing4$Stage<-factor(as.character(regmatches(row.names(MDSsplicing4), gregexpr("^[[:alnum:]]+",row.names(MDSsplicing4)))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))

# First axis: sexes (perpendicular to emb10), second axis: emb10 vs others
ggplot(MDSsplicing4, aes(x=x, y=y, col=Stage, shape=Male))+geom_point(size=10)+scale_color_brewer(type="qual", palette=2)+theme_bw()+scale_shape_manual(values=c("\u2641","\u2642"))
# third axis: stage separation (cannot separate emb18, lar51, separates also adult sexes)
ggplot(MDSsplicing4, aes(x=x, y=z, col=Stage, shape=Male))+geom_point(size=10)+scale_color_brewer(type="qual", palette=2)+theme_bw()+scale_shape_manual(values=c("\u2641","\u2642"))
# fourth axis, gradient across stages, separates sexes in pupae 
ggplot(MDSsplicing4, aes(x=x, y=w, col=Stage, shape=Male))+geom_point(size=10)+scale_color_brewer(type="qual", palette=2)+theme_bw()+scale_shape_manual(values=c("\u2641","\u2642"))

ggplot(MDSsplicing4, aes(x=y, y=w, col=Stage, shape=Male))+geom_point(size=10)+scale_color_brewer(type="qual", palette=2)+theme_bw()+scale_shape_manual(values=c("\u2641","\u2642"))
ggplot(MDSsplicing4, aes(x=z, y=w, col=Stage, shape=Male))+geom_point(size=10)+scale_color_brewer(type="qual", palette=2)+theme_bw()+scale_shape_manual(values=c("\u2641","\u2642"))

MDSsplicing3 <- isoMDS(distsplicing, k=3) # k is the number of dim
MDSsplicing3 # view results

MDSsplicing3<-as.data.frame(MDSsplicing3$points)
colnames(MDSsplicing3)<-c("x", "y", "z")
MDSsplicing3$Male<-grepl("female", row.names(MDSsplicing3))==F
MDSsplicing3$Stage<-factor(as.character(regmatches(row.names(MDSsplicing3), gregexpr("^[[:alnum:]]+",row.names(MDSsplicing3)))),levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))

# First axis: sexes (perpendicular to emb10), second axis: emb10 vs others
ggplot(MDSsplicing3, aes(x=x, y=y, col=Stage, shape=Male))+geom_point(size=10)+scale_color_brewer(type="qual", palette=2)+theme_bw()+scale_shape_manual(values=c("\u2641","\u2642"))

ggplot(MDSsplicing3, aes(x=x, y=z, col=Stage, shape=Male))+geom_point(size=10)+scale_color_brewer(type="qual", palette=2)+theme_bw()+scale_shape_manual(values=c("\u2641","\u2642"))

ggplot(MDSsplicing3, aes(x=y, y=z, col=Stage, shape=Male))+geom_point(size=10)+scale_color_brewer(type="qual", palette=2)+theme_bw()+scale_shape_manual(values=c("\u2641","\u2642"))

library(scatterplot3d)
with(MDSsplicing3, scatterplot3d(x,y,z, main="3D Scatterplot", type="h"))
library(rgl)
with(MDSsplicing3, plot3d(x, y, z, pch=as.numeric(Male)+1, col=as.numeric(Stage), size=3))
with(MDSsplicing4, plot3d(y, z, w, pch=as.numeric(Male)+1, col=as.numeric(Stage), size=3))