# Summarize sex biased and sex specific nodes collapsed at gene level for summaries

# Initialize script
date()
# rm(list=setdiff(ls(),"ODB8_EukOGs_genes"))
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
library(limma)
library(fdrtool)
library(MuMIn)
rm(list = ls())
# initialize output path
newdir<-file.path(getwd(), "Output/sex_biased_genes_summary")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/sex_biased_genes_summary")
dir.create(graphdir)

# Load data
glms2_fdr_coef <- read.csv(file = "./Output/sex_biased_nodes/glms2_fdr_coef.csv")[,-1]
clusterdata <- read.csv(file = "./Output/Results_compiler/clusterdata_full.csv")[,-1]

# Summarize at the gene level
# Extract only transcription nodes or "genes"
genedata <- clusterdata[which(clusterdata$node_type=="transcription"), c("geneID","Name","clusterID","odb8_og_id","Copynumber","splicing_cathegory","n_isoforms","devsexbias","sexbias","spec")]
colnames(genedata)[ncol(genedata)]<-"sexspec"

# Summarize splicing nodes by parent gene ID
isoformsdata <- droplevels(clusterdata[which(clusterdata$node_type=="splicing"),])
isoformsdata <- ddply(.data = isoformsdata, .variables = .(geneID), .fun = summarize, .progress = "text",
                      #                       print(as.character(geneID[1])),
                      #                       print(paste(sexbias, collapse = "_")),
                      #                       print(paste(spec, collapse = "_")),
                      sexbias = if (any(grepl("Male",sexbias), na.rm = T)&any(grepl("Female",sexbias), na.rm = T)) "Confbias" else
                        if (any(grepl("Male",sexbias), na.rm = T)) "Malebias" else 
                          if (any(grepl("Female",sexbias), na.rm = T)) "Femalebias" else "Unbiased",
                      sexspec = if (any(grepl("Male",spec), na.rm = T)&any(grepl("Female",spec), na.rm = T)) "Confspec" else
                        if (any(grepl("Male",spec), na.rm = T)) "Malespec" else 
                          if (any(grepl("Female",spec), na.rm = T)) "Femalespec" else "Aspecific"
)
isoformsdata$sexbias <- as.factor(isoformsdata$sexbias)
isoformsdata$sexspec <- as.factor(isoformsdata$sexspec)
summary(isoformsdata)
dim(isoformsdata)

# Add splicing information to parent gene nodes
genedata <- merge(genedata, isoformsdata, by = "geneID", all.x = T, suffixes = c("_transcription","_splicing"))

# QC table
table(genedata$splicing_cathegory, genedata$sexbias_splicing)
table(genedata$splicing_cathegory, genedata$sexspec_splicing)

## Save as csv
write.csv(x = genedata, file = file.path(newdir, "genedata.csv"))

#### Summarize for paper table

# How many genes have at least one form of sex-linked regulation?
sexregulated_genes <- sum(apply(genedata[,9:10], 1, function(x){any(is.na(x)==F)})|genedata$sexbias_splicing%in%c("Confbias","Femalebias","Malebias")|genedata$sexspec_splicing%in%c("Confspec","Femalespec","Malespec"))
sexregulated_genes
# Which proportion of the final network?
sexregulated_genes/nrow(genedata)

## Genes with only one form of sex-linked regulation
sexregulated_genes <- data.frame(
  spectrans = is.na(genedata$sexspec_transcription)==F, 
  biastrans = is.na(genedata$sexbias_transcription)==F, 
  biasspl = genedata$sexbias_splicing%in%c("Confbias","Femalebias","Malebias"), 
  specspl = genedata$sexspec_splicing%in%c("Confspec","Femalespec","Malespec"))
monosexregulated_genes <- genedata[which(apply(sexregulated_genes, 1, sum)==1),]

summary(monosexregulated_genes)
apply(monosexregulated_genes[,9:10], 2, function(x){sum(is.na(x)==F)})
sum(monosexregulated_genes$sexbias_splicing%in%c("Confbias","Femalebias","Malebias"))
sum(monosexregulated_genes$sexspec_splicing%in%c("Confspec","Femalespec","Malespec"))

pdf(file = file.path(graphdir, "venn_regmodes.pdf"), width = 8.27, height = 5.83)
vennDiagram(vennCounts(sexregulated_genes), main = "\n Genes with Sex-Linked Regulation",cex = 1.2,
            names = c("Specific\nTranscription", "Biased\nTranscription", "Biased\nSplicing", "Specific\nSplicing")
            )
dev.off()
# mosaic(table(sexregulated_genes), shade = T)

# Genes with sexbiased transcription
summary(genedata$sexbias_transcription)
prop.table(table(genedata$sexbias_transcription))
table(is.na(genedata$sexbias_transcription))
prop.table(table(is.na(genedata$sexbias_transcription)))

# Genes with sexspec transcription
summary(genedata$sexspec_transcription)
prop.table(table(genedata$sexspec_transcription))
table(is.na(genedata$sexspec_transcription))
prop.table(table(is.na(genedata$sexspec_transcription)))

# Genes with sexbiased splicing
summary(genedata$sexbias_splicing)
prop.table(table(genedata$sexbias_splicing)[-4])
table(genedata$sexbias_splicing%in%c("Confbias","Malebias","Femalebias"))
prop.table(table(genedata$sexbias_splicing%in%c("Confbias","Malebias","Femalebias")))

# Genes with sexspec transcription
summary(genedata$sexspec_splicing)
prop.table(table(genedata$sexspec_splicing)[-1])
table(genedata$sexspec_splicing%in%c("Confspec","Malespec","Femalespec"))
prop.table(table(genedata$sexspec_splicing%in%c("Confspec","Malespec","Femalespec")))

### Subset to genes with no sex-biased or sex-specific transcription

genedata_notrans <- genedata[which(is.na(genedata$sexbias_transcription)&is.na(genedata$sexspec_transcription)),]
dim(genedata_notrans)
ngenes <- nrow(genedata)

# Genes with sexbiased splicing
summary(genedata_notrans$sexbias_splicing)
prop.table(table(genedata_notrans$sexbias_splicing)[-4])
table(genedata_notrans$sexbias_splicing%in%c("Confbias","Malebias","Femalebias"))
table(genedata_notrans$sexbias_splicing%in%c("Confbias","Malebias","Femalebias"))/ngenes

# Genes with sexspec transcription
summary(genedata_notrans$sexspec_splicing)
prop.table(table(genedata_notrans$sexspec_splicing)[-1])
table(genedata_notrans$sexspec_splicing%in%c("Confspec","Malespec","Femalespec"))
table(genedata_notrans$sexspec_splicing%in%c("Confspec","Malespec","Femalespec"))/ngenes

### Subset to genes with no sex-biased or sex-specific splicing

genedata_nospl <- genedata[which((is.na(genedata$sexbias_splicing)|genedata$sexbias_splicing=="Unbiased")&(is.na(genedata$sexspec_splicing)|genedata$sexspec_splicing=="Aspecific")),]
dim(genedata_nospl)

# Genes with sexbiased transcription
summary(genedata_nospl$sexbias_transcription)
prop.table(table(genedata_nospl$sexbias_transcription))
table(is.na(genedata_nospl$sexbias_transcription))
prop.table(table(is.na(genedata_nospl$sexbias_transcription)))

# Genes with sexspec transcription
summary(genedata_nospl$sexspec_transcription)
prop.table(table(genedata_nospl$sexspec_transcription))
table(is.na(genedata_nospl$sexspec_transcription))
prop.table(table(is.na(genedata_nospl$sexspec_transcription)))

##### Total genes with some form of sex bias or splicing
sexreg_gene<-ddply(.data = genedata, .variables = .(geneID), .progress = "text", summarize, 
                   sexreg = 
                     (sexbias_transcription%in%c("Female","Female, Male","Male")|sexspec_transcription%in%c("Female", "Male")|sexbias_splicing%in%c("Confbias", "Femalebias", "Malebias")|sexspec_splicing%in%c("Confspec","Femalespec","Malespec")
                     ),
                   nsexreg = sum(
                     sexbias_transcription%in%c("Female","Female, Male","Male"),sexspec_transcription%in%c("Female", "Male"),sexbias_splicing%in%c("Confbias", "Femalebias", "Malebias"),sexspec_splicing%in%c("Confspec","Femalespec","Malespec")
                   )
)

summary(sexreg_gene)

table(sexreg_gene$sexreg==1)
prop.table(table(sexreg_gene$sexreg==1))

table(sexreg_gene$nsexreg)
prop.table(table(sexreg_gene$nsexreg))

table(sexreg_gene$nsexreg>1)
prop.table(table(sexreg_gene$nsexreg>1))

table(sexreg_gene[which(sexreg_gene$sexreg==1),]$nsexreg>1)
prop.table(table(sexreg_gene[which(sexreg_gene$sexreg==1),]$nsexreg>1))

## Plot distribution of sebiased genes across development
matrix(unlist(strsplit(x = as.character(genedata$devsexbias), split = NULL)), ncol = 5, byrow = F)

### Check for distrib of clusters with different gonadbias
gonadenr_1 <- clusterdata
gonadenr_1$devsexbias <- as.character(gonadenr_1$devsexbias)
gonadenr_1$devsexbias[which(is.na(gonadenr_1$devsexbias))] <- as.character(gonadenr_1$spec[which(is.na(gonadenr_1$devsexbias))])
gonadenr_1$devsexbias <- as.factor(gonadenr_1$devsexbias)
gonadenr <- ddply(.data = gonadenr_1, .variables = .(devsexbias), summarize, 
                  testes = sum(grepl("Testes", gonadbias)), 
                  ovaries = sum(grepl("Ovaries", gonadbias)),
                  gonads = sum(grepl("Gonad", tissuebias)),
                  adult = sum(grepl("Adult", tissuebias)),
                  meth = sum(adult_female_meth_status=="Methylated", na.rm = T),
                  unmeth = sum(adult_female_meth_status=="Unmethylated", na.rm = T)
)
gonadenr$nnodes <- as.numeric(c(table(gonadenr_1$devsexbias), sum(is.na(gonadenr_1$devsexbias))))
gonadenr$pattern <- as.factor(c(levels(gonadenr_1$devsexbias),"NA"))
# gonadenr$meth <- gonadenr$meth/(table(clusterdata$adult_female_meth_status=="Methylated")[2]/nrow(clusterdata))

ggplot(data = gonadenr, mapping = aes(x = nnodes, y=gonads, label=pattern))+geom_text()+scale_x_log10()+scale_y_log10()+theme_bw()
ggplot(data = gonadenr, mapping = aes(x = nnodes, y = gonads/nnodes, label=pattern))+geom_text()+theme_bw()+scale_x_log10()+scale_y_log10()

ggplot(data = gonadenr, mapping = aes(x = nnodes, y = ovaries/nnodes, label=pattern))+geom_text()+theme_bw()+scale_x_log10()+scale_y_log10()

ggplot(data = gonadenr, mapping = aes(x = (gonads/nnodes)+0.1, y = (ovaries/nnodes)+0.1, label=pattern, col = log(nnodes)))+geom_text()+theme_bw()+scale_x_log10()+scale_y_log10()
ggplot(data = gonadenr, mapping = aes(x = (gonads/nnodes)+0.1, y = (testes/nnodes)+0.1, label=pattern, col = log(nnodes)))+geom_text()+theme_bw()+scale_y_log10()+scale_x_log10()

ggplot(data = gonadenr, mapping = aes(x = nnodes, y = meth/nnodes, label=pattern, col = log(nnodes)))+geom_text()+theme_bw()+scale_x_log10()+geom_hline(yintercept=0.35)+geom_smooth(method = "glm", aes(weight = nnodes), family = "binomial")
# Nodes from genes methylated in adult females are also more frequently female-biased in adult females
# Possible dynamic functional correlation between methylation and upregulation in Nasonia

ggplot(data = gonadenr, mapping = aes(x = nnodes, y = unmeth/nnodes, label=pattern, col = log(nnodes)))+geom_text()+theme_bw()+scale_x_log10()+geom_hline(yintercept=0.5)+geom_smooth(method = "glm", aes(weight = nnodes), family = "binomial")
# Male-biased adult patterns are enriched in nodes from genes which are unmethylated in adult females
# Female-biased adult patterns are depleted in nodes from genes which are unmethylated in adult females
# patterns with no adult sex-bias often but not only show expected levels of methylation

ggplot(data = gonadenr[which(gonadenr$unmeth>0&gonadenr$meth>0),], mapping = aes(x = unmeth/nnodes, y = meth/nnodes, label=pattern, col = log(nnodes)))+geom_text()+theme_bw()+scale_x_log10()+geom_smooth(method="glm", aes(weight=nnodes), family="binomial")

### Check for developmental bias of nodes with gonad-biased expression

devgonad <- clusterdata[which(is.na(clusterdata$gonadbias)==F),c("eigenexonID","geneID", "clusterID","devsexbias","dev_expr", "gonadbias")]
table(devgonad$devsexbias, grepl("Ovaries",devgonad$gonadbias))
table(devgonad$devsexbias, grepl("Testes",devgonad$gonadbias))

## Check which clusters are enriched in testes/ovaries genes
prop.table(table(grepl("Testes",devgonad$gonadbias)))
testesclusters <- prop.table(table(devgonad$clusterID, grepl("Testes",devgonad$gonadbias)), margin = 1)==1
testesclusters <- row.names(testesclusters)[which(testesclusters[,2]==T)]
testesclusters
table(droplevels(clusterdata[which(clusterdata$clusterID=="violetred3"),"devsexbias"]))

prop.table(table(grepl("Ovaries",devgonad$gonadbias)))
ovariesclusters <- prop.table(table(devgonad$clusterID, grepl("Ovaries",devgonad$gonadbias)), margin = 1)>0.01
ovariesclusters <- row.names(ovariesclusters)[which(ovariesclusters[,2]==T)]
ovariesclusters
table(droplevels(clusterdata[which(clusterdata$clusterID==""),"devsexbias"]))

## Genes with conflictual splicing between nodes

table(genedata$sexbias_splicing)
matrix(table(genedata$sexbias_splicing), nrow = 2)
fisher.test(matrix(table(genedata$sexbias_splicing), nrow = 2))
mosaicplot(matrix(table(genedata$sexbias_splicing), nrow = 2), shade = T)

table(genedata$sexspec_splicing)[c(2,3,4,1)]
matrix(table(genedata$sexspec_splicing)[c(2,3,4,1)], ncol = 2)
fisher.test(matrix(table(genedata$sexspec_splicing)[c(2,3,4,1)], nrow = 2))
mosaicplot(matrix(table(genedata$sexspec_splicing)[c(2,3,4,1)], nrow = 2), shade = T)

# Nodes with conflictual sex-bias across development

devsexconflict <- table(grepl("m", clusterdata$devsexbias), grepl("f", clusterdata$devsexbias), clusterdata$node_type, dnn = c("male","female","type"))

fisher.test(devsexconflict[,,"transcription"])
mosaicplot(devsexconflict[,,"transcription"], shade = T)

fisher.test(devsexconflict[,,"splicing"])
mosaicplot(devsexconflict[,,"splicing"], shade = T)

# Check match between transcr and spl bias, excluding dev confl genes

noconfl <- genedata[-which(genedata$sexbias_transcription=="Female, Male"),]

malebiasconfl <- table(noconfl$sexbias_transcription=="Male",noconfl$sexbias_splicing=="Femalebias", dnn = c("MaleBiasTrans","FemaleBiasSpl"))
fisher.test(malebiasconfl)
mosaicplot(malebiasconfl, shade = T)

femalebiasconfl <- table(noconfl$sexbias_transcription=="Female",noconfl$sexbias_splicing=="Malebias", dnn = c("FemaleBiasTrans","MaleBiasSpl"))
fisher.test(femalebiasconfl)
mosaicplot(femalebiasconfl, shade = T)

malespecconfl <- table(noconfl$sexbias_transcription=="Male",noconfl$sexspec_splicing=="Femalespec", dnn = c("MaleBiasTrans","FemaleSpecSpl"))
fisher.test(malespecconfl)
mosaicplot(malespecconfl, shade = T)

femalespecconfl <- table(noconfl$sexbias_transcription=="Female",noconfl$sexspec_splicing=="Malespec", dnn = c("FemaleBiasTrans","MaleSpecSpl"))
fisher.test(femalespecconfl)
mosaicplot(femalespecconfl, shade = T)

######## Outdated code snippets
# isoformsdata <- ddply(.data = isoformsdata[1:5000,], .variables = .(geneID), .fun = summarize, 
#                       anybias = any(is.na(sexbias)==F),
#                       malebias = any(sexbias=="Male"),
#                       femalebias = any(sexbias=="Female"),
#                       conflictbias = any(sexbias=="Female")&any(sexbias=="Male"),
#                       anybias = any(is.na(spec)==F),
#                       malespec = any(spec=="Male"),
#                       femalespec = any(spec=="Female"),
#                       conflictspec = any(spec=="Female")&any(spec=="Male"),
#                       progress = "text"
# )