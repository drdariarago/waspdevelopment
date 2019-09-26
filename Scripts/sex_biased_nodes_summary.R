# Initialize script
## Define sex-biased nodes (LIMMA?)
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
newdir<-file.path(getwd(), "Output/sex_biased_nodes_summary")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/sex_biased_nodes_summary")
dir.create(graphdir)

# Load data
glms2_fdr_coef <- read.csv(file = "./Output/sex_biased_nodes/glms2_fdr_coef.csv")[,-1]
clusterdata <- read.csv(file = "./Output/Results_compiler/clusterdata_full.csv")[,-1]

# Reduce only to genes annotated in set (removes multiple nodes with same pattern within same gene)
glms2_fdr_coef <- glms2_fdr_coef[which(glms2_fdr_coef$node_ID%in%clusterdata$eigenexonID),]
# Set lfdr threshold
threshold<-5^-2
# How many nodes&genes are biased
sum(grepl("f|m", clusterdata$devsexbias))
sum(ddply(.data = clusterdata, .variables = .(geneID), summarize, sexbias2 = any(grepl("f|m",devsexbias)))$sexbias2)
# How many female biased
sum(grepl("f", clusterdata$devsexbias))
sum(ddply(.data = clusterdata, .variables = .(geneID), summarize, sexbias2 = any(grepl("f",devsexbias)))$sexbias2)
# How many male biased
sum(grepl("m", clusterdata$devsexbias))
sum(ddply(.data = clusterdata, .variables = .(geneID), summarize, sexbias2 = any(grepl("m",devsexbias)))$sexbias2)
# How many male and female biased
sum(grepl("m", clusterdata$devsexbias)&grepl("f", clusterdata$devsexbias))
sum(ddply(.data = clusterdata, .variables = .(geneID), summarize, sexbias2 = any(grepl("m",devsexbias))&any(grepl("f",devsexbias)))$sexbias2)

# Plot relative proportions of DE patterns across development
patternfreq<-as.data.frame(table(clusterdata$devsexbias, clusterdata$node_type))
patternfreq<-patternfreq[order(patternfreq$Freq, decreasing = T),]
recast(data = patternfreq, formula = Var1~Var2)

## Quick checking script to assess how many genes contain given pattern
## Patterns of interest:
# single-stage bias "^[.]*[f,m][.]*$"
# bias in >1 stage "[f,m].*[f,m]"
# prepupal only bias "[f,m].*[^f,m]{2,}$"
# prepupal bias "^[.]{,2}[m,f]"
# adut bias "[f,m]$"
# adut only bias "^[^f,m]*[f,m]$"
# any bias "[f,m]"
patternset <- droplevels(clusterdata[grep("[f,m]", clusterdata$devsexbias),])
levels(droplevels(patternset$devsexbias))
length(unique(patternset$geneID))
table(patternset$devsexbias, patternset$node_type)

# How many nodes have sex-bias in any stage?
table(grepl("f|m", clusterdata$devsexbias), clusterdata$node_type)
margin.table(table(grepl("f|m", clusterdata$devsexbias), clusterdata$node_type), 1)
# How many genes have splicing nodes with any sex bias in at least one stage?
length(unique(clusterdata[which(grepl("f|m", clusterdata$devsexbias)&clusterdata$node_type=="splicing"),"geneID"]))
# How many nodes are sexbiased in the pre-pupal stages?
table(grepl("[f,m].{2,}", clusterdata$devsexbias)==T) # 659
659/12612
# How many nodes are sexbiased only in the pre-pupal stages?
table(grepl("[.]{2,}$", clusterdata$devsexbias)&grepl("[f,m]", clusterdata$devsexbias)) # 453
# Proportion vs total sexbiased nodes
453/12612
# Subset nodes with specific pre-pupal sexbias
prepup_only <- clusterdata[which(grepl("[.]{2,}$", clusterdata$devsexbias)&grepl("[f,m]", clusterdata$devsexbias)),]
length(unique(prepup_only$geneID))
# Is thre an interaction with splicing or transcription?
prepup_nodetype <- matrix(c(table(prepup_only$node_type),table(clusterdata$node_type)), ncol=2, dimnames = list(c("Splicing","Transcription"),c("Prepupal","Networkwide")))
prop.table(prepup_nodetype, margin = 2)
fisher.test(prepup_nodetype)
# Is thre an interaction with duplication?
prepup_dupl <- matrix(c(table(prepup_only$Copynumber>1),table(clusterdata$node_type)), ncol=2, dimnames = list(c("Singlecopy","Duplicated"),c("Prepupal","Networkwide")))
prop.table(prepup_dupl, margin = 2)
fisher.test(prepup_dupl)

# How many nodes have female bias in at least one stage?
table(grepl("f", clusterdata$devsexbias), clusterdata$node_type)
margin.table(table(grepl("f", clusterdata$devsexbias), clusterdata$node_type), 1)
# How many genes have splicing nodes with female bias in at least one stage?
length(unique(clusterdata[which(grepl("f", clusterdata$devsexbias)&clusterdata$node_type=="splicing"),"geneID"]))
# How many female biased nodes are not female-biased in adults?
table(grepl("f$", clusterdata$devsexbias)==F&grepl(".*f.", clusterdata$devsexbias)==T) # 53
# Which proportion of the total?
712/7407 # 1%

# How many nodes have male bias in at least one stage?
table(grepl("m", clusterdata$devsexbias), clusterdata$node_type)
margin.table(table(grepl("m", clusterdata$devsexbias), clusterdata$node_type), 1)
# How many genes have splicing nodes with male bias in at least one stage?
length(unique(clusterdata[which(grepl("m", clusterdata$devsexbias)&clusterdata$node_type=="splicing"),"geneID"]))
# How many male biased nodes are not male-biased in adults?
table(grepl("m$", clusterdata$devsexbias)==F&grepl(".*m.", clusterdata$devsexbias)==T) # 1726
# Which proportion of the total?
1726/6111 # 28%
# How many have female adult bias?
table(grepl("f$", clusterdata$devsexbias)==T&grepl(".*m.", clusterdata$devsexbias)==T) # 432
432/6111 # 7%%
# What's the proportion when we remove those?
(1726-432)/6111 # 21%
# Subset nodes with preadult male bias and adult female bias
malepreadult <- clusterdata[which(grepl("m$", clusterdata$devsexbias)==F&grepl(".*m.", clusterdata$devsexbias)==T),]
summary(malepreadult)
prop.table(table(malepreadult$node_type)) # split is almost 50/50, cp 60/40 in unbiased and 40/60 in biased
length(unique(malepreadult$geneID))/nrow(malepreadult) # 20% of malepreadult nodes have another node from the same parent with a similar pattern
levels(droplevels(malepreadult$clusterID))
# Check for proportions of paralogs/isoforms among parent genes
malepreadult_genes <- unique(malepreadult$geneID)
malepreadult_genes <- clusterdata[which(clusterdata$geneID%in%malepreadult_genes),]
malepreadult_genes <- malepreadult_genes[which(malepreadult_genes$node_type=="transcription"),]
table(malepreadult_genes$quality7)
prop.table(table(malepreadult_genes$quality7))
prop.table(table(clusterdata[which(clusterdata$node_type=="transcription"),]$quality7)) #they have MORE paralogs than genomewide, probably not duplication limited
prop.table(table(malepreadult_genes$Copynumber>1))
prop.table(table(clusterdata[which(clusterdata$node_type=="transcription"),]$Copynumber>1)) #they have MORE paralogs than genomewide, probably not duplication limited

## Check for constraints between male and female bias 
# Reduce to genes with more than 1 significant contrast
multibias <- clusterdata[grepl("[f,m].*[f,m]", clusterdata$devsexbias),c("geneID", "eigenexonID", "devsexbias")]
multibias$concordant <- ifelse((grepl("f", multibias$devsexbias)==T&grepl("m", multibias$devsexbias)==F)|(grepl("m", multibias$devsexbias)==T&grepl("f", multibias$devsexbias)==F), "Concordant", "Discordant")

# probability of finding nodes with male sexbias vs probability of finding nodes with female sexbias
table(grepl("m", multibias$devsexbias), grepl("f", multibias$devsexbias), dnn = c("malebias","femalebias"))
prop.table(table(grepl("m", multibias$devsexbias), grepl("f", multibias$devsexbias), dnn = c("malebias","femalebias")))
prop.table(table(grepl("m", multibias$devsexbias), grepl("f", multibias$devsexbias), dnn = c("malebias","femalebias")), 2) # 75% of multi-stage nodes with female bias are also male biased

# Print nodes with rare patterns
rarepatterns <- clusterdata[clusterdata$devsexbias%in%patternfreq[which(patternfreq$Freq<100),1],]
rarepatterns[,c("eigenexonID", "clusterID", "devsexbias")]

# Print early patterns
unique(grep("^[f|m]", x = clusterdata$devsexbias, value = T))
clusterdata[grep("^[f|m]", x = clusterdata$devsexbias), c("eigenexonID","clusterID","devsexbias","dev_expr")]

# Reshape dataset for summary statistics
melt_glms2_fc<-melt(glms2_fdr_coef)
melt_glms2_fc$stage<-factor(as.character(str_extract(melt_glms2_fc$variable, pattern = "^[^_]*")) ,levels=c("emb10", "emb18", "lar51", "pupyel", "adult"))
melt_glms2_fc$sex<-ifelse(grepl("Male", melt_glms2_fc$variable), "Sex", "Stage")
melt_glms2_fc$type<-factor(str_extract(string = melt_glms2_fc$variable, pattern = "[^_]*$"))
melt_glms2_fc$spl<-as.factor(ifelse(grepl("fac", melt_glms2_fc$node_ID), "Splicing", "Transcription"))
glms2_fc<-dcast(data = melt_glms2_fc, formula = node_ID+stage+sex+spl~type)

# tabulate sexbiased nodes per stage and direction (negative are female biased, positive are male biased)
with(glms2_fc[which(glms2_fc$sex=="Sex"),], table(sign(coef), fdr<threshold, stage))

# Plot as volcano
pdf(file = file.path(newdir,"Volcano_plot_dots.pdf"), paper = "a4r")
ggplot(glms2_fc, aes(x=coef, y=fdr, col=spl, alpha=log10(1/fdr)))+geom_point()+geom_hline(aes(yintercept=threshold))+facet_grid(stage~sex+spl)+theme_bw()+scale_y_log10()
dev.off()
# And hexagonally binned volcano
pdf(file = file.path(newdir,"Volcano_plot_hex.pdf"), paper = "a4r")
ggplot(glms2_fc, aes(x=coef, y=fdr))+geom_hex()+geom_hline(aes(yintercept=threshold))+facet_grid(stage~sex+spl)+theme_bw()+scale_y_log10()+scale_fill_gradientn(colours = rainbow(10), trans="log", limits=c(1,1000), breaks = c(1,5,10,25,50,100,250,1000))
dev.off()

### Plot as lines over developmental time
glms_sign_bystage <- ddply(.data = glms2_fc[which(glms2_fc$sex=="Sex"),], .variables = .(stage, sign(coef), spl), summarize, count = sum(fdr<threshold))
names(glms_sign_bystage) <- c("stage","sex","node_type","count")
glms_sign_bystage <- glms_sign_bystage[which(glms_sign_bystage$sex!=0),]
glms_sign_bystage$sex <- as.factor(ifelse(glms_sign_bystage$sex<0, "Female","Male"))
glms_sign_bystage$regtype <- rep(x = "Sex Biased", times=nrow(glms_sign_bystage))
glms_sign_bystage$grouper <- as.factor(apply(X = glms_sign_bystage[,2:3], MARGIN = 1,  paste, collapse="_"))

pdf(file = file.path(graphdir, "SexDE_nodecounts_line.pdf"))
ggplot(data = glms_sign_bystage, aes(x=stage, y=count, col=sex, lty=node_type, group=grouper))+geom_line()+geom_point()+theme_bw()+scale_y_continuous(trans = "log10", breaks=c(5,10,50,100,500,1000,1500,3000), name = "Number of Sex-Biased Nodes", limits = c(1,3500))+scale_x_discrete(name = "Stage", labels = c("Early Embryo","Late Embryo","Larva", "Pupa", "Adult"))+scale_color_discrete(name = "Sex")+scale_linetype_discrete(name = "Node Type")+ggtitle(label = "Number of Sex-Biased Nodes per Stage\n")
dev.off()

## Merge with sexspec and plot as facets side by side
## Plot sex-spec genes per stage expr
spec_bystage <- clusterdata[which(is.na(clusterdata$spec)==F),c("geneID","eigenexonID","spec","node_type","dev_expr")]
spec_stage <- sapply(X = spec_bystage$dev_expr, FUN = str_locate_all, pattern="1")
names(spec_stage) <- spec_bystage$eigenexonID
spec_stage <- ldply(spec_stage)
spec_stage <- spec_stage[,-3]
names(spec_stage) <- c("eigenexonID","stage")
spec_stage$stage <- mapvalues(x = as.factor(spec_stage$stage), from = 1:5, to = c("emb10","emb18","lar51","pupyel","adult"))
spec_bystage <- merge(spec_bystage, spec_stage, by="eigenexonID")
spec_bystage <- ddply(.data = spec_bystage, .variables = .(spec, node_type, stage), summarize, count = length(eigenexonID))
names(spec_bystage) <- c("sex", "node_type", "stage", "count")
spec_bystage$node_type <- mapvalues(x = spec_bystage$node_type, from = c("splicing","transcription"), to = c("Splicing","Transcription"))
spec_bystage$regtype <- rep(x = "Sex Specific", times=nrow(spec_bystage))
spec_bystage$grouper <- as.factor(apply(X = spec_bystage[,1:2], MARGIN = 1,  paste, collapse="_"))

sxsp_sxb_counts<-rbind(spec_bystage[,c("sex","node_type","stage","count","regtype")], glms_sign_bystage[,c("sex","node_type","stage","count","regtype")])
sxsp_sxb_counts$grouper <- as.factor(apply(X = sxsp_sxb_counts[,c("sex","node_type")], MARGIN = 1,  paste, collapse="_"))
sxsp_sxb_counts$node_type <- relevel(x = sxsp_sxb_counts$node_type, ref = "Transcription")

ggplot(data = sxsp_sxb_counts, aes(x=stage, y=count, col=sex, lty=node_type, group = grouper))+geom_line()+geom_point()+theme_bw()+scale_y_continuous(trans = "log10", breaks=c(5,10,50,100,500,1000,1500,3000), name = "Number of Nodes")+scale_x_discrete(name = "Stage", labels = c("Embryo\nEarly","Embryo\nLate","Larva", "Pupa", "Adult"))+scale_color_discrete(name = "Sex")+scale_linetype_discrete(name = "Node Type", labels = c("Gene","Transcript"))+ggtitle(label = "Number of Sex-Biased and Sex-Specific Genes and Isoforms per Stage\n")+facet_grid(.~regtype)
ggsave(filename = file.path(graphdir, "nodecounts_line.pdf"), width = 210, height = 148, units = "mm")

# # Remove sexspec splicing with exactly the same pattern as sexspec transcription
# sexspec_genes <- clusterdata[which(is.na(clusterdata$spec)==F&clusterdata$node_type=="transcription"),c("geneID","dev_expr")]
# sexpec_isoforms <- clusterdata[which(is.na(clusterdata$spec)==F&clusterdata$node_type=="splicing"),c("geneID","eigenexonID","dev_expr")]
# sexspec_isoforms <- merge(sexpec_isoforms, sexpec_genes, by = "geneID", all.x = T)
# 
# length(which(sexspec_isoforms$dev_expr.x==sexspec_isoforms$dev_expr.y))
# sexspec_isoforms[which(sexspec_isoforms$dev_expr.x==sexspec_isoforms$dev_expr.y),]
# length(unique(sexspec_isoforms[which(sexspec_isoforms$dev_expr.x==sexspec_isoforms$dev_expr.y),"geneID"]))
# 
# # Only 140 sex-specific isoforms come from sex-specific genes with the same expression pattern
# 
# sexspec_genes <- clusterdata[which(is.na(clusterdata$spec)==F&clusterdata$node_type=="transcription"),c("geneID")]
# sexpec_isoforms <- clusterdata[which(is.na(clusterdata$spec)==F&clusterdata$node_type=="splicing"),c("geneID","eigenexonID")]
# sexspec_isoforms <- sexspec_isoforms[which(sexspec_isoforms$geneID%in%sexspec_genes==F),]
# summary(sexspec_isoforms)
# dim(sexspec_isoforms)
# 
# # Remove the latter from counts and re-plot
# ## Plot sex-spec genes per stage expr
# spec_bystage <- clusterdata[which(is.na(clusterdata$spec)==F),c("geneID","eigenexonID","spec","node_type","dev_expr")]
# spec_bystage <- spec_bystage[-which(spec_bystage$geneID%in%sexspec_genes&spec_bystage$node_type=="Splicing"),]
# spec_stage <- sapply(X = spec_bystage$dev_expr, FUN = str_locate_all, pattern="1")
# names(spec_stage) <- spec_bystage$eigenexonID
# spec_stage <- ldply(spec_stage)
# spec_stage <- spec_stage[,-3]
# names(spec_stage) <- c("eigenexonID","stage")
# spec_stage$stage <- mapvalues(x = as.factor(spec_stage$stage), from = 1:5, to = c("emb10","emb18","lar51","pupyel","adult"))
# spec_bystage <- merge(spec_bystage, spec_stage, by="eigenexonID")
# spec_bystage <- ddply(.data = spec_bystage, .variables = .(spec, node_type, stage), summarize, count = length(eigenexonID))
# names(spec_bystage) <- c("sex", "node_type", "stage", "count")
# spec_bystage$node_type <- mapvalues(x = spec_bystage$node_type, from = c("splicing","transcription"), to = c("Splicing","Transcription"))
# spec_bystage$regtype <- rep(x = "Sex Specific", times=nrow(spec_bystage))
# spec_bystage$grouper <- as.factor(apply(X = spec_bystage[,1:2], MARGIN = 1,  paste, collapse="_"))
# 
# sxsp_sxb_counts<-rbind(spec_bystage[,c("sex","node_type","stage","count","regtype")], glms_sign_bystage[,c("sex","node_type","stage","count","regtype")])
# sxsp_sxb_counts$grouper <- as.factor(apply(X = sxsp_sxb_counts[,c("sex","node_type")], MARGIN = 1,  paste, collapse="_"))
# sxsp_sxb_counts$node_type <- relevel(x = sxsp_sxb_counts$node_type, ref = "Transcription")
# 
# 
# ggplot(data = sxsp_sxb_counts, aes(x=stage, y=count, col=sex, lty=node_type, group = grouper))+geom_line()+geom_point()+theme_bw()+scale_y_continuous(trans = "log10", breaks=c(5,10,50,100,500,1000,1500,3000), name = "Number of Nodes")+scale_x_discrete(name = "Stage", labels = c("Embryo\nEarly","Embryo\nLate","Larva", "Pupa", "Adult"))+scale_color_discrete(name = "Sex")+scale_linetype_discrete(name = "Node Type", labels = c("Gene","Transcript"))+ggtitle(label = "Number of Sex-Biased and Sex-Specific Genes and Isoforms per Stage\n")+facet_grid(.~regtype)
# ggsave(filename = file.path(graphdir, "nodecounts_line_filtered.pdf"), width = 210, height = 148, units = "mm")


# Save nodes with sexbias in at least 1 stage (and direction) as csv
glms2_sb<-glms2_fc[which(glms2_fc$sex=="Sex"),]
glms2_sb$fdr<-ifelse(glms2_sb$fdr<threshold, "Significant", "Not_Significant")
glms2_sb$coef<-ifelse(glms2_sb$coef>0, "Male_biased", "Female_biased")

# Plot as  histogram
pdf(file = file.path(graphdir, "sexbiased_bystage_histogram.pdf"))
ggplot(glms2_sb[which(glms2_sb$fdr=="Significant"),], aes(x=stage, fill=coef))+geom_histogram(position="dodge")+theme_bw()+scale_y_log10(breaks = rep(1:5,4)*10^(0:3))+facet_wrap(~spl)
dev.off()

# select only sexbiased contrasts
sxbnodes<-glms2_fc[which(glms2_fc$sex=="Sex"&glms2_fc$fdr<threshold),c("node_ID", "stage", "coef", "fdr")]

# How many nodes are sexbiased in development? # currently counting CONTRASTS!
length(unique(sxbnodes$node_ID))
# How many genes are sexbiased in development?
length(unique(str_extract(sxbnodes$node_ID, pattern = "^[^_]*")))
# How many nodes are Female (negative) and Male (positive) biased in each stage?
table(sxbnodes$stage, sign(sxbnodes$coef))
# How many genes are sexbiased before pupation?
length(unique(str_extract(sxbnodes$node_ID[which(sxbnodes$stage%in%c("emb10","emb18","lar51"))], pattern = "^[^_]*")))
# Print pre-pupal sexbiased genes
prepupal_sexbias<-sxbnodes[which(sxbnodes$stage%in%c("emb10","emb18","lar51")),]
prepupal_sexbias$node_type<-as.factor(ifelse(grepl("_fac", prepupal_sexbias$node_ID), "Splicing", "Transcription"))
# And save as csv
write.csv(x = prepupal_sexbias, file = file.path(newdir, "prepupal_sexbias.csv"))

# How many are repeatedly sexbiased?
table(table(sxbnodes$node_ID))
prop.table(table(table(sxbnodes$node_ID)))

# How many are in the same vs opposite direction in different stages?
head(table(sxbnodes$node_ID, sign(sxbnodes$coef)))
table((table(sxbnodes$node_ID, sign(sxbnodes$coef))[,1]>0)+(table(sxbnodes$node_ID, sign(sxbnodes$coef))[,2]>0)>1)
prop.table(table((table(sxbnodes$node_ID, sign(sxbnodes$coef))[,1]>0)+(table(sxbnodes$node_ID, sign(sxbnodes$coef))[,2]>0)>1))
sxbconflict<-which((table(sxbnodes$node_ID, sign(sxbnodes$coef))[,1]>0)+(table(sxbnodes$node_ID, sign(sxbnodes$coef))[,2]>0)>1)

# save list of conflicting nodes as csv file
sxbconflict<-glms2_fdr_coef[which(glms2_fdr_coef$node_ID%in%sxbconflict$node_ID),]
sxbconflict$geneID<-str_extract(string = sxbconflict$node_ID, pattern = "^[^_]*")

length(unique(sxbconflict$node_ID))
length(unique(sxbconflict$geneID))

write.csv(sxbconflict, file=file.path(newdir, "sxbconflict_nodes.csv"))

# Check for nodes with sexbias in pre-pupal stages
sxbconflict[-which(sxbconflict$stage%in%c("adult","pupyel")),]

# Plot distribution of fold changes, split by sex and type
sxbnodes$fac<-ifelse(grepl("fac", sxbnodes$node_ID), "Splicing", "Transcription")
ggplot(sxbnodes, aes(x=abs(coef), fill=as.factor(sign(coef))))+geom_histogram(position="dodge")+theme_bw()+facet_grid(fac~stage)

########## Outdated code
# # Control for enrichment in nodes from duplicated genes and splicing nodes
# library(lme4)
# prepup_clusterdata <- clusterdata[,c("eigenexonID", "geneID", "quality7", "Copynumber", "strata", "node_type", "odb8_og_id")]
# prepup_clusterdata$prepupbias <- as.factor(ifelse(prepup_clusterdata$eigenexonID%in%prepup_only$eigenexonID, "Prepupal_sexbias", "No_Prepupal_sexbias"))
# prepupGLM<- glmer(formula = prepupbias ~ node_type + Copynumber + strata + (1|odb8_og_id), data = prepup_clusterdata, family = binomial(link="logit"))

# Plot one column per sex, sorted by pattern type
# patternfreq_bysex <- patternfreq
# patternfreq_bysex$stages <- unlist(lapply(gregexpr("m|f", patternfreq_bysex$Var1), function(x)paste(x, collapse = "")))
# patternfreq_bysex <- patternfreq_bysex[order(patternfreq_bysex$stages),]
# patternfreq_males <- patternfreq_bysex[grep("m[^f]?$", patternfreq$Var1),]
# patternfreq_females <- patternfreq_bysex[grep("f[^m]?$", patternfreq$Var1),]
# patternfreq_bysex <- merge(patternfreq_females, patternfreq_females, by = "stages", all=T, suffixes = c("_fem","_mal"))
