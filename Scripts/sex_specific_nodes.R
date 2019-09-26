## Annotate sex-specific and sex-biased nodes
# Initialize script
date()
rm(list=ls())
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
# initialize output path
newdir<-file.path(getwd(), "Output/sex_specific_nodes")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/sex_specific_nodes")
dir.create(graphdir)

# load dataset
eigenexon_evalues <- read.csv("./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexon_evalues.csv")

## Version considering expression only if expressed in at least 2 replicates in the same condition
m_eigenexon_evalues<-melt(eigenexon_evalues)
m_eigenexon_evalues$sex<-as.factor(ifelse(grepl("female", m_eigenexon_evalues$variable), "Female", "Male"))
m_eigenexon_evalues$stage<-as.factor(str_extract(m_eigenexon_evalues$variable, pattern = "^[^_]*"))
m_eigenexon_evalues$stage<-factor(m_eigenexon_evalues$stage, levels = c("emb10", "emb18", "lar51", "pupyel", "adult"))

sexspec_rel<-ddply(m_eigenexon_evalues, .(eigenexonID, stage, sex), summarize, expressed=sum(value>0)>1, .progress = "text")

sexspec_rel<-recast(sexspec_rel, eigenexonID+stage~sex)
sexspec_rel<-
  data.frame(eigenexonID=sexspec_rel$eigenexonID,
             geneID=str_extract(string = sexspec_rel$eigenexonID, pattern = "^[^_]*"),
             stage=sexspec_rel$stage,
             spec=revalue(x = as.factor(sexspec_rel$Female-sexspec_rel$Male), replace = c("1"="Female", "-1"="Male", "0"="Aspecific")),
             expr=as.factor((sexspec_rel$Female+sexspec_rel$Male)>0)
  )
sexspec_rel<-sexspec_rel[order(sexspec_rel$geneID, sexspec_rel$eigenexonID, sexspec_rel$stage),]

## Search for genes that are always sex-specific when expressed
sexspec_allstages <- ddply(.data = sexspec_rel[which(sexspec_rel$expr==T),], .variables = .(geneID, eigenexonID), .progress = "text",
                           .fun = summarize,
                           spec = as.factor(if (all(spec=="Female")) "Female" else if (all(spec=="Male")) "Male" else NA)
                           )
# Write to file for Reslts compiler and removal from DE gene analyses
write.csv(sexspec_allstages, file=file.path(newdir, "sexspec_allstages.csv"))

## summarize, all_males, all_females, unbiased, then collapse
sexspec_rel2<-ddply(sexspec_rel, .(geneID, eigenexonID), summarize, 
                    stage_m = paste(stage[which(spec=="Male")], collapse=", "),
                    stage_f = paste(stage[which(spec=="Female")], collapse=", "),
                    spec = as.factor(if (all(spec!="Female")&any(spec=="Male")) "Male" else if (all(spec!="Male")&any(spec=="Female")) "Female" else "Aspecific"),
                    expr=any(expr=T),
                    .progress="text")

# collapse stage_m and stage_f into single column stage_expr
sexspec_rel2$stage_expr<-ifelse(sexspec_rel2$spec!="Aspecific", paste(sexspec_rel2$stage_m, sexspec_rel2$stage_f, sep=""), NA)
sexspec_rel2<-sexspec_rel2[,setdiff(names(sexspec_rel2), c("stage_m","stage_f","expr"))]


### Detect stages in which nodes are expressed, select only stages in which there is expressin (>2 of 3) in males OR females
stagexpr<-ddply(.data = sexspec_rel, .variables = .(eigenexonID, stage), summarize, expressed=as.factor(ifelse(any(expr==T),1,".")), .progress = "text")
# Collapse into short format of zeroes and ones
stagexpr <- recast(data = stagexpr, formula = eigenexonID ~ stage)
stagexpr <- cbind(as.character(stagexpr[,1]), apply(X = stagexpr[,-1], MARGIN = 1, paste, collapse=""))
stagexpr <- as.data.frame(stagexpr)
names(stagexpr) <- c("eigenexonID", "dev_expr")
stagexpr$dev_expr <- as.character(stagexpr$dev.expr)
# Save as csv
write.csv(x = stagexpr, file = file.path(newdir, "stagexpr_nodes_short.csv"))

# Which number of nodes is respectively male and female specific?
summary(sexspec_rel2$spec)
## from how many genes?
# male genes
length(unique(sexspec_rel2$geneID[which(sexspec_rel2$spec=="Male")]))
# female genes
length(unique(sexspec_rel2$geneID[which(sexspec_rel2$spec=="Female")]))

# Save sex-specific nodes annotated with genes
write.csv(sexspec_rel2, file=file.path(newdir, "sexspec_bystage.csv"))

# Save list of female specific genes
write(as.character(unique(sexspec_rel2$geneID[which(sexspec_rel2$spec=="Female")])), file = file.path(newdir, "Female_spec.txt"))
# Save list of male specific genes
write(as.character(unique(sexspec_rel2$geneID[which(sexspec_rel2$spec=="Male")])), file = file.path(newdir, "Male_spec.txt"))

# Which genes possess both male and female specific nodes (splicing mediated sexual conflict)?
spec_levels<-ddply(sexspec_rel2, .(geneID), summarize, 
                   spec_levels=length(levels(droplevels(as.factor(spec)))),
                   any_aspec=any(spec=="Aspecific"),
                   .progress="text")

spec_levels$conflicting<-ifelse((spec_levels$any_aspec==F & spec_levels$spec_levels>1), "Conflicting_double", "not_conflicting")
spec_levels$conflicting<-ifelse((spec_levels$any_aspec==T & spec_levels$spec_levels>2), "Conflicting_triple", spec_levels$conflicting)
spec_levels$conflicting<-as.factor(spec_levels$conflicting)

summary(spec_levels)

# Save nodes with conflicting splicing
write(as.character(spec_levels[which(spec_levels$conflicting=="Conflicting_double"),"geneID"]), file = file.path(newdir, "confl_splicing_double.txt"))
write(as.character(spec_levels[which(spec_levels$conflicting=="Conflicting_triple"),"geneID"]), file = file.path(newdir, "confl_splicing_triple.txt"))
write(as.character(spec_levels[which(spec_levels$conflicting%in%c("Conflicting_double", "Conflicting_triple")),"geneID"]), file = file.path(newdir, "confl_splicing.txt"))

# How many nodes are sex-specific?
prop.table(table(sexspec$spec!="Aspecific"))
# How many genes are sex-specific?
unique(sexspec[which(sexspec$spec!="Aspecific"),c("geneID","spec")]) # problem, what about conflict genes? they will be present in both.
summary(unique(sexspec[which(sexspec$spec!="Aspecific"),c("geneID","spec")])) # conflict genes listed separately for each sex, used in individual gene stats
length(unique(sexspec[which(sexspec$spec!="Aspecific"),c("geneID")])) #  list each sex-specific gene once, used for total stats

# # How many transcription and splicing nodes are sex-specific?
# table(grepl("con", sexspec_rel$eigenexonID),sexspec_rel$spec)
# prop.table(table(grepl("con", sexspec_rel$eigenexonID),sexspec_rel$spec), margin = 2)
# 
# 
# # Plot proportions
# with(droplevels(sexspec_rel[which(sexspec_rel$spec!="Aspecific"),]), table(spec, grepl("fac",eigenexonID)))
# wilcox.test(with(droplevels(sexspec_rel[which(sexspec_rel$spec!="Aspecific"),]), table(spec, grepl("fac",eigenexonID))))
# mosaicplot(with(droplevels(sexspec_rel[which(sexspec_rel$spec!="Aspecific"),]), table(spec, grepl("fac",eigenexonID))), shade = T)

##### Outdated code snippets


## Define sex specific nodes
# # No exceptions
# sexspec_abs<-data.frame(
#   eigenexonID=eigenexon_evalues[,1],
#   geneID=str_extract(eigenexon_evalues[,1], pattern = "^[^_]*"),
#   expr_F=apply(eigenexon_evalues[,which(grepl("_female", colnames(eigenexon_evalues)))],1,function(x)any(x>0)),
#   expr_M=apply(eigenexon_evalues[,which(grepl("_male", colnames(eigenexon_evalues)))],1,function(x)any(x>0))
# )
# sexspec_abs$spec<-sexspec_abs$expr_F-sexspec_abs$expr_M
# sexspec_abs$spec<-revalue(x = as.factor(sexspec_abs$spec), replace = c("-1"="Male", "1"="Female", "0"="Aspecific"))
# # How many transcription nodes are sex-specific?
# table(sexspec_abs[grep("con", sexspec_abs$eigenexonID),]$spec)
# # How many splicing nodes are sex-specific?
# table(sexspec_abs[grep("fac", sexspec_abs$eigenexonID),]$spec)
# 
# with(droplevels(sexspec_abs[which(sexspec_abs$spec!="Aspecific"),]), table(spec, grepl("fac",eigenexonID)))
# wilcox.test(with(droplevels(sexspec_abs[which(sexspec_abs$spec!="Aspecific"),]), table(spec, grepl("fac",eigenexonID))))
# mosaicplot(with(droplevels(sexspec_abs[which(sexspec_abs$spec!="Aspecific"),]), table(spec, grepl("fac",eigenexonID))), shade = T)

# # Considered as expressed only if >0 in at least 2 samples
# sexspec_rel<-data.frame(
#   eigenexonID=eigenexon_evalues[,1],
#   geneID=str_extract(eigenexon_evalues[,1], pattern = "^[^_]*"),
#   expr_F=apply(eigenexon_evalues[,which(grepl("_female", colnames(eigenexon_evalues)))],1,function(x)sum(x>0)>1),
#   expr_M=apply(eigenexon_evalues[,which(grepl("_male", colnames(eigenexon_evalues)))],1,function(x)sum(x>0)>1)
# )
# sexspec_rel$spec<-sexspec_rel$expr_F-sexspec_rel$expr_M
# sexspec_rel$spec<-revalue(x = as.factor(sexspec_rel$spec), replace = c("1"="Female", "-1"="Male", "0"="Aspecific"))

# # Which number of genes is respectively male and female specific?
# sexspec_genes<-ddply(sexspec_rel2, .(geneID), summarize, 
#                      spec = as.factor(if (all(spec!="Female")&any(spec=="Male")) "Male" else if (all(spec!="Male")&any(spec=="Female")) "Female" else "Aspecific"), 
#                      .progress="text")
# summary(sexspec_genes)

# Not working, no idea why
# spec_levels$conflicting<-as.factor(
#   if (spec_levels$any_aspec==T & spec_levels$spec_levels>2) {"Conflicting_triple" }
#   else if (spec_levels$any_aspec==F & spec_levels$spec_levels>1) {"Conflicting_double" }
#   else {"not_conflicting"}
# )
