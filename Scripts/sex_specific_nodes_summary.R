## Analyses for sex-specific nodes section
# Initialize script
date()
rm(list=ls())
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)
# initialize output path
newdir<-file.path(getwd(), "Output/sex_specific_nodes_summary")
dir.create(newdir)
graphdir<-file.path(getwd(), "Graphics/sex_specific_nodes_summary")
dir.create(graphdir)

# Load dataset
clusterdata <- read.csv(file = "./Output/Results_compiler/clusterdata_full.csv")[,-1]
# How many nodes in network are sex-spec vs non sex-spec?
summary(is.na(clusterdata$spec))
# How many are male vs female vs non-spec?
summary(clusterdata$spec)
# From how many genes?
length(unique(clusterdata[which(is.na(clusterdata$spec)==F),"geneID"]))
length(unique(clusterdata[which(clusterdata$spec=="Female"),"geneID"]))
length(unique(clusterdata[which(clusterdata$spec=="Male"),"geneID"]))
##### Compare distributions vs trans/spl (must normalize to proportion of spl to trans nodes in network)
table(clusterdata$spec, grepl("_fac",clusterdata$eigenexonID), dnn = c("spec","spl"))
prop.table(table(clusterdata$spec, grepl("_fac",clusterdata$eigenexonID), dnn = c("spec","spl")))
prop.table(table(clusterdata$spec, grepl("_fac",clusterdata$eigenexonID), dnn = c("spec","spl")), margin = 1)

# How many splicing vs transcription nodes are sexspec?
table(clusterdata$node_type, is.na(clusterdata$spec)==F, dnn = c("Node type","Sex Specific"))
prop.table(table(clusterdata$node_type, is.na(clusterdata$spec)==F, dnn = c("Node type","Sex Specific")), margin = 1)
fisher.test(table(clusterdata$node_type, is.na(clusterdata$spec)==F)) # transcription nodes have lower prob of being sex-specific than splicing nodes
# How many per each sex?
table(clusterdata$node_type, clusterdata$spec, dnn = c("Node type","Sex Specific"))
prop.table(table(clusterdata$node_type, clusterdata$spec, dnn = c("Node type","Sex Specific")), margin = 1)
prop.table(table(clusterdata$node_type, clusterdata$spec, dnn = c("Node type","Sex Specific")), margin = 2)
fisher.test(table(clusterdata$node_type, clusterdata$spec, dnn = c("Node type","Sex Specific")))

### Restrict to transcription nodes
dupl_data_con <- clusterdata[grep("_con",clusterdata$eigenexonID),c("geneID","spec","quality7","Copynumber","odb8_og_id","node_type")]
# Test whether Sexspec nodes are enriched in duplicated genes
spec_par_table_con <- table(is.na(dupl_data_con$spec)==F,dupl_data_con$quality7=="Paralog", dnn = c("spec","paralog"))
fisher.test(x = spec_par_table_con)
prop.table((spec_par_table_con), margin = 1)
margin.table(prop.table(spec_par_table_con), margin = 2)
# Test whether Sexspec nodes are enriched in duplicated genes (restrict to splicing nodes)
dupl_data_fac <- clusterdata[grep("_fac",clusterdata$eigenexonID),c("geneID","spec","quality7")]
spec_par_table_fac <- table(is.na(dupl_data_fac$spec)==F,dupl_data_fac$quality7=="Paralog", dnn = c("spec","paralog"))
fisher.test(x = spec_par_table_fac)
prop.table((spec_par_table_fac), margin = 1)
margin.table(prop.table(spec_par_table_fac), margin = 2)
# Export non-duplicated and duplicated list of sex-specific genes for males and females to test whether non-duplicated ones show more structural proteins
write.csv(x = dupl_data_con[which(is.na(dupl_data_con$spec)==F),], file = file.path(newdir, "Trans_Sexspe_for_enrichment.csv"))

## As above, with ODB8 paralog number in OG
spec_par_table_con <- table(is.na(dupl_data_con$spec)==F,dupl_data_con$Copynumber>1, dnn = c("spec","paralog"))
fisher.test(x = spec_par_table_con)
prop.table((spec_par_table_con), margin = 1)
margin.table(prop.table(spec_par_table_con), margin = 2)

# Test whether Sexspec nodes are enriched in duplicated genes (restrict to splicing nodes)
dupl_data_fac <- clusterdata[grep("_fac",clusterdata$eigenexonID),c("geneID","spec","Copynumber")]
spec_par_table_fac <- table(is.na(dupl_data_fac$spec)==F,dupl_data_fac$Copynumber>1, dnn = c("spec","paralog"))
fisher.test(x = spec_par_table_fac)
prop.table((spec_par_table_fac), margin = 1)
margin.table(prop.table(spec_par_table_fac), margin = 2)

# Check whether duplicated sexspec genes in males and females overlap
dupl_con_OGs <- ddply(.data = clusterdata[which(grepl("_con",clusterdata$eigenexonID)&clusterdata$Copynumber>1),c("odb8_og_id","spec")], .variables = .(odb8_og_id), summarize, malespec = any(spec=="Male"), femalespec = any(spec=="Female"))
dupl_con_OGs$malespec <- ifelse(is.na(dupl_con_OGs$malespec), F, dupl_con_OGs$malespec)
dupl_con_OGs$femalespec <- ifelse(is.na(dupl_con_OGs$femalespec), F, dupl_con_OGs$femalespec)
dupl_con_OGs <- dupl_con_OGs[-nrow(dupl_con_OGs),]
dupl_con_OGs<-droplevels(dupl_con_OGs)

table(dupl_con_OGs$malespec, dupl_con_OGs$femalespec, dnn = c("malespec","femalespec"))
fisher.test(table(dupl_con_OGs$malespec, dupl_con_OGs$femalespec, dnn = c("malespec","femalespec")))
prop.table(table(dupl_con_OGs$malespec, dupl_con_OGs$femalespec, dnn = c("malespec","femalespec")))
dupl_con_OGs[which(dupl_con_OGs$malespec==T&dupl_con_OGs$femalespec==T),]

# Check for transposon/retrotransposon domains, remove and re-test
transposons <- c("EOG87H867","EOG8B8MV0","EOG8CG2R8","EOG8CNT6C","EOG8D2962","EOG8DRCV3","EOG8FJBS3","EOG8HTC74","EOG8M3CQ4","EOG8MGVNF","EOG8NZX77","EOG8VX4M0","EOG8WWV1R")
dupl_con_OGs <- dupl_con_OGs[-which(dupl_con_OGs$odb8_og_id%in%transposons),]
table(dupl_con_OGs$malespec, dupl_con_OGs$femalespec, dnn = c("malespec","femalespec"))
fisher.test(table(dupl_con_OGs$malespec, dupl_con_OGs$femalespec, dnn = c("malespec","femalespec")))
spec_dupl_prop<-prop.table(table(dupl_con_OGs$malespec, dupl_con_OGs$femalespec, dnn = c("malespec","femalespec")))
# compare with expected number of spec in both
spec_dupl_prop[4]/(spec_dupl_prop[2]*spec_dupl_prop[3]) # ~4.5 fold higher than expected
# List all OGs with matching duplicates
write.csv(x = dupl_data_con[which(dupl_data_con$odb8_og_id%in%dupl_con_OGs[which(dupl_con_OGs$malespec==T&dupl_con_OGs$femalespec==T),"odb8_og_id"]),], file = file.path(newdir,"MF_duplicated_OG_con_nodes.csv"))

### What is the sexspec status of parent genes of sexspec splicing nodes?
# Define spec and bias status of parent genes
parentstatus <- clusterdata[grep("_con",clusterdata$eigenexonID), c("geneID","spec","sexbias","devsexbias")]
# add parent status to splicing nodes
splstatus <- merge(clusterdata, parentstatus, by="geneID", suffixes = c("_spl","_trans"))
# Fill NAs with unbiased
splstatus$spec_spl <- as.factor(ifelse(is.na(splstatus$spec_spl), "Aspecific", as.character(splstatus$spec_spl)))
splstatus$spec_trans <- as.factor(ifelse(is.na(splstatus$spec_trans), "Aspecific", as.character(splstatus$spec_trans)))
# Restrict to splicing nodes
splstatus <- splstatus[which(splstatus$node_type=="splicing"),]

# Compare sex-specific splicing nodes to parent's splicing and bias status
table(ifelse(splstatus$spec_trans=="Aspecific","Aspecific","Sexspec"),
      ifelse((is.na(splstatus$devsexbias_trans)|splstatus$devsexbias_trans=="....."),"Unbiased","Sexbiased"), 
      ifelse(splstatus$spec_spl=="Aspecific", "Aspecific","Sexspec"), 
      dnn = c("trans_spec","trans_bias","spl_spec"),
      useNA = "ifany")

prop.table(table(ifelse(splstatus$spec_trans=="Aspecific","Aspecific","Sexspec"), 
                 ifelse((is.na(splstatus$devsexbias_trans)|splstatus$devsexbias_trans=="....."),"Unbiased","Sexbiased"), 
                 ifelse(splstatus$spec_spl=="Aspecific", "Aspecific","Sexspec"), 
                 dnn = c("trans_spec","trans_bias","spl_spec")), 
           margin = 3)


# compare spec nodes to parent's spec
table(splstatus$spec_spl, splstatus$spec_trans, dnn = c("specific","parent"))
prop.table(table(splstatus$spec_spl, splstatus$spec_trans, dnn = c("specific","parent")), margin = 1)
prop.table(table(splstatus$spec_spl, splstatus$spec_trans, dnn = c("specific","parent")), margin = 2)
# Only 5% of fem spec belong to fem spec transc genes, only 14% of malespec belong to malespec trans genes

# Test for association between parent spec and child spec
table(splstatus$spec_spl!="Aspecific", splstatus$spec_trans!="Aspecific", dnn = c("specific","parent_specific"))
prop.table(table(splstatus$spec_spl!="Aspecific", splstatus$spec_trans!="Aspecific", dnn = c("specific","parent_specific")))
prop.table(table(splstatus$spec_spl!="Aspecific", splstatus$spec_trans!="Aspecific", dnn = c("specific","parent_specific")), margin = 1)
# Conclusion: Proportion of sexspec splicing nodes generated by sexspec transcription nodes ~11, similar to overall 10% of sexspec spl nodes in transcriptome

# Compare spec nodes to non-spec parent's bias
table(splstatus$spec_spl!="Aspecific", is.na(splstatus$sexbias_trans)==F, dnn = c("spliced_sexspec","parent_sexbias"))
fisher.test(table(splstatus$spec_spl!="Aspecific", is.na(splstatus$sexbias_trans)==F, dnn = c("spliced_sexspec","parent_sexbias")))
prop.table(table(splstatus$spec_spl!="Aspecific", is.na(splstatus$sexbias_trans)==F, dnn = c("spliced_sexspec","parent_sexbias")), margin = 2)
prop.table(table(splstatus$spec_spl!="Aspecific", is.na(splstatus$sexbias_trans)==F, dnn = c("spliced_sexspec","parent_sexbias")), margin = 1)
# Conclusion: sexbiased parents generate more frequently sexspec splicnig nodes, small effect size (~3% higher, 3 fold higher than the expected by independent probabilities)

# Compare conflict betw parent bias and spec nodes
table(splstatus$spec_spl, splstatus$sexbias_trans, dnn = c("sexspec","parent_sexbias"))
prop.table(table(splstatus$spec_spl, splstatus$sexbias_trans, dnn = c("sexspec","parent_sexbias")), margin = 1)
# Check for non developmentally conflicting genes only
splstatus_MF <- droplevels(splstatus[which(splstatus$sexbias_trans!="Female, Male"),])
table(splstatus_MF$spec_spl, splstatus_MF$sexbias_trans, dnn = c("spliced_sexspec","parent_sexbias"))
fisher.test(table(splstatus_MF$spec_spl, splstatus_MF$sexbias_trans, dnn = c("spliced_sexspec","parent_sexbias")))
prop.table(table(splstatus_MF$spec_spl, splstatus_MF$sexbias_trans, dnn = c("spliced_sexspec","parent_sexbias"))[-1,])
prop.table(table(splstatus_MF$spec_spl, splstatus_MF$sexbias_trans, dnn = c("spliced_sexspec","parent_sexbias"))[-1,], margin = 2)
# Conclusion: Strong association between direction of transcriptional bias and direction of splicing specificity, greater in males


# Check for most frequent expr patterns of parent genes with sexspec spl nodes
table(droplevels(splstatus$devsexbias_trans), droplevels(splstatus$spec_spl), dnn = c("parent_sexbias","spliced_sexspec"))
prop.table(table(splstatus$spec_spl, splstatus$devsexbias_trans, dnn = c("spliced_sexspec","parent_sexbias")), margin = 2)

## Plot sex-spec genes per stage expr
spec_bystage <- clusterdata[which(clusterdata$spec!="Aspecific"),c("geneID","eigenexonID","spec","node_type","dev_expr")]
spec_stage <- sapply(X = spec_bystage$dev_expr, FUN = str_locate_all, pattern="1")
names(spec_stage) <- spec_bystage$eigenexonID
spec_stage <- ldply(spec_stage)
spec_stage <- spec_stage[,-3]
names(spec_stage) <- c("eigenexonID","stage")
spec_stage$stage <- mapvalues(x = as.factor(spec_stage$stage), from = 1:5, to = c("emb10","emb18","lar51","pupyel","adult"))
spec_bystage <- merge(spec_bystage, spec_stage, by="eigenexonID")
spec_bystage <- ddply(.data = spec_bystage, .variables = .(spec, node_type, stage), summarize, count = length(eigenexonID))
spec_bystage$grouper <- as.factor(apply(X = spec_bystage[,1:2], MARGIN = 1,  paste, collapse="_"))

pdf(file = file.path(graphdir, "SexSpec_nodecounts_line.pdf"))
ggplot(data = spec_bystage, aes(x=stage, y=count, col=spec, lty=node_type, group=grouper))+geom_line()+geom_point()+theme_bw()+scale_y_continuous(trans = "log10", breaks=c(5,10,50,100,500,1000,1500,3000), name = "Number of Sex-Specific Nodes", limits=c(3,3500))+scale_x_discrete(name = "Stage", labels = c("Early Embryo","Late Embryo","Larva", "Pupa", "Adult"))+scale_color_discrete(name = "Sex")+scale_linetype_discrete(name = "Node Type")+ggtitle(label = "Number of Sex-Specific Nodes per Stage\n")
dev.off()
