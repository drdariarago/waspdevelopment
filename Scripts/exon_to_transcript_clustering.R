#### script for detection of main gene vs alternative isoforms
## main concept
# import exon dataset, then ddply by geneID and perform hierarchical clustering
# report mean value of each cluster per sample, rank samples by overall gene expression value
# report also which exons are comprised in each transcript
date()
library(plyr)
# import exon dataset
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON
# filter only rows with samples (exclude means)
rm(list=setdiff(ls(),"nasonia_devtesova.exon.r99"))
# nasdevexon<-nasonia_devtesova.exon.r99[14420:19420,grep("[.][123]",names(nasonia_devtesova.exon.r99))]
nasdevexon<-nasonia_devtesova.exon.r99[,grep("[.][123]",names(nasonia_devtesova.exon.r99))]
# Exclude gonads 
nasdevexon<-nasdevexon[,-grep("(ovaries|testes)", colnames(nasdevexon))]
# add gene IDs
library(stringr)
# nasdevexon$geneID<-unlist(regmatches(row.names(nasdevexon), gregexpr("^[^t]*",row.names(nasdevexon))))
# nasdevexon$geneID<-str_extract(row.names(nasdevexon), "^[^t]*") # base version of gene name
nasdevexon$geneID<-str_extract(row.names(nasdevexon), "^[^t]*t") # including the t postscript to disambiguate *b and non b gene names
# add exon ID
# nasdevexon$exonID<-unlist(regmatches(row.names(nasdevexon), gregexpr("[^t]*$",row.names(nasdevexon))))
nasdevexon$exonID<-str_extract(row.names(nasdevexon), "[^t]*$")

##### New workflow idea

## new idea: cluster genes according to binary distances, cut top down when all constitutives are in the same cluster. Then apply kendall clustering within each cluster bottom up and evaluate at each step if there is only one of the subclusters that is always greater than the others. 
# use hcluster from amap package
library(amap)
eigenexons<-list()
for (gID in unique(nasdevexon$geneID)){
  Evalues<-nasdevexon[grep(gID,nasdevexon$geneID),1:30]
  if (length(grep(gID,nasdevexon$geneID))<2) {
    eigenexons[[gID]]<-as.data.frame(matrix(row.names(Evalues), ncol = 1))
    names(eigenexons[[gID]])<-"exonID"
    eigenexons[[gID]]$splicing_cathegory<-"single_exon_gene"
    eigenexons[[gID]]$aspecific<-"sample_aspecific"
    eigenexons[[gID]]$subcluster<-0
    eigenexons[[gID]]$subclusterranks<-1
    eigenexons[[gID]]<-eigenexons[[gID]][c("clusters", "exonID", "aspecific", "subcluster", "subclusterranks", "splicing_cathegory")]
  } else {
    Evalues<-as.data.frame(apply(Evalues, c(1,2), function(x){ifelse(x<=0, 0, x)})) # replace non-expression with zeroes
    expressed_samples<-Evalues[,which(apply(Evalues, 2, function(x){any(x>0)}))]
    if (is.data.frame(expressed_samples)==F|all(is.na(expressed_samples))) {
      eigenexons[[gID]]<-as.data.frame(matrix(row.names(Evalues), ncol = 1))
      names(eigenexons[[gID]])<-"exonID"
      eigenexons[[gID]]$splicing_cathegory<-"no_expressed_samples"
      eigenexons[[gID]]$aspecific<-NA
      eigenexons[[gID]]$subcluster<-0
      eigenexons[[gID]]$subclusterranks<-0
      eigenexons[[gID]]<-eigenexons[[gID]][c("clusters", "exonID", "aspecific", "subcluster", "subclusterranks", "splicing_cathegory")]
      } else {
      tree<-hcluster(Evalues, method = "binary", link = "complete")
      aspecific_exons<-row.names(Evalues[apply(expressed_samples, 1, function(x){all(x>0)}),])
      if (length(aspecific_exons)<1) {
        eigenexons[[gID]]<-as.data.frame(matrix(row.names(Evalues), ncol = 1))
        names(eigenexons[[gID]])<-"exonID"
        eigenexons[[gID]]$splicing_cathegory<-"no_constitutive_exons"
        eigenexons[[gID]]$aspecific<-NA
        eigenexons[[gID]]$subcluster<-0
        eigenexons[[gID]]$subclusterranks<-0
        eigenexons[[gID]]<-eigenexons[[gID]][c("clusters", "exonID", "aspecific", "subcluster", "subclusterranks", "splicing_cathegory")]
      } else {
        Evalues$clusters<-cutree(tree, h = 0.0) # sensitive threshold of acceptable non-co-expression, expressed as proportion of samples
        eigenexons[[gID]]<-Evalues[,"clusters", drop=F]
        eigenexons[[gID]]$exonID<-row.names(eigenexons[[gID]])
        eigenexons[[gID]]$aspecific<-NA
        eigenexons[[gID]]$subcluster<-0
        eigenexons[[gID]]$subclusterranks<-1
        eigenexons[[gID]]$aspecific<-ifelse(row.names(eigenexons[[gID]])%in%aspecific_exons, "sample_aspecific", "sample_specific")
        eigenexons[[gID]]$splicing_cathegory<-"spliced"
        for (cl in unique(Evalues$clusters)){
          clEvalues<-Evalues[which(Evalues$clusters==cl),setdiff(names(Evalues), "clusters")]
          if (nrow(clEvalues)<2) {
            # is already an independently spliced exon
            next
          } else {
            t<-hcluster(clEvalues, method = "kendall", link = "complete")
            for (trans in nrow(clEvalues):1){
              # evaluate relative times it is ranked as first
              eigenexons[[gID]]$subcluster<-cutree(t, k = trans)
              eigenexons[[gID]]$subclusterranks<-ddply(clEvalues, .(subcluster), colwise(median))
              if (nrow(subclusterranks)==1){
                # there are no subranking isoforms
                break}
#               if (max(subclusterranks)>=1) { # other sensitive threshold: determines the proportion of times constitutives are ranked 1st
              if (sum((subclusterranks)>=1)==1) { # other sensitive threshold: external number determines the proportion of times constitutives are ranked 1st, external number controls there are no ties
                matchmaker<-clEvalues[,"subcluster", drop=F] # create table with exon subcluster assignments
                subcluster<-as.data.frame(subclusterranks) # store subcluster ranks
                subcluster$subcluster<-c(1:nrow(subcluster)) #annotate subcluster ranks with subcluster IDs
                matchmaker<-join(matchmaker, subcluster, by="subcluster", type="left") # merge exon with subcluster ID and ranks
                row.names(matchmaker)<-row.names(clEvalues) # add exon names
                eigenexons[[gID]][na.exclude(which(row.names(eigenexons[[gID]])%in%row.names(matchmaker))),"subcluster"]<-matchmaker$subcluster # match row.names to annotation dataset
                eigenexons[[gID]][na.exclude(which(row.names(eigenexons[[gID]])%in%row.names(matchmaker))),"subclusterranks"]<-matchmaker$subclusterranks
                break} else {next}}
          }}}}
  }
# eigenexons[[gID]]<-eigenexons[[gID]][,grep("V1", colnames(eigenexons[[gID]]), invert = T)]
}

# # this version causes problem in R below version 2
eigenexons<-ldply(eigenexons,rbind) 
names(eigenexons)[1]<-"geneID"

# ## begin safe version (messes up with the namings)
# eigenexons<-do.call(rbind.data.frame, eigenexons)
# eigenexons$geneID<-str_extract(row.names(eigenexons), "^[[:alnum:]]*")
## end safe version

# eigenexons$subcluster<-ifelse(is.na(eigenexons$subcluster), 0, eigenexons$subcluster)
# eigenexons$subclusterranks<-ifelse(is.na(eigenexons$subclusterranks), 1, eigenexons$subclusterranks)
eigenexons$constitutive<-ifelse((eigenexons["aspecific"]=="sample_aspecific"&eigenexons["subclusterranks"]=="1"),"constitutive","specific")
# and validate assignments
table(eigenexons$aspecific, eigenexons$subclusterranks==1, eigenexons$constitutive)

# merge cluster-subcluster assignments into unique IDs with designation of constitutiveness 
eigenexons$geneID<-str_extract(eigenexons$geneID, "^[^t]*") # removes final t from gene names
eigenexons$eigenexonID<-paste0(eigenexons$geneID, "_e", eigenexons$clusters, ".", eigenexons$subcluster, ifelse(eigenexons$constitutive=="constitutive", "_con","_fac"), sep="")


# identify genes with no splicing (exon to )
nospliced<-daply(eigenexons, .(geneID), function(x) {length(unique(x$eigenexonID))==1})
nospliced<-names(nospliced[which(nospliced==T)])
eigenexons[which(eigenexons$geneID%in%nospliced&eigenexons$splicing_cathegory=="spliced"),"splicing_cathegory"]<-"not_spliced"

# Save eigenexon assignments
newdir<-file.path(getwd(), "Output/exon_to_transcript_clustering")
dir.create(newdir)
write.csv(eigenexons, file=file.path(newdir,"eigenexons_assignments.csv"), row.names=F)

## average expression values based on unique eigenexon IDs
# create gene expression set
eigenexonE<-nasdevexon[,grep("male", colnames(nasdevexon))]
# eigenexonE<-eigenexonE[which(row.names(eigenexonE)%in%eigenexons$exonID),]
# replace negatives with zeroes
eigenexonE<-as.data.frame(apply(eigenexonE, c(1,2), function(x){ifelse(x<0,0,x)}))
# add exon name and merge with eigenexon assignments
eigenexonE$exonID<-row.names(eigenexonE)
eigenexonE<-merge(eigenexonE, eigenexons, by = "exonID", all.x = T)
# select only columns of interest
eigenexonE<-eigenexonE[,grep("eigenexonID|male|constitutive|splicing", colnames(eigenexonE))]
## select only valid splicing cathegories
eigenexonE<-eigenexonE[grep("single_exon_gene|spliced|not_spliced", eigenexonE$splicing_cathegory),]
## split into constitutives and non
eigenexonE_cons<-eigenexonE[which(eigenexonE$constitutive=="constitutive"),]
eigenexonE_fac<-eigenexonE[which(eigenexonE$constitutive!="constitutive"),]

summary(duplicated(eigenexonE$eigenexonID))
summary(duplicated(eigenexonE_cons$eigenexonID))
summary(duplicated(eigenexonE_fac$eigenexonID))

## average eigenexon values
eigenexonE_cons<-ddply(eigenexonE_cons, .variables = .(eigenexonID), numcolwise(median))
eigenexonE_fac<-ddply(eigenexonE_fac, .variables = .(eigenexonID), numcolwise(median))

summary(duplicated(eigenexonE_cons$eigenexonID))
summary(duplicated(eigenexonE_fac$eigenexonID))

## divide facultatives by respective constitutive value
# duplicates entries for some reasons, adding extra character at the end
eigenexonE_fac2<-apply(eigenexonE_fac, 1, function(x){
  as.numeric(x[-1])/eigenexonE_cons[grep(str_extract(x["eigenexonID"], "^[[:alnum:]]*."), eigenexonE_cons$eigenexonID),-1]
})
eigenexonE_fac2<-ldply(eigenexonE_fac2)
eigenexonE_fac2<-cbind(eigenexonE_fac[,1], eigenexonE_fac2)
names(eigenexonE_fac2)[1]<-"eigenexonID"
  
summary(duplicated(eigenexonE_fac2$eigenexonID))
# which genes have more expression from alternatives than from constitutives?
which(apply(eigenexonE_fac2[,-1], 1, function(x){any(x>1)}))
anomalies<-eigenexonE_fac2[,1][which(apply(eigenexonE_fac2[,-1], 1, function(x){any(x>1)}))]
anomalies2<-eigenexons[which(eigenexons$eigenexonID%in%anomalies),]
# almost all anomalies display a sample-specific expression pattern.
# looks like the co-expression filter enables them to escape the ranking constraints


## stack data-frames
eigenexonE<-rbind(eigenexonE_cons, eigenexonE_fac2)

summary(duplicated(eigenexonE$eigenexonID))

## and save as data.frame
write.csv(eigenexonE, file=file.path(newdir,"eigenexon_evalues.csv"), row.names=F)
save.image(file = "./Output/exon_to_transcript_clustering/workspace.R")

## diagnostics

# calculate reduction in dimensionality
with(ddply(eigenexons,.(geneID), summarize, a=length(.id), b=length(unique(eigenexonID))), sum(b)/sum(a))

# how many genes exceed distributional expectations?
eigenexonE2<-eigenexonE[grep("fac", eigenexonE$eigenexonID),]
eigenexonE2$gene<-str_extract(eigenexonE2$eigenexonID, "^[[:alnum:]]*")
# how many proportions are greater than 1?
# ggplot(melt(eigenexonE2), aes(x=value))+geom_density()+theme_bw()+scale_x_log10()+geom_rug(aes(col=value>1))
prop.table(table(melt(eigenexonE2)$value>1))
# how many sum to more than 1?
# ggplot(melt(ddply(eigenexonE2, .(gene), numcolwise(sum))), aes(x=value))+geom_density()+theme_bw()+scale_x_log10()+geom_rug(aes(col=value>1))
library(reshape2)
prop.table(table(melt(ddply(eigenexonE2, .(gene), numcolwise(sum)))$value>1))

# is sum<1 correlated to the number of exons (if yes, could be cumulative effect of error)
# ggplot(cbind(melt(ddply(eigenexonE2, .(gene), numcolwise(sum))),ddply(eigenexonE2, .(gene),nrow)), aes(x=value>1, y=V1))+geom_boxplot(notch = T, varwidth = T)+scale_y_log10()+xlab(label = "sum greater than zero?")+ylab(label = "number of eigenexons in gene")+theme_bw()
# ggplot(cbind(melt(ddply(eigenexonE2, .(gene), numcolwise(sum))),ddply(eigenexonE2, .(gene),nrow)), aes(y=value, x=V1))+geom_boxplot(aes(group=V1), notch = T, varwidth = T)+scale_y_log10()+ylab(label = "sum of proportions")+xlab(label = "number of eigenexons in gene")+theme_bw()+scale_x_log10()+geom_density2d()+geom_hline(y=1)

with(na.exclude(cbind(melt(ddply(eigenexonE2, .(gene), numcolwise(sum))),ddply(eigenexonE2, .(gene),nrow))), cor(V1, value, method = "kendall"))

# proportion of isoforms with expression greater than their constitutives
prop.table(table(unlist(eigenexonE[,-1])>1))
# proportion of isoforms expressed in absence of their constitutives
prop.table(table(is.na(unlist(eigenexonE[,-1]))))

# Which isoforms are expressed without their constitutives?
eigenexonE[which(apply(eigenexonE, 1, function(x){any(is.na(x))})),]

## there is still a significant amount of correlation across isoforms, which paired with a greater than 1 sum of proportions indicates oversplitting.
##  this will increase dimensionality, but not necessarily bias the network construction


# which gene cathegories are still present?
table(eigenexons[which(eigenexons$.id%in%str_extract(eigenexonE$eigenexonID, "^[[:alnum:]]*")),"splicing_cathegory"])

# what are the non constitutives? do they have low expression?
library(reshape2)
library(ggplot2)
aa<-melt(nasdevexon)
aa<-ddply(aa, .(geneID), function(x){median(x$value)})
bb<-unique(eigenexons[,c("geneID","splicing_cathegory")])
names(bb)<-c("geneID","splicing_cathegory")
cc<-merge(aa, bb)
ggplot(cc, aes(x=V1+2, col=splicing_cathegory))+geom_density()+theme_bw()+scale_color_brewer(type="qual", palette=3)+scale_x_log10()
ggplot(cc, aes(x=V1, fill=splicing_cathegory))+geom_histogram(position="dodge")+theme_bw()+scale_fill_brewer(type = "qual", palette = 3)
ggplot(cc, aes(y=V1, x=splicing_cathegory))+geom_boxplot(notch = T, varwidth = T)+theme_bw()
# indeed: the median is below zero, only 25% is above expression on average. Pattern consistent with almost random expression
# spliced are more expressed than non-spliced, consistent with the dampening effect of alternative exons (when not present bring down overall median)
# single exon genes span the 0 line and have greatest variation, consistent with regulation exclusively trough differential expression

# plot number of exons vs number of clusters
library(ggplot2)
library(stringr)
ggplot(ddply(eigenexons,.(geneID), summarize, a=length(geneID), b=length(unique(eigenexonID))), aes(x=a, y=b))+geom_point(position="jitter")+scale_x_log10()+scale_y_log10()+theme_bw()+geom_smooth()+geom_density2d()


# plot resuidual correlation of faculatative eigenexon in row j
j<-331
aa<-cbind(t(eigenexonE_fac[j,-1]),t(eigenexonE_cons[grep(str_extract(eigenexonE_fac[j,1], "^[[:alnum:]]*"), eigenexonE_cons$eigenexonID),-1]),t(as.numeric(eigenexonE_fac[j,-1])/eigenexonE_cons[grep(str_extract(eigenexonE_fac[j,1], "^[[:alnum:]]*"), eigenexonE_cons$eigenexonID),-1]))
aaa<-cbind(aa, aa[,1]/aa[,2])
aaa[,3]-aaa[,4]
plot(aaa[,2],aaa[,1])

### outdated code snippets

# # different versions for the reduction of transcripts to proportion data
# eigenexonE_fac2<-apply(eigenexonE_fac, 1, function(x){
#   # grep everything before the underscore on the eigenexonID of cons, then divide for matching constitutive exon group
#   #   print(paste0("looking for gene", str_extract(x[1], "^[[:alnum:]]*")))
#   cbind(as.character(x["eigenexonID"]),x[sapply(x, is.numeric)]/eigenexonE_cons[grep(str_extract(x["eigenexonID"], "^[[:alnum:]]*."), eigenexonE_cons$eigenexonID),-1])
# })
# eigenexonE_fac2<-ldply(eigenexonE_fac2)
# names(eigenexonE_fac2)<-names(eigenexonE_fac)


# # begin version based on row.names
# row.names(eigenexonE_fac)<-eigenexonE_fac$eigenexonID
# eigenexonE_fac<-eigenexonE_fac[,-grep("eigenexonID",names(eigenexonE_fac))]
# row.names(eigenexonE_cons)<-eigenexonE_cons$eigenexonID
# eigenexonE_cons<-eigenexonE_cons[,-grep("eigenexonID",names(eigenexonE_cons))]
# 
# eigenexonE_fac2<-apply(eigenexonE_fac, 1, function(x){
#   x/eigenexonE_cons[grep(str_extract(row.names(eigenexonE_fac)[x], "^[[:alnum:]]*."), row.names(eigenexonE_cons))]
# })


# # merge eigenexon assignments with exon information and rename with alphabetcal labels
# nasdevexon<-cbind(ID=row.names(nasdevexon), nasdevexon, ldply(eigenexons)[,-1])
# colnames(nasdevexon)[ncol(nasdevexon)]<-"eigenexonID"
# eigenexonlabels<-sort(apply(expand.grid(c(letters),c(letters)),1,paste0, collapse=""))[1:length(unique(nasdevexon$eigenexonID))]
# nasdevexon$eigenexonID<-as.factor(mapvalues(nasdevexon$eigenexonID, from=unique(nasdevexon$eigenexonID), to=eigenexonlabels))
# # annotate exons present in each eigenexon
# eigenID<-ddply(.data = nasdevexon, .(geneID, eigenexonID), .fun = function(x){paste(x$exonID, collapse = "&")})
# names(eigenID)<-c("geneID","eigenexonID","exonID")
# 
# 
# ## Calculate median eigenexon expression per sample, then rank to find constitutive exons
# nasdeveigenexon<-ddply(.data = nasdevexon, .variables = .(geneID, eigenexonID), .fun = colwise(median, .cols = is.numeric))
# 
# ## sum ranks of eigenexon expression, then rank according to score
# colrank<-colwise(rank, is.numeric)
# eigenranks<-cbind(nasdeveigenexon$eigenexonID, ddply(.data = nasdeveigenexon, .variables = .(geneID), .fun = colrank))
# names(eigenranks)[1:2]<-c("eigenexonID", "geneID")
# eigenranks<-cbind(eigenranks[,c(1,2)],apply(eigenranks[,-c(1,2)], 1, sum))
# names(eigenranks)[3]<-"rankings"
# # remember, those with the LOWEST sum of ranks have higher expression more often, scores of 30 mean that it is first in all samples
# # now select the constitutive eigenexon
# dlply(eigenranks, .(geneID), function(x){
#   print(x[which(x$rankings==min(x$rankings))]$eigenexonID)
#   })
# 


### trial with dlply, for some reason applies clustering to whole dataset
# eigenexons<-dlply(.data = nasdevexon,.variables = .(geneID) ,.fun = function(x){
#   # Calculate dissimilarity matrix and trees
#   tree<-hclust(dist(nasdevexon[,1:30], method="euclidean"), method="ward")
#   # Cut clusters with given height (must pick a valid threshold)
#   as.data.frame(cutree(tree, h =10))
#   tree<<-tree
#   rm(tree)
# })


############ Outdated versions

## clustering based on ranking and exclusion at the same time

# library(amap)
# eigenexons<-list()
# for (gID in unique(nasdevexon$geneID)){
#   if (length(grep(gID,nasdevexon$geneID))<2) {
#     eigenexons[[gID]]<-as.data.frame(matrix(1, nrow = 1), row.names = gID)
#     eigenexons[[gID]]$splicing_cathegory<-"single_exon_gene"
#     eigenexons[[gID]]$constitutives<-NA
#   } else {
#     Evalues<-nasdevexon[grep(gID,nasdevexon$geneID),1:30]
#     Evalues<-as.data.frame(apply(Evalues, c(1,2), function(x){ifelse(x<=0, 0, x)}))
#     expressed_samples<-Evalues[,which(apply(Evalues, 2, function(x){any(x>0)}))]
#     if (is.data.frame(expressed_samples)==F|all(is.na(expressed_samples))) {
#       eigenexons[[gID]]<-as.data.frame(matrix(1, nrow = 1), row.names = gID)
#       eigenexons[[gID]]$splicing_cathegory<-"unexpressed_no_samples"
#       eigenexons[[gID]]$constitutives<-NA
#     } else {
#       tree<-hcluster(Evalues, method = "binary")
#       aspecific_exons<-row.names(Evalues[apply(expressed_samples, 1, function(x){all(x>0)}),])
#       if (length(aspecific_exons)<1) {
#         eigenexons[[gID]]<-as.data.frame(matrix(1, nrow = 1), row.names = gID)
#         eigenexons[[gID]]$splicing_cathegory<-"unexpressed_no_constitutives"
#         eigenexons[[gID]]$constitutives<-NA
#       } else if (length(aspecific_exons)<2) {
#         for (nclust in nrow(Evalues):1){
#           tester<-cutree(tree, k=nclust)
#           const_exon_cluster<-tester[aspecific_exons]
#           tester<-sum(grepl(const_exon_cluster, tester))
#           if (tester==1){next} else {
#             eigenexons[[gID]]<-as.data.frame(cutree(tree, k=nclust+1))
#             eigenexons[[gID]]$splicing_cathegory<-"exons_single_constitutive"
#             eigenexons[[gID]]$constitutives<-NA
#             eigenexons[[gID]][grep(const_exon_cluster, eigenexons[[gID]][,1]),"constitutives"]<-"constitutive"
#             break
#           }}
#       } else {
#         Evalues$clusters<-NA
#         for (nclust in nrow(Evalues):1){
#           Evalues$clusters<-cutree(tree, k = nclust)
#           tester<-ddply(Evalues[which(row.names(Evalues)%in%aspecific_exons),], .(clusters), colwise(median))[,-1]
#           if (nrow(tester)==1){
#             # for some reason assigns constitutive status to multiple clusters
#             eigenexons[[gID]]<-as.data.frame(cutree(tree, k=nclust))
#             eigenexons[[gID]]$splicing_cathegory<-"spliced_cohesive"
#             eigenexons[[gID]]$constitutives<-NA
#             eigenexons[[gID]][grep(paste0(aspecific_exons, collapse="|"), row.names(eigenexons[[gID]])),"constitutives"]<-"constitutive"
#             break}
#           tester<-apply(tester, 2, function(x){rev(rank(x))})
#           tester<-apply(tester, 1, function(x){sum(x==1)})
#           tester<-tester/30
#           # evaluation function measures the proportion of samples in which the best eigenexon is ranked first
#           if (max(tester)==1){
#             eigenexons[[gID]]<-as.data.frame(cutree(tree, k=nclust))
#             eigenexons[[gID]]$splicing_cathegory<-"spliced_ranked"
#             eigenexons[[gID]]$constitutives<-NA
#             # currently selects index among constitutives instead of cluster ID from all clusters
#             eigenexons[[gID]][which(eigenexons[[gID]][,1]%in%which(tester==1)),"constitutives"]<-"constitutive"
#             rm(tree)
#             break} else next}
#       }}}
#   colnames(eigenexons[[gID]])<-c("eigenexonID","cathegory","constitutives")
# }


##### with only 1 constitutive exon split every exon in a transcript

# eigenexons<-list()
# for (gID in unique(nasdevexon$geneID)){
#   if (length(grep(gID,nasdevexon$geneID))<2) {
#     eigenexons[[gID]]<-as.data.frame(matrix(1, nrow = 1), row.names = gID)
#     eigenexons[[gID]]$splicing_cathegory<-"single_exon_gene"
#   } else {
#     Evalues<-nasdevexon[grep(gID,nasdevexon$geneID),1:30]
#     aspecifics<-Evalues[,which(apply(Evalues, 2, function(x){any(x>0)}))]
#     if (is.data.frame(aspecifics)==F) {
#       eigenexons[[gID]]<-as.data.frame(matrix(1, nrow = 1), row.names = gID)
#       eigenexons[[gID]]$splicing_cathegory<-"unexpressed"
#     } else {
#       aspecifics<-row.names(Evalues)[apply(aspecifics, 1, function(x){all(x>0)})]
#       if(length(aspecifics)<2){
#         tree<-hclust(dist(Evalues, method="euclidean"), method="ward")
#         eigenexons[[gID]]<-as.data.frame(cutree(tree, k=nrow(Evalues)))
#         eigenexons[[gID]]$splicing_cathegory<-"exons_single_constitutive"
#         #eigenexons[gID]$constitutives<-aspecifics
#       } else {
#         tree<-hclust(dist(Evalues, method="euclidean"), method="ward")
#         Evalues$clusters<-NA
#         for (nclust in nrow(Evalues):1){
#           Evalues$clusters<-cutree(tree, k = nclust)
#           tester<-ddply(Evalues[which(row.names(Evalues)%in%aspecifics),], .(clusters), colwise(sum))[,-1]
#           if (nrow(tester)==1){
#             eigenexons[[gID]]<-as.data.frame(cutree(tree, k=nclust))
#             eigenexons[[gID]]$splicing_cathegory<-"spliced_cohesive"
#             #eigenexons[[gID]]$constitutives<-aspecifics
#             break
#           }
#           tester<-apply(tester, 2, function(x){rev(rank(x))})
#           tester<-apply(tester, 1, function(x){sum(x==1)})
#           tester<-tester/30
#           # evaluation function measures the proportion of samples in which the best eigenexon is ranked first
#           if (max(tester)==1){
#             eigenexons[[gID]]<-as.data.frame(cutree(tree, k=nclust))
#             eigenexons[[gID]]$splicing_cathegory<-"spliced_ranked"
#             #eigenexons[[gID]]$constitutives<-"boh"
#             rm(tree)
#             break} else next}
#       }}}
#   colnames(eigenexons[[gID]])<-c("eigenexonID","cathegory")
# }



## calculate clustering for each gene (manually define threshold based on known splice variants)

# eigenexons<-list()
# for (gID in unique(nasdevexon$geneID)){
#   if (length(grep(gID,nasdevexon$geneID))<2) {
#     eigenexons[[gID]]<-as.data.frame(matrix(1, nrow = 1), row.names = gID)
#     } else {
#       tree<-hclust(distance(nasdevexon[grep(gID,nasdevexon$geneID),1:30], method="mahalanobis"), method="ward")
#       eigenexons[[gID]]<-as.data.frame(cutree(tree, h=40))
#       rm(tree)}
#   colnames(eigenexons[[gID]])<-"eigenexonID"
#   rm(gID)
# }

## old evaluation mode
## if (max(apply(colrank(ddply(Evalues, .(clusters), numcolwise(sum))),1,function(x){sum(x==nrow(Evalues))}))/30==1){
# if the number of cluster is >1 check if the lowest ranking eigenexon is ranked lowest in all samples. If yes, record result, if no decrease height, if the number of clusters is 1 decrease custer height
# validation formula: 
# min(apply(colrank(ddply(Evalues, .(clusters), numcolwise(sum))),1,sum))>30
# sums expression values of each exon in each eigenexon, then ranks each eigenexon in each sample, then sums the ranks of each eigenexon in every sample, then controls if the eigenexon with the minumal sum of ranks is 30


## p-value based alternative, ends up working as k=2
# t.test(unlist(Evalues[Evalues$clusters==1,-31]),unlist(Evalues[Evalues$clusters!=1,-31]))
## new version, for every assignment pick most expressed eigenexon and t-test against everything else. If p<0.05, use as rank 1
# 
# for (gID in unique(nasdevexon$geneID)){
#   if (length(grep(gID,nasdevexon$geneID))<2) {
#     eigenexons[[gID]]<-as.data.frame(matrix(1, nrow = 1), row.names = gID)
#   } else {
#     Evalues<-nasdevexon[grep(gID,nasdevexon$geneID),1:30]
#     tree<-hclust(dist(Evalues, method="euclidean"), method="ward")
#     Evalues$clusters<-NA
#     for (nclust in 2:nrow(Evalues)){
#       Evalues$clusters<-cutree(tree, k = nclust)
#       tester<-ddply(Evalues, .(clusters), sum)
#       tester<-tester[which(tester$V1==min(tester$V1)),"clusters"]
#       # evaluate if the most expressed cluster is significantly higher than the rest
#       if (t.test(unlist(Evalues[Evalues$clusters==tester,-31]),unlist(Evalues[Evalues$clusters!=tester,-31]))$p.value<0.0001){
#         eigenexons[[gID]]<-as.data.frame(cutree(tree, k=nclust))
#         rm(tree)
#         break} else if (nclust==nrow(Evalues)) {
#           eigenexons[[gID]]<-as.data.frame(matrix(1, nrow = 1), row.names = gID)
#         } else next}
#   }
#   colnames(eigenexons[[gID]])<-"eigenexonID"
#   rm(gID)
# }
# 
# #### version based in proportion of constitutive cluster being ranked 1st
# 
# for (gID in unique(nasdevexon$geneID)){
#   if (length(grep(gID,nasdevexon$geneID))<2) {
#     eigenexons[[gID]]<-as.data.frame(matrix(1, nrow = 1), row.names = gID)
#   } else {
#     Evalues<-nasdevexon[grep(gID,nasdevexon$geneID),1:30]
#     tree<-hclust(dist(Evalues, method="euclidean"), method="ward")
#     Evalues$clusters<-NA
#     for (nclust in nrow(Evalues):1){
#       if (nclust==1){
#         eigenexons[[gID]]<-as.data.frame(matrix(1, nrow = 1), row.names = gID)
#       } else {
#         Evalues$clusters<-cutree(tree, k = nclust)
#         tester<-ddply(Evalues, .(clusters), numcolwise(sum))[,-1]
#         tester<-apply(tester, 2, function(x){rev(rank(x))})
#         tester<-apply(tester, 1, function(x){sum(x==1)})
#         tester<-tester/30
#         # evaluation function measures the proportion of samples in which the best eigenexon is ranked first
#         if (max(tester)==1){
#           eigenexons[[gID]]<-as.data.frame(cutree(tree, k=nclust))
#           rm(tree)
#           break} else next}
#     }}
#   colnames(eigenexons[[gID]])<-"eigenexonID"
#   rm(gID)
# }


# how many ties are there within each gene?
ties<-daply(nasdevexon[,], .(geneID), numcolwise(function(x){(sum(duplicated(x)))}))
sum(unlist(ties))/length(unlist(ties))
ties2<-daply(nasdevexon[,], .(geneID), numcolwise(function(x){(sum(duplicated(round(x, digits=2))))}))
sum(unlist(ties2))/length(unlist(ties2))
ties3<-daply(nasdevexon[,], .(geneID), numcolwise(function(x){(sum(duplicated(round(x, digits=3))))}))
sum(unlist(ties3))/length(unlist(ties3))
ties4<-daply(nasdevexon[,], .(geneID), numcolwise(function(x){(sum(duplicated(round(x, digits=4))))}))
sum(unlist(ties4))/length(unlist(ties4))
ties5<-daply(nasdevexon[,], .(geneID), numcolwise(function(x){(sum(duplicated(round(x, digits=5))))}))
sum(unlist(ties5))/length(unlist(ties5))

