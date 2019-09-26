#### script for detection of main gene vs alternative isoforms
## main concept
# import exon dataset, then ddply by geneID and perform hierarchical clustering
# report mean value of each cluster per sample, rank samples by overall gene expression value
# report also which exons are comprised in each transcript
### CHANGELOG: 2015/07/05, changed noise threshold from 1 to 0, only negative values are now removed

date()
library(plyr)
# import exon dataset
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/nasonia_devtesova.exon.r99.score")
row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON
# filter only rows with samples (exclude means)
rm(list=setdiff(ls(),"nasonia_devtesova.exon.r99"))
# nasdevexon<-nasonia_devtesova.exon.r99[14420:15420,grep("[.][123]",names(nasonia_devtesova.exon.r99))] # comment to run whole dataset
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
    eigenexons[[gID]]$clusters<-0
    eigenexons[[gID]]$clusterranks<-1
    eigenexons[[gID]]<-eigenexons[[gID]][c("clusters", "exonID", "clusterranks", "splicing_cathegory")]
  } else {
    Evalues<-as.data.frame(apply(Evalues, c(1,2), function(x){ifelse(x<=0, 0, x)})) # replace non-expression threshold with zeroes
    expressed_samples<-Evalues[,which(apply(Evalues, 2, function(x){any(x>0)}))]
    if (is.data.frame(expressed_samples)==F|all(is.na(expressed_samples))) {
      eigenexons[[gID]]<-as.data.frame(matrix(row.names(Evalues), ncol = 1))
      names(eigenexons[[gID]])<-"exonID"
      eigenexons[[gID]]$splicing_cathegory<-"no_expressed_samples"
      eigenexons[[gID]]$clusters<-0
      eigenexons[[gID]]$clusterranks<-1
      eigenexons[[gID]]<-eigenexons[[gID]][c("clusters", "exonID", "clusterranks", "splicing_cathegory")]
    } else {
      tree<-hcluster(Evalues, method = "correlation", link = "complete")
      #   aspecific_exons<-row.names(Evalues[apply(expressed_samples, 1, function(x){all(x>0)}),])
      # if (length(aspecific_exons)<1) {
      #   eigenexons[[gID]]<-as.data.frame(matrix(row.names(Evalues), ncol = 1))
      #   names(eigenexons[[gID]])<-"exonID"
      #   eigenexons[[gID]]$splicing_cathegory<-"no_constitutive_exons"
      #   eigenexons[[gID]]$clusters<-0
      #   eigenexons[[gID]]$clusterranks<-0
      #   eigenexons[[gID]]<-eigenexons[[gID]][c("clusters", "exonID", "clusterranks", "splicing_cathegory")]
      # } else {
      for (trans in nrow(Evalues):1){
        # evaluate relative times it is ranked as first
        Evalues$clusters<-cutree(tree, k = trans)
        clusterranks<-ddply(Evalues, .(clusters), colwise(median))
        if (nrow(clusterranks)==1){
          # there are no subranking isoforms
          eigenexons[[gID]]<-as.data.frame(matrix(row.names(Evalues), ncol = 1))
          names(eigenexons[[gID]])<-"exonID"
          eigenexons[[gID]]$clusters<-0
          eigenexons[[gID]]$clusterranks<-1
          eigenexons[[gID]]$splicing_cathegory<-"unspliced"
          eigenexons[[gID]]<-eigenexons[[gID]][c("clusters", "exonID", "clusterranks", "splicing_cathegory")]
          break} else {
            clusterranks<-apply(clusterranks[,-1], 2, function(x){rank(-round(x, digits = 3), ties.method = "min")})
            row.names(clusterranks)<-unique(Evalues$cluster)
            clusterranks<-apply(clusterranks, 1, function(x){sum(x==1)/ncol(clusterranks)})
            #               if (max(clusterranks)>=1) { # other sensitive threshold: determines the proportion of times constitutives are ranked 1st
            if (sum((clusterranks)>=0.96)==1) { # other sensitive threshold: external number determines the proportion of times constitutives are ranked 1st, external number controls there are no ties, the second clause imposes that no other clusters are ranked first in more than 80% of samples
              matchmaker<-Evalues[,"clusters", drop=F] # create table with exon subcluster assignments
              clusters<-as.data.frame(clusterranks) # store subcluster ranks
              clusters$clusters<-c(1:nrow(clusters)) #annotate subcluster ranks with subcluster IDs
              matchmaker<-join(matchmaker, clusters, by="clusters", type="left") # merge exon with subcluster ID and ranks
              row.names(matchmaker)<-row.names(Evalues) # add exon names
              matchmaker$exonID<-row.names(matchmaker)
              eigenexons[[gID]]<-matchmaker
              eigenexons[[gID]]$splicing_cathegory<-"spliced"
              eigenexons[[gID]]<-eigenexons[[gID]][c("clusters", "exonID", "clusterranks", "splicing_cathegory")]
              break} else {next}}}
    }}
}

# # this version causes problem in R below version 2
eigenexons<-ldply(eigenexons,rbind) 
names(eigenexons)[1]<-"geneID"

# ## begin safe version (messes up with the namings)
# eigenexons<-do.call(rbind.data.frame, eigenexons)
# eigenexons$geneID<-str_extract(row.names(eigenexons), "^[[:alnum:]]*")
## end safe version

# assign constitutive/specific status
eigenexons$constitutive<-ifelse((eigenexons$clusterranks>=0.96),"constitutive","facultative")
# and validate assignments
table(eigenexons$clusterranks>=0.96, eigenexons$constitutive)

# merge cluster-subcluster assignments into unique IDs with designation of constitutiveness 
eigenexons$geneID<-str_extract(eigenexons$geneID, "^[^t]*") # removes final t from gene names
eigenexons$eigenexonID<-paste0(eigenexons$geneID, "_e", eigenexons$clusters, ifelse(eigenexons$constitutive=="constitutive", "_con","_fac"), sep="")

# Save eigenexon assignments
newdir<-file.path(getwd(), "Output/exon_to_transcript_clustering_rankingonly")
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
eigenexonE<-eigenexonE[grep("single_exon_gene|spliced|unspliced", eigenexonE$splicing_cathegory),]
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

# summary(duplicated(eigenexonE_fac2$eigenexonID))
# # which genes have more expression from alternatives than from constitutives?
# which(apply(eigenexonE_fac2[,-1], 1, function(x){any(x>1)}))
# anomalies<-eigenexonE_fac2[,1][which(apply(eigenexonE_fac2[,-1], 1, function(x){any(x>1)}))]
# anomalies2<-eigenexons[which(eigenexons$eigenexonID%in%anomalies),]
# # by how much are they higher than their constitutives?
# eigenexonE_fac2[which(apply(eigenexonE_fac2[,-1], 1, function(x){any(x>1)})),]
# # which isoforms go to infinity?
# infinity<-eigenexonE_fac2[which(apply(eigenexonE_fac2[,-1], 1, function(x){any(x==Inf)})),"eigenexonID"]
# eigenexonE_cons[which(str_extract(string = infinity, pattern = "[[:alnum:]]*")%in%str_extract(eigenexonE_cons$eigenexonID, pattern = "[[:alnum:]]*")),]
# # how many genes have a sum of exon proportions greater than 1?
# eigenexonE_fac3<-eigenexonE_fac2
# eigenexonE_fac3$geneID<-str_extract(string = eigenexonE_fac3$eigenexonID, pattern = "[[:alnum:]]*_")
# a<-ddply(eigenexonE_fac3, .(geneID), numcolwise(sum))

# set infinity scores to NA, set scores greater than 1 to 1
# reasons: higher than 1 result from noise on small-scale measurements, therefore possibly complete signal plus noise
# infinity results in constitutive exons ranking below cutoff in samples where nc exons are above cutoff, therefore possibly spurious signal

eigenexonE_fac2[,sapply(eigenexonE_fac2, is.numeric)]<-apply(eigenexonE_fac2[,sapply(eigenexonE_fac2, is.numeric)],c(1,2),function(x){as.numeric(ifelse(x==Inf, NA, x))})
eigenexonE_fac2[,sapply(eigenexonE_fac2, is.numeric)]<-apply(eigenexonE_fac2[,sapply(eigenexonE_fac2, is.numeric)],c(1,2),function(x){as.numeric(ifelse(x>1, 1, x))})

## stack data-frames
eigenexonE<-rbind(eigenexonE_cons, eigenexonE_fac2)

summary(duplicated(eigenexonE$eigenexonID))

## and save as data.frame
write.csv(eigenexonE, file=file.path(newdir,"eigenexon_evalues.csv"), row.names=F)
save.image(file = file.path(newdir,"workspace.R"))

# How many isoforms per gene?
# maxclusters<-daply(eigenexons_assignments, .variables = .(geneID), summarize, max(clusters))
maxclusters<-unlist(maxclusters)
splicat<-ddply(eigenexons_assignments, .variables = .(geneID), summarize, splicing_cathegory[1], .progress = "text")
maxclusters2<-cbind(splicat, maxclusters)
names(maxclusters2)<-c("geneID", "splicing_cathegory", "maxclusters")
