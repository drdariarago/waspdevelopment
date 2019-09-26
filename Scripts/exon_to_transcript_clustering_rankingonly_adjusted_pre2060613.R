#### script for detection of main gene vs alternative isoforms
## main concept
# import exon dataset, then ddply by geneID and perform hierarchical clustering
# report mean value of each cluster per sample, rank samples by overall gene expression value
# report also which exons are comprised in each transcript
### CHANGELOG: 2015/07/05, changed noise threshold from 1 to 0, only negative values are now removed
### CHANGELOG: 2015/07/11, remove exons that exceed noise threshold in less than 3 samples before clustering, add new class (single expressed exon)
### CHANGELOG: 2015/07/13, average duplicate exons (present in multiple transcripts) by using their medians
### CHANGELOG: 2015/07/19, redirect to correct score files, remove duplicated exons (one exon per gene rather than one exon per transcript)
### CHANGELOG: 2015/07/28, remove exons expressed in less than 3 samples before clustering loop, replaced gene matching function grep with %in% (disambiguates between * and *b genes)

## Initialize script
date()
library(plyr)
library(stringr)
library(amap)

# Import exon dataset
nasonia_devtesova.exon.r99 <- read.delim("./Nasonia Score files/2015_05_19_Justin/nasonia_devtesova.exon.r99.score2")
rm(list=setdiff(ls(),"nasonia_devtesova.exon.r99"))

## Create output directories
newdir<-file.path(getwd(), "Output/exon_to_transcript_clustering_rankingonly_adjusted")
dir.create(newdir)
# graphdir<-file.path(getwd(), "Graphics/exon_to_transcript_clustering_rankingonly_adjusted")
# dir.create(graphdir)

# filter only rows with samples (exclude means)
row.names(nasonia_devtesova.exon.r99)<-nasonia_devtesova.exon.r99$EXON
nasdevexon<-nasonia_devtesova.exon.r99[,grep("[.][123]",names(nasonia_devtesova.exon.r99))]
# nasdevexon<-nasonia_devtesova.exon.r99[14420:15420,grep("[.][123]",names(nasonia_devtesova.exon.r99))] # Uncomment to run reduced testing dataset

# Exclude gonads 
nasdevexon<-nasdevexon[,-grep("(ovaries|testes)", colnames(nasdevexon))]
# add gene IDs
nasdevexon$geneID<-str_extract(row.names(nasdevexon), "^[^t]*") # including the t postscript to disambiguate *b and non b gene names

# nasdevexon$geneID<-str_extract(row.names(nasdevexon), "^[^t]*") # base version of gene name
# nasdevexon$geneID<-unlist(regmatches(row.names(nasdevexon), gregexpr("^[^t]*",row.names(nasdevexon)))) # version without stringr

# add exon ID
nasdevexon$exonID<-str_extract(row.names(nasdevexon), "[[:digit:]]*$")
nasdevexon$ExonID<-paste(nasdevexon$geneID, nasdevexon$exonID, sep = "e")

# nasdevexon$exonID<-unlist(regmatches(row.names(nasdevexon), gregexpr("[^t]*$",row.names(nasdevexon)))) # version without stringr

# Retain only one unique exon per gene
nasdevexon<-nasdevexon[which(duplicated(nasdevexon$ExonID)==F),]

# Remove exons expressed in less than 3 samples
nrow(nasdevexon)
length(unique(nasdevexon$geneID))

nasdevexon<-nasdevexon[which(apply(nasdevexon[,1:30],1,function(x){sum(x>0)>2})),]

nrow(nasdevexon)
length(unique(nasdevexon$geneID))

##### Splicing detection

## Cluster genes according to reciprocal correlations, then iteratively cut tree (bottom up) until one cluster is the most expressed or tied for expression across all samples.
# Remove non-expressed exons before calculating correlations
# use hcluster from amap package

eigenexons<-list()
for (gID in unique(nasdevexon$geneID)){
  Evalues<-nasdevexon[which(nasdevexon$geneID%in%gID),setdiff(names(nasdevexon), c("geneID","exonID","ExonID"))]
  if (length(grep(gID,nasdevexon$geneID))<2) {
    eigenexons[[gID]]<-as.data.frame(matrix(row.names(Evalues), ncol = 1))
    names(eigenexons[[gID]])<-"exonID"
    eigenexons[[gID]]$splicing_cathegory<-"single_exon_gene"
    eigenexons[[gID]]$clusters<-0
    eigenexons[[gID]]$clusterranks<-1
    eigenexons[[gID]]<-eigenexons[[gID]][c("clusters", "exonID", "clusterranks", "splicing_cathegory")]
  } else {
    Evalues<-as.data.frame(apply(Evalues, c(1,2), function(x){ifelse(x<=0, 0, x)})) # replace non-expression threshold with zeroes
    Evalues2<-Evalues[apply(Evalues, 1, function(x){sum(x>0)>=3}),] # remove exons expressed in less than three samples
    expressed_samples<-Evalues2[,which(apply(Evalues2, 2, function(x){any(x>0)}))]
    if (is.data.frame(expressed_samples)==F|all(is.na(expressed_samples))|nrow(expressed_samples)<1) {
      eigenexons[[gID]]<-as.data.frame(matrix(row.names(Evalues), ncol = 1))
      names(eigenexons[[gID]])<-"exonID"
      eigenexons[[gID]]$splicing_cathegory<-"no_expressed_samples"
      eigenexons[[gID]]$clusters<-0
      eigenexons[[gID]]$clusterranks<-1
      eigenexons[[gID]]<-eigenexons[[gID]][c("clusters", "exonID", "clusterranks", "splicing_cathegory")]
    } else if (nrow(expressed_samples)==1) {
      eigenexons[[gID]]<-as.data.frame(matrix(row.names(Evalues2), ncol = 1))
      names(eigenexons[[gID]])<-"exonID"
      eigenexons[[gID]]$splicing_cathegory<-"single_expressed_exon"
      eigenexons[[gID]]$clusters<-0
      eigenexons[[gID]]$clusterranks<-1
      eigenexons[[gID]]<-eigenexons[[gID]][c("clusters", "exonID", "clusterranks", "splicing_cathegory")]
    } else {
      tree<-hcluster(Evalues2, method = "correlation", link = "complete")
      for (trans in nrow(Evalues2):1){
        # evaluate relative times it is ranked as first
        Evalues2$clusters<-cutree(tree, k = trans)
        clusterranks<-ddply(Evalues2, .(clusters), colwise(median))
        if (nrow(clusterranks)==1){
          # there are no subranking isoforms
          eigenexons[[gID]]<-as.data.frame(matrix(row.names(Evalues2), ncol = 1))
          names(eigenexons[[gID]])<-"exonID"
          eigenexons[[gID]]$clusters<-0
          eigenexons[[gID]]$clusterranks<-1
          eigenexons[[gID]]$splicing_cathegory<-"unspliced"
          eigenexons[[gID]]<-eigenexons[[gID]][c("clusters", "exonID", "clusterranks", "splicing_cathegory")]
          break} else {
            clusterranks<-apply(clusterranks[,-1], 2, function(x){rank(-round(x, digits = 3), ties.method = "min")})
            row.names(clusterranks)<-unique(Evalues2$cluster)
            clusterranks<-apply(clusterranks, 1, function(x){sum(x==1)/ncol(clusterranks)})
            #               if (max(clusterranks)>=1) { # other sensitive threshold: determines the proportion of times constitutives are ranked 1st
            if (sum((clusterranks)>=0.96)==1) { # other sensitive threshold: external number determines the proportion of times constitutives are ranked 1st, external number controls there are no ties, the second clause imposes that no other clusters are ranked first in more than 80% of samples
              matchmaker<-Evalues2[,"clusters", drop=F] # create table with exon subcluster assignments
              clusters<-as.data.frame(clusterranks) # store subcluster ranks
              clusters$clusters<-c(1:nrow(clusters)) #annotate subcluster ranks with subcluster IDs
              matchmaker<-join(matchmaker, clusters, by="clusters", type="left") # merge exon with subcluster ID and ranks
              row.names(matchmaker)<-row.names(Evalues2) # add exon names
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
# ## end safe version

# assign constitutive/specific status
eigenexons$constitutive<-ifelse((eigenexons$clusterranks>=0.96),"constitutive","facultative")
# and validate assignments
table(eigenexons$clusterranks>=0.96, eigenexons$constitutive)

# merge cluster assignments into unique IDs with designation of constitutiveness 
# eigenexons$geneID<-str_extract(eigenexons$geneID, "^[^t]*") # removes final t from gene names
eigenexons$eigenexonID<-paste0(eigenexons$geneID, "_e", eigenexons$clusters, ifelse(eigenexons$constitutive=="constitutive", "_con","_fac"), sep="")

# Save eigenexon assignments
write.csv(eigenexons, file=file.path(newdir,"eigenexons_assignments.csv"), row.names=F)

## Average expression values based on unique eigenexon IDs
# create gene expression set
eigenexonE<-nasdevexon[,grep("male", colnames(nasdevexon))]
# eigenexonE<-eigenexonE[which(row.names(eigenexonE)%in%eigenexons$exonID),]

# Replace negatives with zeroes
eigenexonE<-as.data.frame(apply(eigenexonE, c(1,2), function(x){ifelse(x<0,0,x)}))
# Add exon name and merge with eigenexon assignments
eigenexonE$exonID<-row.names(eigenexonE)
eigenexonE<-merge(eigenexonE, eigenexons, by = "exonID", all.x = F, all.y=T)
# select only columns of interest
eigenexonE<-eigenexonE[,grep("eigenexonID|male|constitutive|splicing", colnames(eigenexonE))]
## select only expressed nodes
eigenexonE<-eigenexonE[grep("single_exon_gene|spliced|unspliced|single_expressed_exon", eigenexonE$splicing_cathegory),]
## split into constitutives and facultatives
eigenexonE_cons<-eigenexonE[which(eigenexonE$constitutive=="constitutive"),]
eigenexonE_fac<-eigenexonE[which(eigenexonE$constitutive!="constitutive"),]

summary(duplicated(eigenexonE$eigenexonID))
summary(duplicated(eigenexonE_cons$eigenexonID))
summary(duplicated(eigenexonE_fac$eigenexonID))

## average eigenexon values
eigenexonE_cons<-ddply(eigenexonE_cons, .variables = .(eigenexonID), numcolwise(median), na.rm = T)
eigenexonE_fac<-ddply(eigenexonE_fac, .variables = .(eigenexonID), numcolwise(median), na.rm = T)

summary(duplicated(eigenexonE_cons$eigenexonID))
summary(duplicated(eigenexonE_fac$eigenexonID))

## divide facultatives by respective constitutive value
# duplicates entries since * genes match to *b genes as well, adding extra character at the end
eigenexonE_fac2<-apply(eigenexonE_fac, 1, function(x){
  as.numeric(x[-1])/eigenexonE_cons[grep(str_extract(x["eigenexonID"], "^[[:alnum:]]*."), eigenexonE_cons$eigenexonID),-1]
})
eigenexonE_fac2<-ldply(eigenexonE_fac2)
eigenexonE_fac2<-cbind(eigenexonE_fac[,1], eigenexonE_fac2)
names(eigenexonE_fac2)[1]<-"eigenexonID"

# set infinity and not a number scores to 0, set scores greater than 1 to 1
# reasons: higher than 1 result from noise on small-scale measurements, therefore possibly complete signal plus noise
# infinity results in constitutive exons ranking below cutoff in samples where nc exons are above cutoff, therefore possibly spurious signal. NaN arise from 0/0 errors.

eigenexonE_fac2[,sapply(eigenexonE_fac2, is.numeric)]<-apply(eigenexonE_fac2[,sapply(eigenexonE_fac2, is.numeric)],c(1,2),function(x){as.numeric(ifelse(is.na(x), 0, x))})
eigenexonE_fac2[,sapply(eigenexonE_fac2, is.numeric)]<-apply(eigenexonE_fac2[,sapply(eigenexonE_fac2, is.numeric)],c(1,2),function(x){as.numeric(ifelse(x>1, 1, x))})

## stack data-frames
eigenexonE<-rbind(eigenexonE_cons, eigenexonE_fac2)
summary(duplicated(eigenexonE$eigenexonID))

# Remove nodes expressed in less than 3 samples
goodRows<-apply(eigenexonE[,-1], 1, function(x){sum(na.exclude(x)>0)>2})
summary(goodRows)
eigenexonE[which(goodRows==F),]
eigenexonE<-eigenexonE[which(goodRows==T),]

## and save as data.frame
write.csv(eigenexonE, file=file.path(newdir,"eigenexon_evalues.csv"), row.names=F)
save.image(file = file.path(newdir,"workspace.R"))

###

# # Detect transcripts with NA values
# 
# eigenexonE[which(apply(eigenexonE, 1, function(x){any(is.na(x))})),]
# eigenexonE[which(apply(eigenexonE, 1, function(x){any(is.na(x))})),"eigenexonID"]
# 
# eigenexons[which(eigenexons$eigenexonID%in%eigenexonE[which(apply(eigenexonE, 1, function(x){any(is.na(x))})),"eigenexonID"]),]
# 
# # filter only genes with NAs
# eigenexons<-eigenexons[which(eigenexons$geneID%in%eigenexons[which(eigenexons$eigenexonID%in%eigenexonE[which(apply(eigenexonE, 1, function(x){any(is.na(x))})),"eigenexonID"]),"geneID"]),]
