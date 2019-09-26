## Draft for analysis of clustering of sex-de genes, based on unidimensional spatial distances

## Ideal workflow:
# Create vector of gene positions, identify number of genes between any pair of gene (distance matrix, also possible in kb or cM), index all sex DE genes and calculate distance to nearest sex DE gene. Look at distribution
# Model to be analyzed as poisson/negative binomial (discrete occurrences) or random walk

## Final method workflow: import positions as vector, calculate differences and store in matrix, calculate closeness (1/abs(d)) of next de gene (non-de excluded via multiplication by zero). Store position of closest match, absolute distance and number of intervening non DE genes (distance<next de distance)
## MUST either add correction for genes on margin of scaffold OR remove them from test set
## Create vector of gene positions
a<-sample(230000, 1000, replace=F)
# create vector of sex de
de<-rbinom(1000,1,0.5)

# Create distance matrix between any pair of genes
d<-sapply(a, function(x){x-a})
# and remove diagonal
d[which(d==0)]<-NA
# for all sex de find closest other sex de, distance and number of genes between the two
closest<-data.frame(genepos=rep(NA, length(de)), distance=rep(NA, length(de)), between=rep(NA, length(de)))
for (i in 1:length(de)){
  if (de[i]==1){
    # find closest by max 1/distance, remove non de multiplying by zero
    closest$genepos[i]<-which.max((1/abs(d[i,]))*de)
    closest$distance[i]<-d[i,closest$genepos[i]]
    # If closest is positive, find lower distances (greater than zero)
    if (closest$distance[i]>0){
      closest$between[i]<-sum(d[i,which(d[i,]>0)]<closest$distance[i])}
    # If closest is negative, find higher distances (lower than zero)
    else {
      closest$between[i]<-sum(d[i,which(d[i,]<0)]>closest$distance[i])}
    closest$distance[i]<-abs(closest$distance[i])
  }
}
closest

library(lattice)
densityplot(closest$distance)
densityplot(log10(closest$distance))
densityplot(closest$between)
densityplot(log10(closest$between))
plot(log10(closest$distance), closest$between)


### Remaining problems: different chromosomes, infinite distances. Solutions
# a: create list of positions and apply to list, then merge results
# b: store also chromosome information. if chromosomes differ, set distance to NA

# Test of patchyness: compare with random distribution of betweenness (poisson, lambda=observed proportion of DE)
# Expected betweenness: clustered<1/p(DE)<overdispersed
# Expected relationship between distance and betweennes: between genes=(1-p(DE))*genes/bp*bp (does not account for gene size)

## IMPORTANT NOTE: these additional metrics (absolute distance and between genes) can be used in the binomial model to test for medium-scale spatial correlation of gene expression (i.e. chromatine remodelling?) or bidirectional enhancers

## Multivariate alternative: compare topological overlap for distance matrix and gene expression correlation matrix