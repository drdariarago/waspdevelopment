# The clustering coefficient measures the cliquishness of a gene. Many references use this concept.
# For our definition of the clustering coefficient in weighted networks consult Zhang and Horvath
#(2005) and Dong and Horvath (2007) 

# The function ClusterCoef.fun computes the cluster coefficients.
# Input is an adjacency matrix 
if(exists("ClusterCoef.fun")) rm(ClusterCoef.fun) ; ClusterCoef.fun=function(adjmat1) {
  diag(adjmat1)=0
  no.nodes=dim(adjmat1)[[1]]
  computeLinksInNeighbors <- function(x, imatrix){x %*% imatrix %*% x}
  nolinksNeighbors <- c(rep(-666,no.nodes))
  total.edge <- c(rep(-666,no.nodes))
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
  if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
    nolinksNeighbors <- apply(adjmat1, 1, computeLinksInNeighbors, imatrix=adjmat1)
    plainsum  <- apply(adjmat1, 1, sum)
    squaresum <- apply(adjmat1^2, 1, sum)
    total.edge = plainsum^2 - squaresum
    CChelp=rep(-666, no.nodes)
    CChelp=ifelse(total.edge==0,0, nolinksNeighbors/total.edge)
    CChelp}
} # end of function