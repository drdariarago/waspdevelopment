

KODAMA = function(data,
                  n         = 100,
                  m         = 20,
                  FUN_VAR   = function(x) {ceiling(ncol(x))},
                  FUN_SAM   = function(x) {ceiling(nrow(x)*0.75)},
                  bagging   = FALSE,
                  FUN       = KNN.CV,
                  f.par     = list(kn=10),
                  verbose   = TRUE,
                  startv    = NULL,
                  constrain = NULL,
                  fix       = rep(F,nrow(data)),
                  epsilon   = 0.05,
                  shake     = TRUE){


   require("e1071")  
   require("entropy")
   require("impute")
   data = as.matrix(data)
   nva=ncol(data)
   nsa=nrow(data)
   ma=matrix(0,ncol=nsa,nrow=nsa)
   normalization=matrix(0,ncol=nsa,nrow=nsa)   
   FUN_VAR=FUN_VAR(data)
   FUN_SAM=FUN_SAM(data)
   vect_acc=matrix(nrow=n,ncol=m)
   accu=NULL
   #########
   whF=which(!fix)
   whT=which(fix)   
   FUN_SAM=FUN_SAM-length(whT)
   
   data=t(impute.knn(t(data),rowmax = 1, colmax = 1)$data)


   for(k in 1:n){
      if(verbose==T)   print(k)
      
      
        sva=.Internal(sample(nva,FUN_VAR, FALSE, NULL))      
        ssa=c(whT,sample(whF,FUN_SAM, bagging, NULL))

        x=data[ssa,sva]

   
      xva=ncol(x)
      xsa=nrow(x)
      if(is.vector(startv))    clbest=startv[ssa]
      if(is.function(startv))  {clbest=startv(x)}else{ 
      if(is.null(startv)) {    
      	      

         clbest=1:xsa
          
     
      }else{ 
         if(any(is.na(startv))) {
            unw=unique(startv)
            ull=length(unw)
            unw=unw[-which(is.na(unw))]     
            ghg=is.na(startv)
            startv[ghg]=(1:(sum(ghg)+length(unw)))[-unw]
            clbest=startv[ssa]
         }
      }
      }

      
      
      yatta=core(x,clbest,m=m,FUN=FUN,f.par=f.par,constrain=constrain[ssa],fix=fix[ssa],shake)
      clbest=yatta$c
      accu=yatta$a
      knk=length(yatta$v)
      vect_acc[k,1:knk]=yatta$v
      
      uni=unique(clbest)
      nun=length(uni)
      for(ii in 1: nun)
         ma[ssa[clbest==uni[ii]],ssa[clbest==uni[ii]]] = ma[ssa[clbest==uni[ii]],ssa[clbest==uni[ii]]]+1
      normalization[ssa,ssa]=normalization[ssa,ssa]+1
   }

#   ma=ma/diag(ma)

   ma=ma/normalization

   Edist=as.matrix(dist(data))

   
   ma[ma<epsilon]=0

   mam=(1/ma)*Edist

   mam=allShortestPaths(mam)$length

   prox=Edist/mam
   diag(prox)=1
   prox[is.na(prox)]=0
   me=entropy(prox)
   prox[is.na(prox)]=0
   mam2=mam;mam2[is.na(mam2)]=max(mam2,na.rm=T)
   return(list(kodama=mam,acc=accu,entropy=me,prox=prox,dis=mam2,v=vect_acc))
}

#############################################################                          

core = function(x,              #matrix
                clbest,         #starting vector
                m=20,           #number of cycles
                FUN=KNN.CV,
                f.par=list(kn=10),
                constrain=NULL,
                fix=NULL,
                shake=TRUE)
{
   if(is.null(constrain))   constrain=1:length(clbest)
   if(is.null(fix))         fix=rep(FALSE,length(clbest))
   xsa=nrow(x)
   nconc=length(unique(constrain))   
   clbest=as.factor(clbest)  
   cvpred=clbest
   cvpredbest=do.call(FUN, c(list(x,clbest,constrain),f.par))    
   levels(cvpredbest)=levels(clbest)  

   if(shake==TRUE){
      accbest=mean(as.numeric(clbest==cvpredbest))   
   }else{
      accbest=0	
   }
   success=FALSE
   j=0
   vect_acc=NULL
   while(j<m & !success){
      j=j+1
      cl=clbest
      temp=.Internal(sample(nconc, 1, FALSE, NULL))
      ss=.Internal(sample(nconc, temp, FALSE, NULL))
      for(k in ss) {       
         sele=(constrain==k & !fix)     
         if(sum(sele!=0)) {
            soso=sample(unique(cvpredbest[sele]),1)
            cl[sele]=soso     
         } 
      }
      cvpred=do.call(FUN, c(list(x,(cl),constrain),f.par))
      levels(cvpred)=levels(cl)
      accTOT=mean(as.numeric(cl==cvpred))
      if(accTOT>accbest){
         cvpredbest=cvpred
         clbest=cl       
         accbest=accTOT
      }
      vect_acc[j]=accbest
      if(accTOT==1)   success=TRUE
   }
   return(list(c=clbest,a=accbest,v=vect_acc))
}

#############################################################                          

PCA.CA.KNN.CV = function(x,cl,constrain,kn=10,variance=0.9){
   require("e1071")
  # require("knnflex")
   cl = as.factor(cl)
   uni= unique(cl)
   xsa=nrow(x)
   if(length(uni)==1) return(rep(uni,xsa))
   pr=prcomp(x)
   va=pr$sdev^2
   va=va/sum(va)
   cum=which(cumsum(va)<variance)
   if(length(cum)<2) cum=1:2
   x=pr$x[,cum]
   ncomp=ncol(x)
   ca.out <- cancor(x, transformy(cl))
   nn <- length(ca.out$cor)
   mm <- ncol(ca.out$xcoef)
   if (nn == 1 && mm > 1)  nn <- nn + 1
   if (mm < ncomp)   ca.out$xcoef = rbind(ca.out$xcoef,matrix(0,ncol=mm,nrow=ncomp-mm))
   x <- x %*% ca.out$xcoef[,1:nn,drop=F]
   cvpred=cl
   fold = kfold(constrain) 
   xdist=knn.dist(x)
   for(i in 1:10){
      ff=fold==i
      w1=which(ff)
      w9=which(!ff)
      cvpred[w1]=knn.predict(train=w9,test=w1,cl,xdist,k=kn)
   }
   cvpred
}


#############################################################    

KNN.CV = function(x,cl,constrain,kn=10){
  # require("knnflex")
   xsa=nrow(x)
   uni = unique(cl)
   cl=as.factor(cl)
   if(length(uni)==1) return(rep(uni,xsa))
   cvpred = cl
   fold = kfold(constrain)  
   xdist=knn.dist(x)
   for(i in 1:10){
      ff=fold==i
      w1=which(ff)
      w9=which(!ff)
      cvpred[w1]=knn.predict(train=w9,test=w1,cl,xdist,k=kn)
   }
   cvpred
}

#############################################################                          
 
PLS.SVM.CV = function(x,cl,constrain,ncomp=5,...){
   require("e1071")
   require("plsgenomics")
   xsa=nrow(x)
   pls.svm.predict = function (train,test,cl,data,ncomp=5, ...){
      Xtrain=data[train,]
      Xtest=data[test,,drop=F]
      Ytrain=cl[train]
      uni=unique(Ytrain)
      if(length(uni)==1) return(rep(uni,length(test)))
      ntrain <- nrow(Xtrain)
      nn <- min(dim(Xtrain))-1
      if ( ncomp == 0 | ncomp > nn) ncomp <- nn
      Xtrain <- scale(Xtrain, center=T, scale=F)
      Xtest <-  scale(Xtest, center=attr(Xtrain,"scaled:center"), scale=F)
      red.out = pls.regression(Xtrain = Xtrain, transformy(Ytrain), Xtest = NULL, ncomp = ncomp)
      new.train <- matrix(red.out$T[,1:ncomp], nrow=ntrain, ncol=ncomp)
      new.test <- Xtest %*% red.out$R[,1:ncomp]

      cl.out <- svm(x=new.train, y=Ytrain,type="C-classification", ... )
      predict(object = cl.out, newdata = new.test)
   }
   cvpred=NULL
   fold = kfold(constrain)
   for(i in 1:10){
      ff=fold==i
      w1=which(ff)
      w9=which(!ff)
      if(length(unique(cl[w9]))>1){
         cvpred[w1]=pls.svm.predict(train=w9,test=w1,cl,x,ncomp=ncomp,...) 
      }else{
         cvpred[w1]=cl[w1]
      }
      
   }
   as.factor(cvpred)
} 




###################################

kfold = function(constrain,k=10){
   constrain=as.numeric(as.factor(constrain))
   xsa=max(constrain)
   v=sample(1:xsa)
   (v[constrain]%%k+1)
}
#############################################################

transformy = function (y) {
   y  =  as.numeric (as.factor( y ))
   n  = length(y)
   nc = max(y,na.rm=T)
   Y  = matrix(0, n, nc)
   for(k in 1:nc) 
      Y[, k] <- as.numeric(y == k)
   Y
}

#############################################################


gauss = function(n=c(100,100,100),dims=2){
   require("Matrix")
   require("MASS")
   clusters=length(n)
   
   data=NULL
   zz=sample(c(10,20,30,50,100),1)
   for(i in 1:clusters){
      Sigma <- as.matrix(as.dist(matrix(runif(dims^2)*2-1,dims,dims)))
      diag(Sigma)=1
      Sigma=nearPD(Sigma)$mat
      data=rbind(data,mvrnorm(n=n[i], runif(dims)*(zz/sqrt(dims)), Sigma))
   }
   data
}

spirals = function (n=c(100,100,100),sd=c(0,0,0)){
   clusters=length(n)
   x=NULL
   y=NULL
   for(i in 1:clusters){
      t=seq(1/(4*pi),1,length.out=n[i])^0.5*2*pi
      a=rnorm(n[i],sd=sd[i])
      x=c(x,cos(t+(2*pi*i)/clusters)*(t+a))
      y=c(y,sin(t+(2*pi*i)/clusters)*(t+a))
   }
   cbind(x,y)   
}



Kindex=function(u,pp,knn=5){
   dd=as.matrix(dist(u))
   n=nrow(u)
   v=NULL
   for(i in 1:n)
      v[i]=sum(sort(dd[i,which(pp[i]==pp)])[c(2:(knn+1))])<sum(sort(dd[i,-which(pp[i]==pp)])[c(1:knn)])
   mean(v)
}


mean.na = function (x, ...) 
{
    mean(x[!(is.na(x) | is.infinite(x))])
}

order.na = function (x, na.last = TRUE) 
{
    y <- order(x)
    n <- sum(is.na(x))
    tmp <- (length(x) - n + 1):length(x)
    if (!is.na(na.last)) {
        if (na.last) 
            res <- y
        if (!na.last) {
            if (n == 0) 
                res <- y
            else res <- c(y[tmp], y[-tmp])
        }
    }
    if (is.na(na.last)) {
        warning("NA's discarded")
        res <- y[-tmp]
    }
    res
}

sum.na=function (x, ...) 
{
    res <- NA
    tmp <- !(is.na(x) | is.infinite(x))
    if (sum(tmp) > 0) 
        res <- sum(x[tmp])
    res
}

var.na =function (x) 
{
    res <- NA
    tmp <- !(is.na(x) | is.infinite(x))
    if (sum(tmp) > 1) 
        res <- var(x[tmp])
    res
}

###################################################










#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
# KNNDIST
# majority.R

`majority` <- function(x){
  x <- as.factor(x)
  n <- nlevels(x)
  votes <- rep(0,n)
  for (i in 1:length(x)) votes[as.integer(x[i])] <- votes[as.integer(x[i])]+1
  levels(x)[order(votes,decreasing=TRUE,sample(1:n,n))[1]]
  }

#########################################################

# knn.probability.R

`knn.probability` <-
function(train, test, y, dist.matrix, k=1, ties.meth="min") {

#number of predictions to make
n<-length(test)

#sort the indexes for the training and test sets
if (is.unsorted(train)) train<-sort(train)
if (is.unsorted(test)) test<-sort(test)

#only need the rows for the test data and columns
#for the training data
d<-dist.matrix[test,train]

#ensure y is a factor
y<-as.factor(y)

#only need the responses for the training data
if (length(y)>length(train)) y<-y[train]

#calculate closest neighbors and
#return aggregate response for the k closest neighbors
if (n==1) {
  d<-rank(d, ties.method = ties.meth)
  x<-classprob(y[d <= k])
  x<-data.frame(x)
  names(x)<-test
  row.names(x)<-levels(y)
  return(x)
  }
else {
  d<-t(apply(d,1,function(x) rank(x,ties.method=ties.meth)))
  x<-apply(d,1,function(x) classprob(y[x<=k]))
  row.names(x)<-levels(y)
  return(x)
  }
}


#########################################################

# knn.predict.R

`knn.predict` <-
function(train, test, y, dist.matrix, k=1,
    agg.meth=if (is.factor(y)) "majority" else "mean",
    ties.meth="min") {

#number of predictions to make
n<-length(test)

#sort the indexes for the training and test sets
if (is.unsorted(train)) train<-sort(train)
if (is.unsorted(test)) test<-sort(test)

#only need the rows for the test data and columns
#for the training data
d<-dist.matrix[test,train]

#only need the responses for the training data
if (length(y)>length(train)) y<-y[train]

#calculate closest neighbors and
#return aggregate response for the k closest neighbors
if (n==1) {
  d <- rank(d, ties.method = ties.meth)
  x <- apply(data.frame(y[d <= k]), 2, agg.meth)
  names(x) <- test
  return(x)
  }
else {
  d<-t(apply(d,1,function(x) rank(x,ties.method=ties.meth)))
  apply(d,1,function(x) apply(data.frame(y[x<=k]),2,agg.meth))
  }
}


#########################################################

# knn.dist.R

`knn.dist` <-
function(x, dist.meth="euclidean", p=2) {
#create a distance matrix using all values in the data
d<-as.matrix(dist(x,dist.meth,p))
#fix for some small high persision errors
round(d,digits=15)
}

#########################################################

# classprob.R


`classprob` <- function(x){
  x <- as.factor(x)
  n <- nlevels(x)
  votes <- rep(0, n)
  for (i in 1:length(x)) votes[as.integer(x[i])] <- votes[as.integer(x[i])]+1
  votes/length(x)
  }
