TshrinkPlus = function(Data, classes, meanPlus = NULL, variancePlus = NULL, verbose = TRUE,nn = 0.7){

if(verbose == TRUE){
  if(sum(colnames(Data) != names(classes))>0 |is.null(names(classes))){
    cat('colnames of Data do not match names of classes')
  }
  }

if(verbose == TRUE){
  if(sum(rowMeans(Data)<5)>0){
    cat('You might want to consider filtering out genes with small average expression from your data')
  }
  }


if(verbose == TRUE){
  if(sum((meanPlus)<5)>0&sum((meanPlus)>100)>0){
    cat('You might want to consider filtering out genes with small average expression from your additional data')
  }
}

int = rownames(Data)
if(!is.null(meanPlus)){
varStandard = StandardiseVariance(meanPlus,variancePlus)
int = intersect(rownames(Data), names(which(!is.na(varStandard ))))
}

if(is.null(meanPlus)&!is.null(variancePlus)){
  varStandard = variancePlus
  int = intersect(rownames(Data), names(which(!is.na(varStandard ))))
}

design = matrix(0,length(classes),length(unique(classes)))
design[,1] = 1
for(i in 2:length(unique(classes))){
design[classes == unique(classes)[i],i] = 1
}

if(length(unique(classes))>2){cat('Please use only two unique class labels... pretty please');break}

udes = unique(design)
designUniq = design
for(i in 1:nrow(udes)){
  designUniq[,i] = rowSums(design==matrix(udes[i,],dim(design)[1],dim(design)[2],byrow = TRUE))==ncol(design)
}
colUse = colSums(designUniq)>1
if(sum(colUse)==0)stop("Need more replication or simpler design")

if(is.null(meanPlus)){
  
  PREDVAR = as.matrix(Data)%*%designUniq
  VARIANCE = as.matrix(Data)%*%designUniq
  for(i in which(colUse)){
    rowM = rowMeans(Data[,designUniq[,i]==1])
    variance = apply(Data[,designUniq[,i]==1],1,var)
    VARIANCE[,i] = variance
    avgExpr = rowM
    predVar = avgExpr
    u = avgExpr>0
    fit = locfit((variance[u])~lp(log(avgExpr[u]),nn = nn),family = 'gamma')
    predVar4 = (fitted( fit )) 
    predVar[u] = predVar4
    predVar1 = pmax(predVar,avgExpr)
    PREDVAR[,i] = predVar1
  }
      
}


if(!is.null(meanPlus)){
  
  
  PREDVAR = as.matrix(Data[int,])%*%designUniq
  VARIANCE = as.matrix(Data[int,])%*%designUniq
  for(i in which(colUse)){
    rowM = rowMeans(Data[int,designUniq[,i]==1])
    variance = apply(Data[int,designUniq[,i]==1],1,var)
    VARIANCE[,i] = variance
    avgExpr = rowM
    predVar = avgExpr
    u = avgExpr>0
     qV = varStandard [int][u]
  
    fit = locfit((variance[u])~lp(log(avgExpr[u]),qV,nn = nn,scale = TRUE),maxk=300,family = 'gamma')
    predVar4 = (fitted( fit ))
    predVar[u] = predVar4
    predVar1 = pmax(predVar,avgExpr)
    PREDVAR[,i] = predVar1
  }
  }

u = names(which(rowMeans(Data[int,])>50))
N = matrix(colSums(designUniq),length(u),dim(designUniq)[2],byrow = TRUE)
lambda = 2/mean((N-1)*(VARIANCE[u,]/PREDVAR[u,]-1)^2,na.rm =TRUE)
N = matrix(colSums(designUniq),dim(VARIANCE)[1],2,byrow = TRUE)

VAR =  ((rowSums(PREDVAR*(N-1))/(sum(designUniq)-2))*lambda+ (rowSums(VARIANCE*(N-1))/(sum(designUniq)-2))*(1-lambda))



n2 = Inf
n1 = sum(designUniq)-2
df = ((1)*lambda+ ((rowSums(VARIANCE*(N-1))/(sum(designUniq)-2))/(rowSums(PREDVAR*(N-1))/(sum(designUniq)-2)))*(1-lambda))^2/(((rowSums(VARIANCE*(N-1))/(sum(designUniq)-2))/(rowSums(PREDVAR*(N-1))/(sum(designUniq)-2)))^2*(1-lambda)^2/n1+ lambda^2/n2)

df = df[!is.na(df)]
df = mean(df)

T = (rowSums(Data[int,classes == unique(classes)[1]]) - rowSums(Data[int,classes == unique(classes)[2]]))^2/(VAR*length(classes))
P = 2*pt(sqrt(T),df,lower.tail = FALSE)
P = P[!is.na(P)]

PP = rep(NA,dim(Data)[1])
names(PP) == rownames(Data)
PP[int] = P
pvalue = PP
list(pvalue = pvalue,df = df,lambda = lambda )

}

StandardiseVariance = function(Mean,Variance,nn = 0.7){
  use = Mean>0
  if(max(Mean)>30&min(Mean)>=0)Mean = log(Mean)
  # use = Mean>0#&Variance>0
  VARS = Variance[use]
  AVG = Mean[use]
  fit = locfit(VARS~lp(AVG,nn = nn),family = 'gamma')
  VARS = log(VARS/fitted(fit))
  Variance[!use] = NA
  Variance[use] = VARS
  Variance
}

