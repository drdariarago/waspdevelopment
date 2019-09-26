pMim = function(DataMi,DataG,classes,targets,pathways,Zmi=NULL,Zg=NULL,stat = 'cor',corP = 'Fisher',nInt = 2,randomTarg=FALSE, randomPath=FALSE,verbose =TRUE){

if(verbose == TRUE){
  if(sum(colnames(DataG) != names(classes))>0 |is.null(names(classes))){
    cat('colnames of Gene Data do not match names of classes')
  }
  }
if(verbose == TRUE){
  if(sum(colnames(DataG) != names(classes))>0 |is.null(names(classes))){
    cat('colnames of miRNA Data do not match names of classes')
  }
}

if(verbose == TRUE){
  if(sum(rowMeans(DataG)<5)>0){
    cat('You might want to consider filtering out genes with small average expression from your data')
  }
  }

if(verbose == TRUE){
  if(sum(rowMeans(DataMi)<5)>0){
    cat('You might want to consider filtering out miRNA with small average expression from your data')
  }
}

  y = names(classes)
  X = t(DataG)
  Y = t(DataMi)
  

#Estimate correlations
if(stat=='cor'){
  corYX = cor(Y[y,],X[y,])
  n = length(y)
  
#Estimate "probability" of correlations.

  if(corP=='estimate'){
    TR = corYX
    TR[] = ((rank(corYX)-0.5)/length(corYX))
    TR = qnorm(TR)
  }else{
    if(n>3){
      TR <- sqrt(n-3)*atanh(corYX)
    }else{
      cat('Need more replicates or use cor="estimate"')
      break
    }
  }
}  

#Estimate z-scores
if(is.null(Zmi)){
  require(limma)
  design = model.matrix(~classes)
  Dat = DataMi[,names(classes)]
  fit = lmFit(Dat,design)
  ordinary.t <- fit$coef / fit$stdev.unscaled / pmax(pmax(fit$sigma,sqrt(rowMeans(Dat))),1)
  Tmi = ordinary.t[,2]
  Zmi = qnorm(pt((Tmi),fit$df.res)) 
}

if(stat=='de'){
  if(is.null(Zmi)){
    require(limma)
    design = model.matrix(~classes)
    Dat = DataG[,names(classes)]
    fit = lmFit(Dat,design)
    ordinary.t <- fit$coef / fit$stdev.unscaled / pmax(pmax(fit$sigma,sqrt(rowMeans(Dat))),1)
    Tg = ordinary.t[,2]
    Zg = qnorm(pt((Tg),fit$df.res)) 
  }
}


  #Create mapMat and pathMat, matrix versions of the target and pathways

  mi = colnames(Y)
  msgr = intersect(unlist(pathways),unlist(targets))
  mapMat = matrix(0,length(mi),length(msgr))
  rownames(mapMat) = mi
  colnames(mapMat) = msgr
   for(i in rownames(mapMat)){
    if(randomTarg==TRUE){mapMat[i,sample(colnames(mapMat),sum(msgr%in%targets[[i]]))] = 1
    }else{mapMat[i,msgr%in%targets[[i]]] = 1}
    
  }
  use = names(which(rowSums(mapMat)>0))
  mapMat =mapMat[use,]
  mi = use
  
  pathname = names(pathways)
  pathMat = matrix(0,length(pathname),length(msgr))
  rownames(pathMat) = pathname
  colnames(pathMat) = msgr
   for(j in 1:length(pathname)){
    if(randomPath==TRUE){pathMat[j,sample(colnames(mapMat),sum(msgr%in%pathways[[j]]))] = 1
    }else{pathMat[j,msgr%in%pathways[[j]]] = 1}
  }  
  
  N = mapMat%*%t(pathMat) 
  
  test = which(N>=nInt,2)
  
  
  
  G = matrix(NA, length(mi),length(pathname))
  rownames(G) = mi
  colnames(G) = pathname
  
  Tmat = corYX
  Tmat[] = NA
  for(k in 1:nrow(test)){
    i = mi[test[k,1]]
    j = test[k,2]
    pathway = pathMat[j,]
    tr = TR[i,]
    binding = mapMat[i,]
    
   
    if(tolower(stat) == 'de'){t = Zg[names(TR[i,])]}
    if(tolower(stat) == 'cor'){t = TR[i,]}
    
    Tmat[i,]=t
    
    g = sum(t[names(which(binding*pathway==1))])/sqrt(sum(binding*pathway))
    
    G[i,j] = pnorm(g)    
    
  }
  
  G = t(apply(G,1,function(x){x[!is.na(x)] = p.adjust(x[!is.na(x)],'fdr');return(x)}))
  Score = pchisq(-2*(log(1-G) + log(1-matrix(pnorm(-abs(Zmi))[rownames(G)],dim(G)[1],dim(G)[2],byrow = FALSE))),4)
  
pval= Score

updown = matrix(sign(Zmi)[rownames(G)],dim(G)[1],dim(G)[2],byrow = FALSE)
updown[updown ==1] = 'Up'
updown[updown ==-1] = 'Down'
use = which(!is.na(pval),2)

res = cbind(rownames(use),updown[use],pval[use],colnames(pval)[use[,2]])
res = res[order(res[,3]),]

res2 = NULL
for(i in unique(res[,3])){
  use = res[,3]==i 
  pathways = paste(res[use,4],collapse = ', ')
  miRNA = res[use,1][1]
  direction = res[use,2][1]
  res2 = rbind(res2,c(miRNA,direction,signif(as.numeric(i),2),pathways))
}


colnames(res2) = c('miRNA','direction', 'Score', 'Pathways')


  list(Results = res2,Scores = Score)
}
