 ### Unify main TOM
# Initialize script
newdir<-file.path(getwd(), "Output/dTOMmat")
dir.create(newdir)
# load matrix list
load(file="./Output/WGCNA_clustering_biweight/dTOM")
# Create empty dTOM with appropriate dimensions
dTOMmat<-matrix(c(1),nrow=(nrow(dTOM[[1]])+nrow(dTOM[[2]])),ncol=(ncol(dTOM[[1]])+ncol(dTOM[[2]])))
# add matrices to main matrix
dTOMmat[1:nrow(dTOM[[1]]),1:ncol(dTOM[[1]])]<-dTOM[[1]]
dTOMmat[(nrow(dTOM[[1]])+1):nrow(dTOMmat),(ncol(dTOM[[1]])+1):ncol(dTOMmat)]<-dTOM[[2]]
# Add names to matrices
colnames(dTOMmat)<-c(colnames(dTOM[[1]]),colnames(dTOM[[2]]))
rownames(dTOMmat)<-c(rownames(dTOM[[1]]),rownames(dTOM[[2]]))
# save main dTOM
save(dTOMmat, file = file.path(newdir, "dTOMmat"), compress = T)
# save reduced dTOM for testing
dTOMmat_test<-dTOMmat[1:1000,1:1000]
save(dTOMmat_test, file=file.path(newdir, "dTOMmat_test"), compress = T)