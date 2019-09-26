## extract embryonal sex biased nodes, paste onto eigenexon scores and compute their nearest neighbours
library(WGCNA)

# import node expression data
CCRE_nasdevgeneexon <- read.csv(file = "./Output/CCREs/CCRE_nasdevgeneexon.csv")
row.names(CCRE_nasdevgeneexon) <- CCRE_nasdevgeneexon$X
CCRE_nasdevgeneexon <- CCRE_nasdevgeneexon[,-1]

# filter nodes with early expression
devsexbiased_nodes <- read.csv(file = "./Output/sex_biased_nodes_CCRE/sexbiased_nodes_short_5E-2.csv")[,-1]
names(devsexbiased_nodes)<-c("nodeID","devsexbias")
earlynodes <- devsexbiased_nodes[grep(pattern = "^.?[m,f]", x = devsexbiased_nodes$devsexbias),"nodeID"]
earlynodes <- CCRE_nasdevgeneexon[which(row.names(CCRE_nasdevgeneexon)%in%earlynodes),]

# Import and merge cluster eigenexons
load(file = "./Output/WGCNA_CCRE_clustering_biweight/CCRE_MEList.RData")
ME_expr <- MEList$eigengenes
row.names(ME_expr) <- colnames(CCRE_nasdevgeneexon)
ME_expr <- t(ME_expr)
newnet <- rbind(ME_expr, earlynodes)

# calculate adjacency matrix and plot
newnet_mat <- bicor(x = t(newnet), use = "pairwise", quick = 0)
dim(newnet_mat)

# Which are the nearest nodes to Feminizer (Nasvi2EG005321_e2_con)?
FemRank <- nrow(newnet_mat)-rank(newnet_mat[grep(pattern = "Nasvi2EG005321_e2_con", row.names(newnet_mat)),])
FemRank <- FemRank[order(FemRank)]
FemRank[1:20]
# Closest cluster tan4 (17th rank), malebiased in adults, closest related Nasvi2EG010980 (doublesex, .ffff) and Nasvi2EG023025_e3_con (reverse transcriptase, .m...), 
# Nasvi2EG013953 (unknown putative DS repair protein f....) also interesting