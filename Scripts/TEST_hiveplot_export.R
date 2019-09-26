## export sample adjacency matrix as network
# load libraries and data
library(ndtv)
testnet<-adj[[1]][1:100,1:100]
testnet<-apply(X = testnet, MARGIN = c(1,2), FUN = function(x){ifelse(x<quantile(testnet, .5), 0,1)})
row.names(testnet)<-1:nrow(testnet)
# convert matrix to network
testnet<-as.network.matrix(x = testnet, matrix.type = "adjacency", directed = F, ignore.eval = T)
# save network as .dot file
export.dot(testnet,file="./Output/example_adj_mat.dot")


testnet<-adj2HPD(M = testnet, type = "2D")

# most options in the HPD object are columns in the hpd$nodes and hpd$edges data.frames, which are also indexed by node IDs (id of the two nodes for edges)
# we can calculate the network parameters in WGCNA and then either paste then in the HPD object via merge, or coerce te data.frame into an HPD object
# the remaining parameters are desc (description, character) axis cols (color of the axes, character) and type (2D or 3D, character)

# add degree (radius) to each node
testnet<-mineHPD(testnet, option = "rad <- tot.edge.count")
# assign to axis based on degree of node
testnet$nodes$axis<-as.integer(sapply(testnet$nodes$radius, function(x){if(x<25){1}else if (x<50){2}else {3}}))
