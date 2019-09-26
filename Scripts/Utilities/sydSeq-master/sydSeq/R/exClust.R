exClust = function(simpleGRanges,cut= 1.5){
	geneList <- IRanges::split(simpleGRanges,as.character(elementMetadata(simpleGRanges)$gene_id))
	simpleList <- IRanges::lapply(geneList,Clust,cut)

anno = eval(parse(text=(paste('c(',paste(paste('simpleList[[',1:length(simpleList),']]',sep = ''),collapse=','),')',sep = ''))))
anno
}


Clust =function(simpleGRanges, cut) {  
  
  major_group <- function(grouping) {
    group_id = names(table(grouping))[table(grouping) == 
                                        max(table(grouping))][1]
    return(names(grouping)[grouping == group_id])
  }
  counts = (elementMetadata(simpleGRanges)[, grep("counts", 
                                                  IRanges::colnames(elementMetadata(simpleGRanges)))])
  if (is.null(counts)) {
    cat("need counts matrix/n")
    break
  }
  IRanges::colnames(counts) = gsub("counts.", "", IRanges::colnames(counts))
  counts = IRanges::as.data.frame(counts)
  ncount = dim(counts)[2]
  counts = as.matrix(counts, ncol = ncount)
  NROW = nrow(counts)
  rownames(counts) = 1:nrow(counts)
  ir = ranges(simpleGRanges)
  width = width(ir)
  if (dim(counts)[1] > 2 & sum(rowSums(counts) > 5) > 1) {
    x.gs = colSums(counts)
    x.ge = rowSums(counts)
    x.g = sum(counts)
    MU = x.ge %*% t(x.gs)/x.g
    MU[MU < 0] = 0
    Z = (2 * sqrt(counts + 3/8) - 2 * sqrt(MU + 3/8))
    corZ = (Z[, ] %*% t(Z[, ]))/(dim(Z)[2] - 1)
    width = width[rowSums(counts) != 0]
    Z = Z[rowSums(counts) != 0, ]
    corZ = corZ[rowSums(counts) != 0, rowSums(counts) != 
                  0]
    counts = counts[rowSums(counts) != 0, ]
    probe_corr = corZ
    probe_corr[is.na(probe_corr)] = 0
    probe_corr_dist <- as.dist(1 - probe_corr)
    probe_tree <- hclust(probe_corr_dist, method = "complete")
    cut_trees = cutree(probe_tree, h = cut)
    countsW = counts/matrix(width, nrow(counts), ncol(counts))
    big = rowMeans(apply(countsW, 2, rank) * (rowSums(counts) > 
                                                0) * (width > 10))
    usetree = cut_trees
    tree = usetree
    bigsum2 = function(tree, big) {
      lu = length(unique(tree))
      SS = rep(0, lu)
      for (i in 1:lu) {
        bb = big[tree == i]
        bb = bb[bb != 0]
        SS[i] = mean(head(-sort(-bb), 2))
      }
      which(SS == max(SS,na.rm = TRUE))
    }
    bs = bigsum2(usetree, big)
    major_cluster <- names(usetree[usetree %in% bs])
    CC = simpleGRanges[(1:NROW) %in% major_cluster]
  }
  else CC = simpleGRanges
  CC
}