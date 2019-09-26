geneCounts = function(anno){
columns = grep('counts.',IRanges::colnames(elementMetadata(anno)))
counts = IRanges::as.data.frame(elementMetadata(anno)[,columns])
geneSplit = IRanges::split(counts,as.character(elementMetadata(anno)$gene_id))
sumCounts = IRanges::lapply(geneSplit,function(x)colSums(x))
counts = do.call('rbind',sumCounts)
colnames(counts) = gsub("counts.","",colnames(counts))
as.matrix(counts)
}
