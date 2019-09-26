mergeCounts = function(c,GRanges){
	counts <- c
	merged <- GRanges(seqnames = seqnames(GRanges),
		ranges = ranges(GRanges),
		strand = strand(GRanges), gene_id = GRanges@elementMetadata$gene_id,
		counts = as.data.frame(counts))
      	names(merged) <- GRanges@elementMetadata$gene_id
      	merged
}
