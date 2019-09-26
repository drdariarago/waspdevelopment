simplifyGFeatures = function( gFeature ){
require(GenomicFeatures)
require(GenomicRanges)
	cat("Loading data from database...\n")
 	GR <- exons( gFeature ,columns = c('exon_id','gene_id'))
	cat("Simplifying annotation...\n")
	ranges <- ranges(GR)
	names(ranges) <- IRanges::unlist(elementMetadata(GR)$gene_id)
 	chrList <- IRanges::split(ranges,as.character(seqnames(GR)))
	simpleList <- IRanges::lapply(chrList,simplifyIRanges)
	simpleAnno <- NULL
	for(i in 1:length(simpleList)){
		simpleAnno = rbind(cbind(ensembl_gene_id=names(simpleList[[i]]), start = as.numeric(simpleList[[i]]@start), end = as.numeric(simpleList[[i]]@start + simpleList[[i]]@width -1),chr = names(simpleList[i]), strand = '*'),simpleAnno)
		}
	simpleRanges <- GRanges(seqnames = simpleAnno[,'chr'],
		ranges = IRanges(start = as.numeric(simpleAnno[,'start']),
		end = as.numeric(simpleAnno[,'end'])),
		strand = simpleAnno[,'strand'], gene_id = simpleAnno[,'ensembl_gene_id'])
      names(simpleRanges) = simpleAnno[,'ensembl_gene_id']
	simpleRanges
	}


simplifyIRanges = function(iRange){
disjointRanges <- IRanges::disjoin(iRange)
un <- unique(names(iRange))
nam <- rep('non',length(disjointRanges))
for(i in 1:length(un)){
use = which(!is.na(IRanges::findOverlaps(disjointRanges,iRange[which(names(iRange)%in%un[i])],select = 'first')))
nam[use][nam[use]!= 'non'] <- 'multi'
nam[use][nam[use]== 'non'] <- as.character(un[i])
}
names(disjointRanges) <- nam
simpleRanges <- disjointRanges[names(disjointRanges)!='non'&names(disjointRanges)!='multi']
simpleRanges
}
