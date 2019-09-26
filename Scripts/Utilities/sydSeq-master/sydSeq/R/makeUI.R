makeUI = function( simpleGRanges, gFeature ){
require(Biobase)
cat("Loading data from database...\n")
 	GR <- exons( gFeature ,columns = c('tx_name','exon_id','exon_name','gene_id'))
geneID <- IRanges::unlist(elementMetadata(GR)$gene_id)
names(geneID) <- IRanges::unlist(elementMetadata(GR)$exon_name)
exonsG <- IRanges::unlist(exonsBy(gFeature, by=c("tx"),use.names=TRUE))
transcriptID <- substr(names(exonsG[,]),1,nchar(strsplit(names((exonsG)[,]),'.1')[[1]]))
names(transcriptID) <- elementMetadata(exonsG)$exon_name

tranSplit <- split(transcriptID,geneID[names(transcriptID)])
noTran <- unlist(lapply(tranSplit,function(x)length(unique(x))))

findUI = function(x){
useSimple <- seqnames((simpleGRanges)) == x
useGenomic <- seqnames(exonsG) == x
geneNames <- elementMetadata(simpleGRanges[useSimple])$gene_id
simpleGRanges[useSimple][countOverlaps(ranges(simpleGRanges[useSimple]),ranges(exonsG[useGenomic]))==noTran[geneNames]]
}

chrNames <- as.character(unique(seqnames(simpleGRanges)))
names(chrNames) <- chrNames
chrNames <- as.list(chrNames)
UIanno <- lapply(chrNames,findUI)
eval(parse(text=(paste('c(',paste(paste('UIanno[[',1:length(UIanno),']]',sep = ''),collapse=','),')',sep = ''))))
}

