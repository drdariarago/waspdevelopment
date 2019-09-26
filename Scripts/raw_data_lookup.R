## Raw data lookup
# Load data
nasdevexon <- read.csv(file = "./Output/exon_to_transcript_clustering_rankingonly_adjusted/nasdevexon_prefiltering.csv")
exonass <- read.csv(file = "./Output/exon_to_transcript_clustering_rankingonly_adjusted/eigenexons_assignments.csv")[,-1]

# Reshape as descriptors~expr
a<-melt(nasdevexon, id.vars = c("X","geneID","exonID","ExonID"))
# add sex annotation
a$sex <- ifelse(grepl(pattern = "female", x = a$variable),"female","male")
# add stage
a$stage <- factor(x = substr(a$variable, start = 0, stop = 5), levels = c("emb10","emb18","lar51","pupye","adult"))
# add assignments
a$exonID <- a$X
b <- merge(exonass[,c("eigenexonID","exonID")], a[,-1], by = "exonID", all.y = F)
b <- droplevels(b)
# add eigenexon assignments


# plot (Transformer/feminizer)
ggplot(data = a[which(a$geneID=="Nasvi2EG005321"),], mapping = aes(x = stage, y = value, col = sex)) + geom_point() + geom_smooth(aes(group = sex))
ggplot(data = a[which(a$geneID=="Nasvi2EG005322"),], mapping = aes(x = stage, y = value, col = sex)) + geom_point() + geom_smooth(aes(group = sex))
# DSX
# ggplot(data = a[which(a$geneID=="Nasvi2EG010976"),], mapping = aes(x = stage, y = value, col = sex)) + geom_point() + geom_smooth(aes(group = sex)) + facet_grid(eigenexonID~.)
ggplot(data = b[which(b$geneID=="Nasvi2EG010980"),], mapping = aes(x = stage, y = value, col = sex)) + geom_point() + geom_smooth(aes(group = sex)) + facet_grid(eigenexonID~.)
# Sex-lethal
ggplot(data = b[which(b$geneID=="Nasvi2EG000104"),], mapping = aes(x = stage, y = value, col = sex)) + geom_point() + geom_smooth(aes(group = sex)) + facet_grid(eigenexonID~.)


# generic function
ggplot(data = b[which(b$geneID=="Nasvi2EG017727"),], mapping = aes(x = stage, y = value, col = sex)) + geom_point() + geom_smooth(aes(group = sex)) + facet_grid(eigenexonID~.)
