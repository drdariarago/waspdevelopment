## Generate a list of geneIDs per each cluster
# Initialize script
library(stringr)
library(plyr)
newdir<-file.path(getwd(), "./Output/WGCNA_GO_input")
dir.create(newdir)

# Load dataset
clusterassignments_simple <- read.csv("./Output/WGCNA_clustering_biweight/clusterassignments_simple.csv")

# Create gene IDs
clusterassignments_simple$geneID<-str_extract(string = clusterassignments_simple$transcriptID, pattern = "^[^_]*")

# Specify first transcript (only for preliminary GO enrichment)
clusterassignments_simple$geneID<-paste(clusterassignments_simple$geneID, "t1", sep = "")

# Generate one list per each cluster
clusterassignments_list<-dlply(.data = clusterassignments_simple, .variables = .(clusterID), .fun = summarize, geneID=geneID)
clusterassignments_list<-lapply(X = clusterassignments_list, FUN = unique)
clusternames<-names(clusterassignments_list)

# Save gene IDs from each list in their separate files

mapply(
  x=clusterassignments_list,
  y=clusternames,
  FUN = function(x,y){
    write(unlist(x), file = file.path(newdir, paste(y, ".txt")))
  }
)

# Save reference gene list (same gene IDs, all genes included in each cluster)

write(unlist(clusterassignments_list), file=file.path(newdir, "OGS2_cluster_refset.txt"))
