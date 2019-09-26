# pick duplicates from eigenexons clustering (after merging with eigenexonE)
dup1<-eigenexonE2[duplicated(eigenexonE2$exonID),"eigenexonID"]
# pick duplicates from eigenexons clustering before network construction
dup2<-nasdevgeneexon[duplicated(nasdevgeneexon$eigenexonID),"eigenexonID"]
# compare the two
dup2%in%dup1
dup1%in%dup2
# now compare via base gene
library(stringr)
str_extract(dup2, "[[:alnum:]]*")%in%str_extract(dup1, "[[:alnum:]]*")
# all of the extra entries in the main dataset come from the same genes as the duplicated ones in the annotation dataset
str_extract(dup1, "[[:alnum:]]*")%in%str_extract(dup2, "[[:alnum:]]*")
# and vice versa

# this happens when the subranking algorithm reaches 1 bottom up. Must test more genes
# problem in Nasvi2EG001724 happens because the geneID includes the appendix b, but the eigenexonID doesn't (total of 22 genes annotated with b in OGS2, present in different locations of the genome than the non-b form)
# the loop selects rows via grep, therefore the *b genes are included in the non b form as well. Problem solved by using id that incudes a postscript (t), which separates genes in aaat and aaabt unique IDs