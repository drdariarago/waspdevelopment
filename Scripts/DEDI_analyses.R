## Check clusters with both DE and DI

DEDI <- clusterdata[which(clusterdata$cluster_devsexbias!="....."&is.na(clusterdata$DiffIntegrated)==F),]
# Remove only MF cluster
DEDI <- DEDI[-grep(pattern = "mf", x = DEDI$cluster_devsexbias),]

# Count convergent vs divergent for mDE fDI, fDE mDI
table(grepl(pattern = "m", x = DEDI$cluster_devsexbias), grepl(pattern = "Female", x = DEDI$DiffIntegrated))
table(grepl(pattern = "f", x = DEDI$cluster_devsexbias), grepl(pattern = "Male", x = DEDI$DiffIntegrated))


sum(table(grepl(pattern = "f", x = DEDI$cluster_devsexbias), grepl(pattern = "Male", x = DEDI$DiffIntegrated)))
fisher.test(table(grepl(pattern = "f", x = DEDI$cluster_devsexbias), grepl(pattern = "Male", x = DEDI$DiffIntegrated)))

DEDI$divergent <- xor(grepl(pattern = "f", x = DEDI$cluster_devsexbias),grepl(pattern = "Male", x = DEDI$DiffIntegrated))==F
