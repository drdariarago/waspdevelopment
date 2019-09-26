## Comparison between adult-based sex-bias and developmental sex-bias
# Run Basestats

# What is the proportion of convergent genes in each set?
prop.table(table(dev=Fdrfitall[,10]<0.0001, adult=Fdrfitallgonads[,1]<0.0001))
# What is the correlation between the two sets? Here we use the same permissive FDR threshold
cor(Fdrfitall[,10]<0.0001, Fdrfitallgonads[,1]<0.0001, method="spearman")
# Note: increasing stringency to 1 False Discovery in each set yelds better correlation
cor(Fdrfitall[,10]<0.00001, Fdrfitallgonads[,1]<0.0001, method="spearman")
# Now comparing correlation of q-values (almost perfect)
cor(Fdrfitall[,10], Fdrfitallgonads[,1], method="spearman")
# plotting correlation
plot(log10(Fdrfitall[,10]), log10(Fdrfitallgonads[,1]), xlab="Developmental set q-values", ylab="Adult-Gonad set q-values", 
     col=as.factor((Fdrfitall[,10]<0.00001)+2*(Fdrfitallgonads[,1]<0.0001)),
     main="Comparison between sex-specific genes in adults\n between the Developmental and the Gonadal set"
     )
