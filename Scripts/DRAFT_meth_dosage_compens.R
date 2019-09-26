### Check whether nodes that are female biased and methylated in adult females escape dosage compensation by checking whether they have 2-fold increases in abundance
# Load data
glmsexbiased <- read.csv(file="./Output/sex_biased_nodes/glms2_fdr_coef.csv")[,-1]

fem_adult_bias_met<-clusterdata[which(clusterdata$adult_female_meth_status=="Methylated"&grepl("f$", clusterdata$devsexbias)),"eigenexonID"]
fem_adult_bias_unmet<-clusterdata[which(clusterdata$adult_female_meth_status=="Unmethylated"&grepl("f$", clusterdata$devsexbias)),"eigenexonID"]
fem_adult_bias <- clusterdata[grepl("f$", clusterdata$devsexbias),"eigenexonID"]

library(lattice)
fem_adult_bias_met <- glmsexbiased[glmsexbiased$node_ID%in%fem_adult_bias_met, c("node_ID", "adult_Male_fdr", "adult_Male_coef")]
fem_adult_bias_met <- fem_adult_bias_met[which(fem_adult_bias_met$adult_Male_fdr<0.05),]

fem_adult_bias_unmet <- glmsexbiased[glmsexbiased$node_ID%in%fem_adult_bias_unmet, c("node_ID", "adult_Male_fdr", "adult_Male_coef")]
fem_adult_bias_unmet <- fem_adult_bias_unmet[which(fem_adult_bias_unmet$adult_Male_fdr<0.05),]

fem_adult_bias <- glmsexbiased[glmsexbiased$node_ID%in%fem_adult_bias, c("node_ID", "adult_Male_fdr", "adult_Male_coef")]
fem_adult_bias <- fem_adult_bias[which(fem_adult_bias$adult_Male_fdr<0.05),]

densityplot(fem_adult_bias_met$adult_Male_coef)
densityplot(fem_adult_bias_unmet$adult_Male_coef)
densityplot(fem_adult_bias$adult_Male_coef)

median(fem_adult_bias_met$adult_Male_coef)
median(fem_adult_bias_unmet$adult_Male_coef)
median(fem_adult_bias$adult_Male_coef)

## Sexbiased nodes in adults that are methylated in adult females are almost exclusively female-biased

## Compare adult expression and sex-bias of genes with and without female adult meth, expected increase in expr only for highly expr genes
# pattern opposite from expected: methylated genes are more female-biased even at low expression and their divergence with male expressed genes decreases with increasing expression (hitting upper limit?)

metunmeth <- clusterdata[grep("_con", clusterdata$eigenexonID),c("eigenexonID","geneID","adult_female_meth_status","devsexbias")]
metunmeth <- merge(metunmeth, glmsexbiased[,c("node_ID", "adult_coef","adult_Male_coef")], by.x = "eigenexonID", by.y = "node_ID")
ggplot(data = metunmeth, aes(y = adult_coef-adult_Male_coef, x = (adult_Male_coef+adult_coef)/2, col=adult_female_meth_status))+geom_point(alpha=0.3)+theme_bw()+geom_smooth()
ggplot(data = metunmeth, aes(y = adult_Male_coef, x = (adult_Male_coef+adult_coef)/2, col=adult_female_meth_status))+geom_point(alpha=0.3)+theme_bw()+geom_smooth()
ggplot(data = metunmeth[grep("f$", metunmeth$devsexbias),], aes(x = adult_coef, y = adult_Male_coef, col=adult_female_meth_status))+geom_point()+theme_bw()+geom_smooth()
ggplot(data = metunmeth[grep("[m|f]$", metunmeth$devsexbias),], aes(x = adult_coef, y = adult_Male_coef, col=adult_female_meth_status))+geom_point()+theme_bw()+geom_smooth()

ggplot(data = metunmeth, aes(x = adult_Male_coef, col=adult_female_meth_status))+theme_bw()+geom_density()
ggplot(data = metunmeth, aes(x = adult_Male_coef*(3.39-adult_coef), col=adult_female_meth_status))+theme_bw()+geom_density()
ggplot(data = metunmeth[grep("[m|f]$", metunmeth$devsexbias),], aes(x = adult_Male_coef, col=adult_female_meth_status))+theme_bw()+geom_density()

ggplot(data = metunmeth[grep("[m|f]$", metunmeth$devsexbias),], aes(x = adult_coef, y = adult_Male_coef, col=adult_female_meth_status))+theme_bw()+geom_density2d()
