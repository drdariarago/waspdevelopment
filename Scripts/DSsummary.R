## stats for poster

# tabulate stage-specific and sex-specific genes
stspec<-apply(fit2_splicing$fdr_splicing, 1, function(x){any(x[1:5]<0.00001)})
sxspec<-apply(fit2_splicing$fdr_splicing, 1, function(x){any(x[6:10]<0.00001)})
addmargins(table(stspec, sxspec),margin=c(1,2))
prop.table(table(stspec, sxspec))
prop.table(table(stspec, sxspec), margin=1)
# is there any interaction betwen the two?
fisher.test(table(stspec, sxspec))


prop.table(table(stspec, sxspec), margin=2)


addmargins(table(stspec, sxspec), margin=1)

# How many genes are DS in any number of stages?
table(apply(fit2_splicing$fdr_splicing[,1:5],1,function(x) sum(x<0.00001)))

# How many genes are sex DE in any number of stages?
table(apply(fit2_splicing$fdr_splicing[,6:10],1,function(x) sum(x<0.00001)))

# tabulate single-stage specific expression for sex
stspec1<-apply(fit2_splicing$fdr_splicing, 1, function(x){sum(x[1:5]<0.00001)==1})
sxspec1<-apply(fit2_splicing$fdr_splicing, 1, function(x){sum(x[6:10]<0.00001)==1})
addmargins(table(stspec1, sxspec1),margin=c(1,2))
prop.table(table(stspec1, sxspec1))
prop.table(table(stspec1, sxspec1), margin=1)

# How many stage specific genes are also single stage specific?
addmargins(table(stspec, stspec1),margin=c(2))
prop.table(table(stspec, stspec1),margin=c(2))
# How many sex specific genes are also single-stage specific?
addmargins(table(sxspec, sxspec1),margin=c(2))
prop.table(table(sxspec, sxspec1),margin=c(2))