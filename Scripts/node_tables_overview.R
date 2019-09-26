## Check proportions of nodes with sexbias in early vs late stage
transcriptdata <- read.csv(file = "./Finaltables/transcriptdata_full.csv")[,-1]
nodedata <- read.csv(file = "./Finaltables/nodedata.full.csv")[,-1]

## For generic nodes
# emb10 only
prop.table(table(
  grepl(pattern = "^[m,f]", x = nodedata$devsexbias), 
  grepl(pattern = "[m,f]$", x = nodedata$devsexbias),
  dnn = c("emb10bias","adultbias")
), margin = 1)
fisher.test(
  table(
    grepl(pattern = "^[m,f]", x = nodedata$devsexbias), 
    grepl(pattern = "[m,f]$", x = nodedata$devsexbias)
  ))
# all preadult
prop.table(table(
  grepl(pattern = "[m,f].+", x = nodedata$devsexbias), 
  grepl(pattern = "[m,f]$", x = nodedata$devsexbias),
  dnn = c("preadultbias","adultbias")
), margin = 1)
fisher.test(
  table(
    grepl(pattern = "[m,f].+", x = nodedata$devsexbias), 
    grepl(pattern = "[m,f]$", x = nodedata$devsexbias)
  ))



## Check if there is convergence in sex bias direction in early vs late
biases = data.frame( 
  preadultbias = ifelse(grepl(pattern = "m.+", transcriptdata$devsexbias), "M", ifelse(grepl(pattern = "f.+", x = transcriptdata$devsexbias), "F", "Unbiased")),
  adultbias = ifelse(grepl(pattern = "m$", transcriptdata$devsexbias), "M", ifelse(grepl(pattern = "f$", x = transcriptdata$devsexbias), "F", "Unbiased")),
  emb10bias = ifelse(grepl(pattern = "^m", transcriptdata$devsexbias), "M", ifelse(grepl(pattern = "^f", x = transcriptdata$devsexbias), "F", "Unbiased"))
)

table(biases$preadultbias, biases$adultbias, dnn = c("preadult","adult"))
prop.table(table(biases$preadultbias, biases$adultbias, dnn = c("preadult","adult")), margin = 2)

fisher.test(table(biases$preadultbias, biases$adultbias, dnn = c("preadult","adult"))[c(1,2),c(1,2)])
prop.table(table(biases$preadultbias, biases$adultbias, dnn = c("preadult","adult"))[c(1,2),c(1,2)], margin = 1)

table(biases$emb10bias, biases$adultbias, dnn = c("emb10bias","adult"))
prop.table(table(biases$emb10bias, biases$adultbias, dnn = c("emb10bias","adult"))[c(1,2),c(1,2)], margin = 1)
fisher.test(table(biases$emb10bias, biases$adultbias, dnn = c("emb10bias","adult"))[c(1,2),c(1,2)])


## Cp emb10 vs rest for ortholog/paralog and transcript/splicing
orthotransbias <- data.frame(
  ortho = relevel(nodedata$quality7, ref = "Ortholog"),
  trans = nodedata$event_type,
  bias = grepl(pattern = "^[m,f]", x = nodedata$devsexbias),
  prepupal = grepl(pattern = "[m,f].{2,}", x = nodedata$devsexbias)
)
orthotransbias

library(MuMIn)

glm1 <- glm(formula = bias ~ ortho*trans, family = binomial, data = orthotransbias)
dredge(glm1)
plot(glm1)
summary(glm1)

glm2 <- glm(formula = prepupal ~ ortho*trans, family = binomial, data = orthotransbias)
dredge(glm2)
plot(glm2)
summary(glm2)

prop.table(table(orthotransbias$bias, orthotransbias$ortho), margin = 2)
prop.table(table(orthotransbias$bias, orthotransbias$trans), margin = 2)

prop.table(table(orthotransbias$prepupal, orthotransbias$ortho), margin = 2)
prop.table(table(orthotransbias$prepupal, orthotransbias$trans), margin = 2)
