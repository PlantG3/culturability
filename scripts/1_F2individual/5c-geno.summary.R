setwd(".")
source("./genocall.fromGATK.R")
#install.packages("dplyr")#http://www.datasciencemadesimple.com/drop-variables-columns-r-using-dplyr/
library(stats)
library(dplyr)#to drop column

snps <- read.delim("../2o-MT2B73Ref4.select.txt", stringsAsFactors = F)

nrow(snps)
colnames(snps)
head(snps); nrow(snps)
#
### NA to 0
for (i in 6:ncol(snps)) {
  snps[is.na(snps[, i]), i] <- 0
}

ncol(snps)
snps1<-snps[rowSums(snps[,18:261])>0,]
nrow(snps1);ncol(snps1)
snps1
snps<-snps1


### classification
nsnps <- nrow(snps)
nsnps

### geno all
calli <- colnames(snps)[grep("XT.*REF", colnames(snps))]
calli.samples <- gsub("_REF", "", calli)
for (ecs in calli.samples) {
  cat(ecs, "\n")
  snps[, ecs] <- apply(snps[, paste0(ecs, c("_REF", "_ALT"))], 1, genocall,
                       homo=4, homo.perc=0.9, heter=2, heter.perc=0.2,
                       overall.perc=0.9, ref.code="B", alt.code="A",
                       heter.code="H", missing.code="0", nocall.code = "0")
}

###
### output
###
snps.sub <- snps[, c(1:4, grep("XT.*[0-9]$", colnames(snps)))]
head(snps.sub)
nrow(snps.sub)
ncol(snps.sub)
colnames(snps.sub)
#snps.sub.1<-snps.sub[!(rowSums(snps.sub[,5:123])==0),]

write.table(snps.sub, "1o-snps.genotype.txt", row.names = F, quote = F, sep = "\t")
