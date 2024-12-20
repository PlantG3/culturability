setwd(".")
source("./genocall.fromGATK.R")
#install.packages("dplyr")#http://www.datasciencemadesimple.com/drop-variables-columns-r-using-dplyr/
library(stats)
library(dplyr)#to drop column


### loading data
snps <- read.delim("../2o-MT2B73Ref4.select.txt", stringsAsFactors = F)
head(snps); nrow(snps)

### NA to 0
for (i in 6:ncol(snps)) {
  snps[is.na(snps[, i]), i] <- 0
}

snps<-snps[rowSums(snps[,18:257])>0,]
nrow(snps1);ncol(snps)
#snps<-snps1
# A188-1
snps$A188R1 <- apply(snps[, c("A188.1_REF", "A188.1_ALT")], 1, genocall,
                     homo=8, homo.perc=0.95, heter=2, heter.perc=0.2,
                     overall.perc=0.9, ref.code="B", alt.code="A",
                     heter.code="H", missing.code="-", nocall.code = "-")

## A188-2
snps$A188R2 <- apply(snps[, c("A188.2_REF", "A188.2_ALT")], 1, genocall,
                     homo=8, homo.perc=0.95, heter=2, heter.perc=0.2,
                     overall.perc=0.9, ref.code="B", alt.code="A",
                     heter.code="H", missing.code="-", nocall.code = "-")


# B73R1
snps$B73R1 <- apply(snps[, c("B73R1_REF", "B73R1_ALT")], 1, genocall,
                    homo=8, homo.perc=0.95, heter=2, heter.perc=0.2,
                    overall.perc=0.9, ref.code="B", alt.code="A",
                    heter.code="H", missing.code="-", nocall.code = "-")

## B73R2
snps$B73R2 <- apply(snps[, c("B73R2_REF", "B73R2_ALT")], 1, genocall,
                    homo=8, homo.perc=0.95, heter=2, heter.perc=0.2,
                    overall.perc=0.9, ref.code="B", alt.code="A",
                    heter.code="H", missing.code="-", nocall.code = "-")

### summary
table(snps$A188R1)
table(snps$A188R2)
table(snps$B73R1)
table(snps$B73R2)

### filter
snps <- snps[snps$A188R1 == "A" & snps$A188R2 == "A" &
               snps$B73R1 == "B" & snps$B73R2 == "B", ]

### classification
nsnps <- nrow(snps)
nsnps

### geno all
calli <- colnames(snps)[grep("HiII[AB]_REF", colnames(snps))]
calli.samples <- gsub("_REF", "", calli)
for (ecs in calli.samples) {
  cat(ecs, "\n")
  snps[, ecs] <- apply(snps[, paste0(ecs, c("_REF", "_ALT"))], 1, genocall,
                       homo=4, homo.perc=0.9, heter=1, heter.perc=0.2,
                       overall.perc=0.9, ref.code="B", alt.code="A",
                       heter.code="H", missing.code="0", nocall.code = "0")
}

###
### output
###
snps.sub <- snps[, c(1:4, grep("HiII[AB]$", colnames(snps)))]
head(snps.sub)
nrow(snps.sub)
snps.sub<-snps.sub[snps.sub$HiIIA!=0 | snps.sub$HiIIB!=0,]
write.table(snps.sub, "1o-snps.genotype.HiII.txt", row.names = F, quote = F, sep = "\t")
snps.sub<-read.delim("1o-snps.genotype.HiII.txt",stringsAsFactors = F)
snps.sub<-snps.sub[snps.sub$HiIIA!=0 | snps.sub$HiIIB!=0,]
nrow(snps.sub)

###################################################################################################
###filter
##calculate the missing rate per marker - mm 
##calculate the missing rate per individual- mi -not finished yet
par(mfrow=c(1,3))
dfm<-snps.sub[,5:(ncol(snps.sub))]
colnames(dfm)<-NULL
rownames(dfm)<-NULL

mm<-c()
y<-ncol(dfm)
for (n in 1:nrow(dfm)){
  dftest<-dfm[n,]
  dftest<-as.vector(t(dftest))
  x<-sum(as.integer(dftest=="0"))
  
  Mm<-x/y
  mm<-c(mm,Mm)
}
snps[15,]

hist(mm,main="Histogram of Missing rate per markers")



dfm<-snps.sub[,5:(ncol(snps.sub))]
colnames(dfm)<-NULL
rownames(dfm)<-NULL

mi<-c()
y<-nrow(dfm)
for (n in 1:ncol(dfm)){
  dftest<-dfm[,n]
  x<-sum(as.integer(dftest=="0"))
  
  Mi<-x/y
  mi<-c(mi,Mi)
}

hist(mi,main="Histogram of Missing rate per individual")






##calculate the heterozygote ratio
dfm<-snps.sub[,5:(ncol(snps.sub))]
colnames(dfm)<-NULL
rownames(dfm)<-NULL
dfmh<-dfm
dfma<-dfm
h<-c()
for (n in 1:ncol(dfm)){
  dfmh[dfm[,n]=="H",n]<-1
  dfmh[dfm[,n]!="1",n]<-0
  dfma[dfm[,n]!="0",n]<-1
  dfma[dfm[,n]=="0",n]<-0
  h<-x/y
}
x<-rowSums(dfmh)
y<-rowSums(dfma)

length((h))
hist(h,main="Histogram of heterozygous genotype")


dev.off()


pdf("genotype-distribution-before-filtering.pdf",width=12,height=4)
par(mfrow=c(1,3))
hist(mm,main="Histogram of Missing rate per markers")
hist(mi,main="Histogram of Missing rate per individual")
hist(h,main="Histogram of heterozygous genotype")


dev.off()


#######################################################################################  
dfmn<-snps.sub
dfmnS<-dfmn[h==0,]
nrow( dfmnS)
dfmn<-dfmn[mm<0.98,]
nrow(dfmn)
ncol(dfmn)
nrow(snpsecs.all)
dfmnp<-merge(snpsParents,dfmn, by=c("CHR","POS","REF","ALT"),all.y=T)
nrow(dfmnp)

write.table(dfmn, "1o-snps.genotype.filter.HiII.txt", row.names = F, quote = F, sep = "\t")

















