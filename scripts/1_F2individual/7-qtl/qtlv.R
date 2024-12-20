setwd(".")
options(bitmapType = 'cairo')
library(qtl)
######Read the data
qtld <-read.cross(format = "csvr", file = "qtldata.csv",na.strings = "-",genotypes =  c("AA","AB","BB") ,alleles = c("A","B"))

###########Check data
summary(qtld)
plot(qtld)
plotMissing(qtld)


png(filename = "genotypes2.png", pointsize = 12)
par(mfrow = c(1, 2), las = 1)
plot(ntyped(qtld), ylab = "No. typed markers", main = "No. genotypes by individual")
plot(ntyped(qtld, "mar"), ylab = "No. typed individuals", main = "No. genotypes by marker")
dev.off()


markerstodrop<-names(ntyped(qtld, "mar")[ntyped(qtld, "mar")<80])
qtld<-drop.markers(qtld,markerstodrop)
qtld <- calc.errorlod(qtld, error.prob=0.01)
qtld <- calc.genoprob(qtld, step=0, error.prob=0.01)

####method hk Haley-Knott regression (Haley and Knott 1992)
phenoName<-c("CallusType")
for (i in 1:length(phenoName)){
  pn<-phenoName[i]
  m<-"em"
  imqtl <- scanone(qtld,method=m,model="binary",pheno.col = i)
  imqtl
  imqtl.nperm <- scanone(qtld, pheno.col = i,method=m, model="binary", n.perm = 1000)
  #plotMap(qtld)
  ## Doing permutation in batch mode ...
  chrset<- imqtl[!duplicated(imqtl$chr),] 
  chrset<-chrset$chr
  imqtl[order(imqtl$pos),]
  qtl.plot(imqtl, main.text  = paste0("QTL mapping on pheno ",pn,"(method=",m,")"),
           plot.filename = paste0("2o-",pn,"(method=",m,")",".pdf"),
           chr.set = chrset,cexaxis = 1.2,nolines = F,saveplot = T, 
           plot.heigh=3.2,h = summary(imqtl.nperm, alpha=0.05))
  h = as.vector(t(summary(imqtl.nperm, alpha=0.05)))
  hilod<-imqtl[imqtl$lod >= h,]
  x<-genopheno[genopheno$geno.marker %in% row.names(hilod),]
  y<-genopheno[i,]
  df<-rbind(y,x)
  write.table(df,file = paste0("2o-",pn,"(method=",m,")","-hilod.txt"),quote=F,sep="\t")
  
}
