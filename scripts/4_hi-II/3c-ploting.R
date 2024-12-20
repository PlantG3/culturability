setwd(".")

####output code as (1,2,3(hetero),0)
#install.packages("ggplot2")
library("ggplot2")
#install.packages("RColorBrewer")
library(RColorBrewer)
display.brewer.all()


#############################################################################################################
#plot
source("../sourcecode/segplot_HiII.R")


out.prefix="Genotype of chromsome segments of "
width=5
height=3
pdf(paste0(out.prefix, "Hi-II.withoutSNP.pdf"), width=width, height=height,  pointsize = 4)

segplot.ind.HiII.nosnp(segfile="../2o-geno.Hi-II.seg.txt",
                 snpfile="../1o-snps.genotype.HiII.txt",
                 chrsizefile="path_to_your/chr.size.refgen4.txt",
                 chr.names=1:10, out.prefix="Genotype of ", geno.code=c("A","B","H"),geno.col=brewer.pal(n = 8, name = "Dark2")[c(1,2,3)],
                 out.fmt="pdf",width=5,height=5,cex.main=1,cex.axis = 1, cex.lab = 1.2,snp.plot=F,saveplot=F)
dev.off()

