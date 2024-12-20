setwd(".")
### threshold
res <- read.delim("mapping.result.txt", stringsAsFactors = F)

pval.cutoff <- 0.05 / sum(!is.na(res$pval))
cat("bonferroni correction pvalue cutoff", pval.cutoff, "\n")
-log10(pval.cutoff)
ymax <- 15
pval.cutoff <-7.337295e-07  ### bonferroni correction pvalue cutoff 7.337295e-07 

source("../pvalbsr.oneplot.R")
chrsize <- read.delim("path_to_your/chr.size.refgen4.txt")


pdf("3o-mapping.BSA-4c-mapping.pdf", width=6, height=4)
pvalbsr.oneplot(input1=res[res$log2odd < 0, c(1, 2, 14)],input2 =res[res$log2odd > 0, c(1, 2, 14)],  chr.set=1:10,
                chr.size=chrsize,cols2=c("olivedrab4", "lightsalmon3"),
                pvalue.cutoff=pval.cutoff, pch=19, cex=0.5, cex.line.width=0.5, ylab.text="-log10(p value)",
                main.text="BSR mapping of callus types", mar.set=c(4,4,4,1), cexaxis=1.5, plot.cex.lab=1.5)
dev.off()
