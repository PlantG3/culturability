### setup the working directory:
source("../sourcecode/starRCmerger.R")
source("../sourcecode/starALNsummary.R")
source("../sourcecode/size.factor.DESeq.R")

### read counts
datapath <- "../1C_star"
rc <- starRCmerger(datapath = datapath, suffix = "ReadsPerGene.out.tab")

### normalization factors:
sizefactor <- size.factor.DESeq(rc[, -1])
head(sizefactor)

### normalization
norm <- rc
for (i in 2:ncol(rc)) {
  norm[, i] <- round(rc[, i] / sizefactor[i - 1], 0)
}

### output
write.table(rc, "1D.1o-F2callus.raw.read.counts", row.names = F, sep = "\t", quote = F)
write.table(norm, "1D.1o-F2callus.norm.read.counts", row.names = F, sep = "\t", quote = F)
