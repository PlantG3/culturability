setwd(".")
source("../sourcecode/binomial.R")
infile <- "4o-snp.variants.select.selected.txt"
d <- read.delim(infile)
datacols <- 6:13
minSum <- 50
stattest.data <- d[rowSums(d[, datacols]) >= minSum, ]
dx2 <- apply(stattest.data[, 6:13], 1, binomial.od, cat1pos = c(1, 3, 5, 7),
			cat2pos = c(2, 4, 6, 8), geno = c("I", "I", "II", "II"),
			lambda = 3.138632)

dx2t <- t(dx2)
colnames(dx2t) <- c("I1_REF", "I1_ALT", "I2_REF", "I2_ALT", "II1_REF", "II1_ALT", "II2_REF", "II2_ALT", "log2odd", "pval")
out <- cbind(stattest.data[, 1:4], dx2t)
write.table(out,  "3o-mapping.binomial.result.txt", quote = F, row.names = F, sep = "\t")

