setwd(".")
source("../sourcecode/genosegF2.R")

### load data
gd <- read.delim("1o-snps.genotype.txt", stringsAsFactors = F)
head (gd)

library("genomap")
library("DNAcopy")


genocolumns <- grep("XT.", colnames(gd))
genocolumns1 <- grep("XT.II", colnames(gd))
genocolumns2 <- grep("XT.I_", colnames(gd))

length((genocolumns));length(genocolumns1);length(genocolumns2)
### generate segments
segres <- genosegF2(geno = gd, genocols = genocolumns, chrname = "CHR", posname = "POS",
                    output.common = "2o-geno", chromosomes = paste0( "",1:10), min.seg.size = 1000000,
                    allele1.name = "B", allele2.name = "A", hetero.name = "H", missing.name = 0,
                    seg.mean.cutoffs = c(-0.75, -0.3, 0.3, 0.75))


###
seg2score(seg.input = "2o-geno.seg.txt", segscore.output = "3o-geno.segmarker.txt", missing.data.code = 0, binsize = 1000000)



