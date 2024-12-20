#########################segmentation
library("genoseg")
library("DNAcopy")
setwd(".")
source("../sourcecode/genosegF2.R")

### load data
gd <- read.delim("1o-snps.genotype.HiII.txt", stringsAsFactors = F)
head (gd)



genocolumns <- grep("Hi", colnames(gd))

length((genocolumns));length(genocolumns1);length(genocolumns2)
### generate segments ###need create function for RIL Hi-II will use RIL
segres <- genosegF2(geno = gd, genocols = genocolumns, chrname = "CHR", posname = "POS",
                    output.common = "2o-geno.Hi-II", chromosomes = paste0( "",1:10), min.seg.size = 200000,
                    allele1.name = "A", allele2.name = "B", hetero.name = "H", missing.name = 0)



###
seg2score(seg.input = "2o-geno.Hi-II.seg.txt", segscore.output = "3o-geno.segmarker.hi-II.txt", missing.data.code = 0, binsize = 1000000)
