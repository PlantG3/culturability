################################################################################
### Gene comparison
### Sanzhen Liu
### 11/14/2020
################################################################################
setwd(".")
### required packages
library(DESeq2)

### log file
log_file <- "1E.1o-DE.log"
cat("DE comparison\n", file = log_file)

### gene information
genefun <- read.delim("path_to_your/B73Ref4.ensembl46.functional.annotation.txt")

### RC data that have removed mRNAs that average RC larger than 2 ###
datapath="../1D_readcounts"
d <- read.delim(paste0(datapath, "/1D.1o-F2callus.raw.read.counts"), check.names=F)

# experimental design:
expdesign <- read.delim("expdesign.txt", stringsAsFactors=F)

### comparison information (subject to change)
group1 <- "II"
group2 <- "I"
comp <- paste0(group1, "_", group2)

### subset
d2 <- d[, c("Gene", as.character(expdesign$Sample))]
d2 <- d2[rowMeans(d2[, -1]) >= 5, ]

### sample design information in the DESeq format
sample.info <- data.frame(row.names=expdesign$Sample,
                          type=as.factor(expdesign$type))

### DESeq analysis
genes <- d2[, 1] ### mRNA sequence
dds <- DESeqDataSetFromMatrix(countData=as.matrix(d2[, -1]),
                              colData=sample.info,
                              formula(~type))

### walt test result
dds <- DESeq(dds)

### test result
res <- results(dds, independentFiltering=T, contrast = c("type", group1, group2))
res <- data.frame(res)
res$Gene <- as.character(genes)
res <- res[, c("Gene", "log2FoldChange", "pvalue", "padj")]
colnames(res) <- c("Gene", paste0(comp, ".", c("log2fc", "pval", "qval")))

out <- merge(d2, res, by = "Gene", all = T)

#### filter
fdr.cutoff <- 0.1
sig <- !is.na(res[, paste0(comp, ".qval")]) & res[, paste0(comp, ".qval")] <= fdr.cutoff
nsig <- sum(sig)
nup <- sum(sig & res[, paste0(comp, ".log2fc")] > 0)
ndown <- sum(sig & res[, paste0(comp, ".log2fc")] < 0)
cat(group1, " vs ", group2, ": ", nsig, "(up=", nup, "; down=", ndown, ")", sep = "", file = log_file, append = T)
cat("\n", file = log_file, append = T)

#out <- merge(normd, res, by = "Gene", all = T)

### output
output_file <- paste0("1E.1o-", group1, "_", group2, ".DE")
write.table(out, output_file, quote=F, row.names=F, sep="\t")


