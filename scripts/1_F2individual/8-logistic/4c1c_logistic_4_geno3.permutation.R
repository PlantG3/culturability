setwd(".")

########################################################################################
# data
########################################################################################
data <- read.delim("3o-geno.segmarker.txt", stringsAsFactors=F)

########################################################################################
# function to perform logistic regression
########################################################################################
logist <- function(d, marker_col=1, geno_cols, missing_code="0") {
  ### logistic regression to identify the impact of genotypes on callus types (I and II)
  marker <- d[marker_col]
  geno <- as.character(d[geno_cols])
  pheno <- as.character(names(d)[geno_cols])
  pheno <- gsub("XT.", "", gsub("_.*", "", pheno))
  geno <- as.character(geno)
  valid <- (geno != missing_code)
  geno <- geno[valid]
  pheno <- pheno[valid]
  df <- data.frame(geno, pheno)
  logit.full <- glm(factor(pheno)~factor(geno), data=df, family=binomial(link="logit"))
  logit.reduced <- glm(factor(pheno)~1, data=df, family=binomial(link="logit"))
  aov <- anova(logit.full, logit.reduced, test="LRT")
  pval <- aov[2,5]
  pval
}

########################################################################################
# perform logistic regression
########################################################################################
indata <- data[, -c(1,3,4,5)]
geno_cols <- grep("^XT", colnames(indata))
pvals <- apply(indata, 1, logist, geno_cols=geno_cols)

##################################################################
perm_minPval <- NULL
marker.name <- colnames(indata)[1]
original.names <- colnames(indata)[-1]
for (i in 1:1000) {
  new.names <- sample(original.names, length(original.names))
  permindata <- indata
  colnames(permindata) <- c(marker.name, new.names)
  permpvals <- apply(permindata, 1, logist, geno_cols=geno_cols)
  perm_minPval <- c(perm_minPval, min(permpvals, na.rm=T))
  if (i/100 == round(i/100)) {
    cat(i, ":", min(permpvals, na.rm=T), "\n")
  }
}
#save(perm_minPval,file="perm_minPval.RData")
load("perm_minPval.RData")
perm_pval_cutoff <- quantile(c(perm_minPval, min(pvals)), 0.05)
perm_pval_cutoff <- quantile(c(perm_minPval, min(pvals)), 0.10)
perm_pval_cutoff <- quantile(c(perm_minPval, min(pvals)), 0.15)
cat("permutation pcutoff=", perm_pval_cutoff, "\n")
##################################################################
#perm_pval_cutoff <- 6.460798e-05  # from 1000 permutations

perm_pval_cutoff<-0.0002222313 
########################################################################################
# pvalue cutoffs (FDR)
########################################################################################
qvals <- p.adjust(pvals, method="BH")

qval_cutoff_strict <- 0.01
sum(qvals<qval_cutoff_strict)
pval_cutoff_strict <- max(pvals[qvals<qval_cutoff_strict])
-log10(pval_cutoff_strict)

qval_cutoff_loose <- 0.05
sum(qvals<qval_cutoff_loose)
pval_cutoff_loose <- max(pvals[qvals<qval_cutoff_loose])
-log10(pval_cutoff_loose)

pval_cutoff <- c(pval_cutoff_loose, pval_cutoff_strict)

out <- cbind(data[, c("Chr", "Pos")], pvals)
out2 <- cbind(data[, c("Chr", "Pos", "Marker")], pvals)
write.table(out2, "4c1o_logistic.test.txt", quote=F, row.names=F, sep="\t")

########################################################################################
# plotting
########################################################################################
source("../sourcecode/pval.oneplot.R")
chrsize <- read.delim("sourcecode/chr.size.refgen4.txt")
out<-read.delim("./4c1o_logistic.test.txt",stringsAsFactors = F)
out<-out[,-3]


perm_pval_cutoff <- quantile(c(perm_minPval, min(pvals)), 0.05)
perm_pval_cutoff10 <- quantile(c(perm_minPval, min(pvals)), 0.10)
perm_pval_cutoff15 <- quantile(c(perm_minPval, min(pvals)), 0.15)
pdf("4c1o_logistic.chrplot-permuqvalue02092021.newname.pdf", width=6, height=4)
pval.oneplot(input=out, chr.set=1:10, chr.size=chrsize, cols=c("olivedrab4", "lightsalmon3"),
             pvalue.cutoff=perm_pval_cutoff, pch=19, cex=0.5, cex.line.width=1, ylab.text="-log10(p value)",xlab.text = "chromosome",
             main.text="Callus types QTL mapping (logistic regression)", mar.set=c(4,4,4,1), cexaxis=1.2, plot.cex.lab=1.5)
abline(h = -log10(pval_cutoff_strict), col = "dark grey", lwd = 1, lty = 2)
#abline(h = -log10(perm_pval_cutoff10), col = "dark grey", lwd = 1, lty = 2)

dev.off()