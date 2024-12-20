fdr_cutoff <- 0.1
genes <- read.delim("path_to_your/B73Ref4.ensembl46.functional.annotation.txt", stringsAsFactors=F)
dt <- read.delim("2E.1o-fast_slow.DE", stringsAsFactors=F)
head(dt)
go <- read.delim("path_to_your/B73Ref4.gene2go.v1.0.txt", stringsAsFactors=F)
source("../sourcecode/goseq.R")
source("../sourcecode/goplot.R")

# go enrichment
dt4go <- dt[!is.na(dt$fast_slow.qval), ]
sig <- dt4go$fast_slow.qval<=fdr_cutoff
nosig <- dt4go$fast_slow.qval>fdr_cutoff
dt4go$Sig <- "no"
dt4go$Sig[sig] <- "yes"
nrow(dt4go)

fc_colname <- "fast_slow.log2fc"
system("mkdir 2E.2o-GO")
goup <- goseq.auto(data=dt4go, godb=go, geneheader="Gene", sigcolname="Sig",
                   nsampling=100000, rawdatacol=c(2:5),
                   log2colname=fc_colname, padjmethod="BH", qvalcutoff=NULL,
                   up.down="up", pvalcutoff=0.01,
                   outpath="2E.2o-GO")

godn <- goseq.auto(data=dt4go, godb=go, geneheader="Gene", sigcolname="Sig",
                   nsampling=100000, rawdatacol=c(2:5),
                   log2colname=fc_colname, padjmethod="BH", qvalcutoff=NULL,
                   up.down="down", pvalcutoff=0.01,
                   outpath="2E.2o-GO")

gosig <- goseq.auto(data=dt4go, godb=go, geneheader="Gene", sigcolname="Sig",
                   nsampling=100000, rawdatacol=c(2:5),
                   log2colname=fc_colname, padjmethod="BH", qvalcutoff=NULL,
                   pvalcutoff=0.01, outpath="2E.2o-GO")

# significant GO in up
pdf("2E.2o-enrichGO.Up.pdf", width=6, height=4)
goup <- read.delim("2E.2o-GO/Sig.up.enriched.goseq.txt", stringsAsFactors=F)
goup4plot <- goup[goup$qvalue<0.05, ]
goup4plot$Term <- gsub(",.*", "", goup4plot$Term)
goplot(go.out=goup4plot, main.space=2, order.by = "pvals",
       main="Up DEGs in fast versus slow",  pval.cex=0.6,
       term.space=13.5, pval.cutoff = 0.05, xlim = c(0, 150), add.total=F)
dev.off()

# significant GO in dn
pdf("2E.2o-enrichGO.Dn.pdf", width=6, height=4.5)
godn <- read.delim("2E.2o-GO/Sig.down.enriched.goseq.txt", stringsAsFactors=F)
godn4plot <- godn[godn$qvalue<0.05, ]
godn4plot$Term <- gsub(",.*", "", godn4plot$Term)
goplot(go.out=godn4plot, main.space=2, order.by = "pvals",
       main="Down DEGs in fast versus slow", pval.cex=0.6,
       term.space=13, pval.cutoff=0.05, xlim=c(0, 120), add.total=F)
dev.off()

# significant GO in both up&down
pdf("2E.2o-enrichGO.Up_Dn.pdf", width=6, height=8)
gosig <- read.delim("2E.2o-GO/Sig.up_down.enriched.goseq.txt", stringsAsFactors=F)
gosig4plot <- gosig[gosig$pvalue<0.05, ]
gosig4plot$Term <- gsub(",.*", "", gosig4plot$Term)
gosig4plot <- gosig4plot[!duplicated(gosig4plot$Term), ]
goplot(go.out=gosig4plot, main.space=2, order.by = "pvals",
       main="DEGs in fast versus slow",  pval.cex=0.6,
       term.space=18, pval.cutoff = 0.05, xlim = c(0,545), add.total=F)
dev.off()

