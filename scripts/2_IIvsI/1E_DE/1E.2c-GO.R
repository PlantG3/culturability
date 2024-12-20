fdr_cutoff <- 0.1
genes <- read.delim("path_to_your/B73Ref4.ensembl46.functional.annotation.txt", stringsAsFactors=F)
dt <- read.delim("1E.1o-II_I.DE", stringsAsFactors=F)
head(dt)
go <- read.delim("path_to_your/B73Ref4.gene2go.v1.0.txt", stringsAsFactors=F)
source("../sourcecode/goseq.R")
source("../sourcecode/goplot.R")

# go enrichment
dt4go <- dt[!is.na(dt$II_I.qval), ]
sig <- dt4go$II_I.qval<=fdr_cutoff
nosig <- dt4go$II_I.qval>fdr_cutoff
dt4go$Sig <- "no"
dt4go$Sig[sig] <- "yes"
nrow(dt4go)

fc_colname <- "II_I.log2fc"
system("mkdir 1E.2o-GO")
goup <- goseq.auto(data=dt4go, godb=go, geneheader="Gene", sigcolname="Sig",
                   nsampling=100000, rawdatacol=c(2:5),
                   log2colname=fc_colname, padjmethod="BH", qvalcutoff=NULL,
                   up.down="up", pvalcutoff=0.01,
                   outpath="1E.2o-GO")

godn <- goseq.auto(data=dt4go, godb=go, geneheader="Gene", sigcolname="Sig",
                   nsampling=100000, rawdatacol=c(2:5),
                   log2colname=fc_colname, padjmethod="BH", qvalcutoff=NULL,
                   up.down="down", pvalcutoff=0.01,
                   outpath="1E.2o-GO")

gosig <- goseq.auto(data=dt4go, godb=go, geneheader="Gene", sigcolname="Sig",
                   nsampling=100000, rawdatacol=c(2:5),
                   log2colname=fc_colname, padjmethod="BH", qvalcutoff=NULL,
                   pvalcutoff=0.01, outpath="1E.2o-GO")

# significant GO in up
pdf("1D.2o-enrichGO.Up.pdf", width=4.5, height=3)
goup4plot <- goup
goup4plot$Term <- gsub(",.*", "", goup4plot$Term)
goplot(go.out=goup4plot, main.space = 0.5, order.by = "pvals", 
       term.space = 11, pval.cutoff = 0.05, xlim = c(0, 25), add.total=F)
dev.off()

# significant GO in dn
pdf("1D.2o-enrichGO.Dn.pdf", width=4.5, height=3)
godn4plot <- godn
godn4plot$Term <- gsub(",.*", "", godn4plot$Term)
goplot(go.out=godn4plot, main.space = 0.5, order.by = "pvals", 
       term.space = 11, pval.cutoff = 0.05, xlim = c(0, 20), add.total=F)
dev.off()

# significant GO in both up&down
pdf("1D.2o-enrichGO.Up_Dn.pdf", width=4.5, height=4)
gosig4plot <- gosig
gosig4plot$Term <- gsub(",.*", "", gosig4plot$Term)
goplot(go.out=gosig4plot, main.space = 0.5, order.by = "pvals", 
       term.space = 11, pval.cutoff = 0.05, xlim = c(0, 40), add.total=F)
dev.off()

