# function 5: all chromosome in the same figure:
pval.oneplot <- function (input, yrange = NULL, xlab.text="chromosome",
                      mar.set =c(4.5, 4.5, 0.5, 0.5), mgp.set = c(2.5, 0.8, 0),
                      cols = c("grey16", "grey32"), pvalue.cutoff = NULL,
                      ylab.text = "-log10(p-values)", na.rm = T, pval_nonzero = T,
                      main.text="", axisline.width=1.5, chr.set, plot.cex.lab = 1.3,
                      chr.size, order.by.chrsize=F, label.rm=NULL, cexaxis=1.3,
                      cex.line.width = 0.5, saveplot=FALSE, plot.width=3, plot.heigh=2.5, 
                      plot.path=".", plot.filename="default.png", ...) {
  
  # input should have three columns, "Chr", "Pos", "Pval"
  # size contains length information for every chromosome or contigs
  # size has two columns: chr and size
  plot.path <- gsub("/$", "", plot.path)
  colnames(input) <- c("Chr", "Pos", "Pval")
  input$Chr <- as.character(input$Chr)
  input <- input[input$Chr %in% chr.set, ]
  
  if (na.rm) {
    input <- input[!is.na(input$Pval), ]
  }
  
  if (pval_nonzero) {
    input[input$Pval == 0, "Pval"] = min(input[input$Pval != 0, "Pval"])
  }
  
  input$neglogP <- -log10(input$Pval)
  # chromosome length:
  colnames(chr.size) <- c("Chr", "Size")
  chr.size$Chr <- as.character(chr.size$Chr)
  chr.size <- chr.size[chr.size$Chr %in% chr.set, ]
  
  if (order.by.chrsize) {
    chr.size <- chr.size[order(chr.size$Size, decreasing=T), ]
  }
  
  # to judge the number is odd or even:
  odd <- function(x) {
    if ((round(x/2,0) - i/2)==0) {
      y <- 0
    } else {
      y <- 1
    }
    return(y)
  }
  
  # plotting
  accum <- 0
  all.col <- NULL
  all.chr <- NULL
  centers <- NULL
  gap <- sum(chr.size$Size)/100
  if (is.null(yrange)) {
    ymax <- max(input$neglogP)
    yrange <- c(0, ymax)
  } else {
    ymax <- max(yrange)
  }
  
  if (saveplot) {
    #pdf(paste(plot.path, plot.filename, sep="/"), width=plot.width, heigh=plot.heigh)
    png(paste(plot.path, plot.filename, sep="/"), pointsize = 4, res = 2000, units = "in", 
        width=plot.width, heigh=plot.heigh)
  }
  
  ###
  ### plot
  ###
  
  par(mar = mar.set, mgp = mgp.set)
  plot(NULL, NULL, ylim = yrange,
       xlim = c(0, gap * nrow(chr.size) + sum(chr.size$Size)),
       xaxt = "n", yaxt = "n",
       xlab = xlab.text, ylab = ylab.text,
       main = main.text, cex.lab = plot.cex.lab,
       bty = "n", cex.axis = cexaxis, lwd = cex.line.width)
  
  axis(side = 2, cex.axis = cexaxis, lwd = cex.line.width)
  box(lwd = cex.line.width)
  
  
  all.accum <- NULL
  
  # cols of points
  if (length(cols) < nrow(chr.size)) {
    cols <- rep(cols, ceiling(nrow(chr.size) / length(cols)))
  }

  # plot points for each chromosome
  for (i in 1:(nrow(chr.size))) {
    all.accum <- c(all.accum, accum)
    pre.accum <- accum
    chr <- chr.size[i, "Chr"]
    len <- chr.size[i, "Size"]
    plot.col <- cols[i]
    
    pos <- input[input$Chr==chr, "Pos"]
    prob <- input[input$Chr==chr, "neglogP"]
    #lines(c(accum, accum+len), c(-(ymax/50), -(ymax/50)), col=plot.col, lwd=axisline.width, lend=1)
    points(accum+pos, prob, col=plot.col, lwd=cex.line.width, ...)
    accum <- accum + len + gap
    center.point <- (pre.accum + accum - gap)/2
    all.col <- c(all.col, plot.col)
    all.chr <- c(all.chr, chr)
    centers <- c(centers, center.point)
  }
  if (!is.null(label.rm)) {
    for (each.label.rm in label.rm) {
      all.chr <- gsub(each.label.rm, "", all.chr)
    }
  }
  print(all.chr)
  axis(side = 1, at = centers, labels=all.chr, tick=F, cex.axis=cexaxis, lwd = cex.line.width)
  
  # pvalue cutoff
  if (!is.null(pvalue.cutoff)) {
	abline(h = -log10(pvalue.cutoff), col = "red", lwd = cex.line.width, lty = 2)#darkorange2
  }

  if (saveplot) dev.off()
  
  names(all.accum) <- chr.set
  return(all.accum)
}
