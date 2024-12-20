
##############################################################################################################################
#function:plot Hi-II AB in a single plot without SNP
segplot.ind.HiII.nosnp<-function(segfile="4b_geno3/4b2o_F2.geno.seg.filt.txt",
                           snpfile="4b_geno3/4b1o_F2.geno.txt",
                           chrsizefile="path_to_your/chr.size.refgen4.txt",
                           chr.names=1:10, out.prefix="1o.seg.plot.B73Ref4", geno.code=c(1,2,3),geno.col=c("dark green","blue","red"),
                           out.fmt="pdf",width=5,height=5,cex.main=3,cex.axis = 2, cex.lab = 2,snp.plot=F,saveplot=T){
  
  #install.packages("RColorBrewer")
  library(RColorBrewer)
  
  
  #read files
  geno.seg<-read.delim(segfile,stringsAsFactors = F)
  snpgeno<-read.delim(snpfile,stringsAsFactors = F)
  #snpgeno$CHR<-snpgeno$CHROM #bug potential
  
  #####label the qtls 
  
  segm<-read.delim("path_to_your/4b2o_F2.geno.seg.filt.score.txt")
  segm$Marker<-paste0(segm$Chr,"_",segm$Pos)
  ajpeak<-read.delim("path_to_your/peaks.adjust.txt")
  sig.m<-ajpeak$GBS_seg_marker
  
  
  ##plot function
  seg.lines <- function(segdata, genocode, segcol, line.cex = 1) {
    eseg.geno <- segdata[segdata$Genotype == genocode, ]
    if (nrow(eseg.geno) > 0) {
      for (k in 1:nrow(eseg.geno)) {
        lines(x = eseg.geno[k, 3:4]/1000000, y = rep((shift), 2), lend = 2, cex = line.cex, col = segcol,lwd=line.cex) 
      }
    }
  }
  
  
  ###plot by individuals
  chrsize <- read.delim(chrsizefile,stringsAsFactors = F, header=F)
  colnames(chrsize)<-c("Chr","Size")
  chrsize<-chrsize[chrsize$Chr %in% chr.names,]
  a<-geno.seg$Individual
  a<-a[!duplicated(a)]
  length(a)
  allsamples<-a
  par(mfrow=c(1,2))
  if (saveplot==T){
    if (out.fmt=="pdf"){
      pdf(paste0(out.prefix, "Hi-II.pdf"), width=width, height=height,  pointsize = 4)
    }
    if (out.fmt=="png"){
      png(paste0(out.prefix, "Hi-II.png"),width=width, height=height,   pointsize = 4)
    }
  }
  for (esample in allsamples) {
    
    ####need to be edited
    # par(mar = c(5, 4, 3, 2))
    
    plot(NULL,NULL, ylim=c(0, 10), xlim=c(0, as.numeric(max(chrsize$Size))/1000000),
         bty="n", type="n",yaxt="n", cex.main=cex.main, main=paste0(out.prefix,esample), 
         xlab="Physical coordinate (Mb)", ylab="", cex.axis = cex.axis, cex.lab = cex.axis)
 #   legend(250,10.35,legend = geno.code,col = geno.col,lty=1, cex=cex.lab*(1-0.5))
 legend(200,10.35,legend = c("homozygous A188","homozygous B73","heterozygous"),col = geno.col,lty=1, cex=cex.lab*(1-0.5))
    
    count <- 0 
    
    for (j in order(as.numeric(chrsize$Chr),decreasing = T)) { # 10 chromosomes
      eseg<- geno.seg[geno.seg$Chr ==paste0("", j,"") & geno.seg$Individual == esample, ]
      
      shift <- 1 * count
      
      
      ### seg draw:
      
      seg.lines(segdata  = eseg, genocode = geno.code[1], segcol = geno.col[1], line.cex = 0.5*5)
      seg.lines(segdata  = eseg, genocode = geno.code[2], segcol = geno.col[2], line.cex = 0.5*5)
      seg.lines(segdata  = eseg, genocode = geno.code[3], segcol = geno.col[3],line.cex = 0.5*5)
      if(j %in% ajpeak$Chromosome){
        rect(xleft=ajpeak[ajpeak$Chromosome==j,"pm.start"]/1000000,ybottom = shift-0.4,ytop =shift+0.4,xright =ajpeak[ajpeak$Chromosome==j,"pm.end"]/1000000, lty = 2,border = brewer.pal(n = 8, name = "Dark2")[c(2)],col = NA  )
        points(x = ajpeak[ajpeak$Chromosome==j,"peak.pos"]/1000000,y=shift-0.3, col=brewer.pal(n = 8, name = "Dark2")[c(2)],pch=17,cex=1.5)
      }
      
      if (snp.plot==TRUE){
        dp<-snpgeno[snpgeno$CHR==paste0("", j,""),c("CHR", "POS", esample)]
        dp1 <- dp[dp[, 3] == geno.code[1], ]
        dp2<- dp[dp[, 3] == geno.code[2], ]
        dp3 <- dp[dp[, 3] == geno.code[3], ]
        points(dp1$POS/1000000, rep(1.5/10 + shift, nrow(dp1)), pch = "|", cex = 0.5, col = geno.col[1])
        points(dp2$POS/1000000, rep(3/10 + shift, nrow(dp2)), pch = "|", cex = 0.5, col =geno.col[2])
        points(dp3$POS/1000000, rep(4.5/10 + shift, nrow(dp3)), pch = "|", cex = 0.5, col = geno.col[3])
      } 
      
      ### label
      text(x = -22, y = shift, labels = paste0("Chr ",j,""), cex = cex.lab*(1-0.2), xpd = T)
      count <- count + 1
    }
    if(saveplot==T){
      dev.off()}
  }
  
}


