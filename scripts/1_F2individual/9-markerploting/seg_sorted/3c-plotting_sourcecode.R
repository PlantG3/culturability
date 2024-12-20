library(RColorBrewer)
segplot.chr<-function(segfile="4b_geno3/4b2o_F2.geno.seg.filt.txt",
                      chrsizefile="path_to_your/chr.size.refgen4.txt",
                      chr.names=1:10, out.prefix="1o.seg.plot.B73Ref4", geno.code=c(1,2,3),
                      geno.col=c("dark green","blue","red"),main.lab="Genotype segregation of chrosome ",
                      out.fmt="pdf",width=5,height=3,cex.main=3,cex.axis = 2, cex.lab = 2){
  
  #read files
geno.seg<-read.delim(segfile,stringsAsFactors = F)
geno.seg1<-geno.seg[,c(2:4,10:12)]
geno.seg2<-geno.seg[,c(2:4,18:20)]

  #geno.seg$Phenotype<-geno.seg$Individual
  
  
  
  ##plot function
  seg.lines <- function(segdata,  segcol) {
    eseg.geno <- segdata
    if (nrow(eseg.geno) > 0) {
      for (k in 1:nrow(eseg.geno)) {
        lines(x = eseg.geno[k, 3:4]/1000000,
              y = rep((shift+1), 2), lend = 2, 
              cex = line.cex, col = segcol,lwd=line.cex) 
      }
    }
  }
  
  
  ##plot function
  seg.rects <- function(segdata,  segcol) {
    eseg.geno <- segdata
      for (k in 1:nrow(eseg.geno)) {
       
        rect(xleft =eseg.geno[k, 2]/1000000, ybottom = 0 ,
             ytop = eseg.geno[k, 6], xright = eseg.geno[k, 3]/1000000, 
             col=segcol[3],border =  segcol[3]) 
        
        rect(xleft =eseg.geno[k, 2]/1000000, ybottom = eseg.geno[k, 6] ,
             ytop = eseg.geno[k, 6]+eseg.geno[k, 5], xright = eseg.geno[k, 3]/1000000, 
             col=segcol[2],border =  segcol[2]) 
        
        rect(xleft =eseg.geno[k, 2]/1000000, ybottom = eseg.geno[k, 6]+eseg.geno[k, 5] ,
             ytop = eseg.geno[k, 6]+eseg.geno[k, 5]+eseg.geno[k, 4], 
             xright = eseg.geno[k, 3]/1000000, 
             col=segcol[1],border =  segcol[1]) 
      
    }
  }
  
  
  
  
  
  segm<-read.delim("path_to_your/4b2o_F2.geno.seg.filt.score.txt")
  segm$Marker<-paste0(segm$Chr,"_",segm$Pos)
  ajpeak<-read.delim("path_to_your/peaks.adjust.txt")
  sig.m<-ajpeak$GBS_seg_marker
  
  
  ###plot by individuals
  ### plots for each individual
  chrsize <- read.delim(chrsizefile,stringsAsFactors = F, header=F)
  colnames(chrsize)<-c("Chr","Size")
  chrsize<-chrsize[chrsize$Chr %in% chr.names,]
  chrsize$Size<-as.numeric(chrsize$Size)
  
  for (j in 1:10) { # 10 chromosomes
    if (out.fmt=="pdf"){
      pdf(paste0(out.prefix, j, ".pdf"), width=width, height=height,  pointsize = 4)
    }
    if (out.fmt=="png"){
      png(paste0(out.prefix, j, ".png"),width=width, height=height,   pointsize = 4)
    }
    ####need to be edited
    
      
    par(mar = c(0, 5, 4, 1),mfrow=c(2,1))
    
    plot(NULL,NULL, ylim=c(0, 1*1.01), xlim=c(0, 1.2*chrsize[chrsize$Chr==j,"Size"]/1000000),
          type="n", cex.main=cex.main, main=paste0(main.lab, j),xaxt="n",bty="n",
         xlab="", ylab="genotype frequency", cex.axis = cex.axis, cex.lab = cex.axis) #
    legend(chrsize[chrsize$Chr==j,"Size"]/1000000*1.003,1.05,
           legend = c("A188","B73","Hetero"),col = geno.col,lty=1, cex=cex.axis*0.8)
    

      eseg<- geno.seg1[geno.seg1$Chr ==paste0("", j,"") , ]
     
      ### seg draw:
      
      seg.rects(segdata  = eseg,  segcol = geno.col)
     
      
      #axis(2, at =  shift, labels = "seg",pos=(1), cex.axis = 1,las=2,tick = F); # left side of axis
      
 
    ### label
      a<-colnames(geno.seg1)[ncol(geno.seg1)]
      a<-gsub('.[A-Za-z]+$',"",a)
    text(x = chrsize[chrsize$Chr==j,"Size"]/1000000*1.1, y = 0.5, 
         labels = a, cex = cex.axis,pos = 2, xpd = T)
    if(j %in% ajpeak$Chromosome){
      rect(xleft=ajpeak[ajpeak$Chromosome==j,"pm.start"]/1000000,ybottom = -3,ytop =1.2,xright =ajpeak[ajpeak$Chromosome==j,"pm.end"]/1000000, lty = 2,border = "red",col = NA  )
      abline(v=ajpeak[ajpeak$Chromosome==j,"peak.pos"]/1000000,col="gray")
    }
    
    par(mar = c(4, 5,0, 1))
    
    
    
    plot(NULL,NULL, ylim=c(0, 1.01*1), xlim=c(0, 1.2*chrsize[chrsize$Chr==j,"Size"]/1000000),
      type="n", cex.main=cex.main, main="", bty="n",
         xlab="Physical coordinate (Mb)", ylab="genotype frequency", cex.axis = cex.axis, cex.lab = cex.axis)
    # legend(300,10,legend = geno.code,col = geno.col,lty=1, cex=0.8)
    
    
    eseg<- geno.seg2[geno.seg2$Chr ==paste0("", j,"") , ]
    
    ### seg draw:
    
    seg.rects(segdata  = eseg,  segcol = geno.col)
      
      #axis(2, at =  shift, labels = "seg",pos=(1), cex.axis = 1,las=2,tick = F); # left side of axis
      
      ### label
    a<-colnames(geno.seg2)[ncol(geno.seg2)]
    a<-gsub('.[A-Za-z]+$',"",a)
    text(x = chrsize[chrsize$Chr==j,"Size"]/1000000*1.1,
         y = 0.5, labels = a, cex = cex.axis,pos = 2, xpd = T)
    if(j %in% ajpeak$Chromosome){
      rect(xleft=ajpeak[ajpeak$Chromosome==j,"pm.start"]/1000000,ybottom = 0-1,ytop =1.2,xright =ajpeak[ajpeak$Chromosome==j,"pm.end"]/1000000, lty = 2,border = "red",col = NA  )
      abline(v=ajpeak[ajpeak$Chromosome==j,"peak.pos"]/1000000,col="gray")
    }
    dev.off()
  }
  
}







