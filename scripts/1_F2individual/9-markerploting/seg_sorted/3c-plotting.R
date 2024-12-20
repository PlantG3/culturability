setwd(".")

chrsize <- read.delim("path_to_your/chr.size.refgen4.txt")

seg<-read.delim("2o-seg.marker.geno.chisq.txt",stringsAsFactors = F)

colnames(seg)
geno.seg<-seg
geno.seg1<-geno.seg[,c(2:4,10:12)]
geno.seg2<-geno.seg[,c(2:4,17:19)]

library(RColorBrewer)

source("3c-plotting_sourcecode.R")

segplot.chr(segfile="2o-seg.marker.geno.chisq.txt",
            chrsizefile="path_to_your/chr.size.refgen4.txt",
            chr.names=1:10, out.prefix="3o.Genotype_segregation_of_chrosome ", geno.code=c("2","1","3"),
            geno.col=brewer.pal(n = 8, name = "Dark2")[c(1,2,3)],main.lab="Genotype segregation of chrosome ",
            out.fmt="pdf",width=6,height=4,cex.main=2,cex.axis = 2, cex.lab = 2)



###############################################################################################
## plot in 1 png
###old
seg<-read.delim("2o-seg.marker.geno.chisq.txt",stringsAsFactors = F)
seg<-seg[,c(1:5,14:21,6:13,22,23,24)]
XT.II<-seg[,c(1:5,6:13,22,24)]
XT.I<-seg[,c(1:5,14:22),23]
colnames(seg)

pdf(paste0("3o-all-chr.seg.distribution.pdf"), width=13, height=10, pointsize = 2)####need to be edited
#par(mar = c(5, 4, 3, 2),mfrow=c(2,5),oma=c(0,0,0,0)) 


m<-matrix(c(1,1,1,1,1,2:11),nrow = 3,byrow = TRUE)
layout(mat=m,heights=c(0.1, 0.45, 0.45))

par(mar = c(0.5, 3, 4, 2))
plot(1, type = "n", axes=FALSE, xlab="", ylab="",main="Genotype Distribution of B73/A188 F2 Individuals",cex.main=4)
#plot_colors <- c("blue","black", "green", "orange", "pink")
#legend(x = "top",inset = 0, legend = c("A188","Hetero","B73","Chromosome Size"),pch=c(15,15,15,15), xjust=0.5,#yjust = 1,
#      col=c("limegreen","gray87","blue"),horiz=TRUE, pt.cex = 2,cex = 0.8,box.lty = 0)
legend(x = "center",inset = 0, legend = c("A188","B73","Heterozygous/Both"),#yjust = 1,
       col=c("limegreen","blue","gray87"),horiz=TRUE, lwd=5, cex=1,border = F)
#legend(x = "top",inset = 0,
#      legend = c("Fabricated Metal", "Iron and Steel", "Paper","Beverages", "Tobacco"), 
#col=plot_colors,  horiz = TRUE)


par(mar = c(2, 5, 3, 0))

j<-1
plot(NULL,NULL, ylim=c(0,270), xlim=c(0, chrsize$Size[chrsize$ChrNum==j]/1000000),
     bty="n", type="n",yaxt="n", cex.main=2.0, main=paste0("Chr",j), 
     xlab="", ylab="", cex.axis = 1, cex.lab = 1) 


dp<-seg[seg$Chr==j,]
#
p<-dp$pvalue
bin<-dp$Pos
#-log10(0.004)
#-log10(p[which.min(p)])
points(bin/1000000, 200+80*(-log10(p)/1.5), pch=20, cex =0.4)
lines(x=c(-2,-2),y=c(200,280),lwd=0.5)
lines(x=c(-2,-1),y=c(200,200),lwd=0.5)
lines(x=c(-2,-1),y=c(200+80,200+80),lwd=0.5)
lines(x=c(-2,-1),y=c(200+2*80,200+2*80),lwd=0.5)
lines(x=c(-2,-1),y=c(200+3*80,200+3*80),lwd=0.5)
lines(x=c(-2,-1),y=c(200+4*80,200+4*80),lwd=0.5)
lines(x=c(-2,-1),y=c(200+5*80,200+5*80),lwd=0.5)

lines(x=c(0,chrsize$Size[chrsize$ChrNum==j]/1000000),y=c(200+80*(-log10(0.05)),200+80*(-log10(0.05))), col="red", lwd=0.5,lty=2)
text(x =-5.5,y=200+2,"0",cex=0.7)
text(x =-5.5,y=200+80,"1",cex=0.7)
text(x =-5.5,y=200+2*80,"2",cex=0.7)
text(x =-5.5,y=200+3*80,"3",cex=0.7)
text(x =-5.5,y=200+4*80,"4",cex=0.7)
text(x =-5.5,y=200+5*80-2,"5",cex=0.7)
#text(x=-15,y=200+80,"-log10(p)/1.5",cex=1,srt=90)
axis(2, at =200+2*80, labels = "-log10(p)/1.5" ,pos=(2), cex.axis = 1,las=2,tick = F)

count<-0
for (s in c("XT.II","XT.I")){####Plot by callus types
  #for (s in c("seg_I","seg_II")){####Plot by callus types
  # dwinseg<-s  
  
  count<-count+1
  sy<-200-100*count
  startcol<-6*count+2*(count-1)
  endcol<-6*count+2*(count-1)+7
  dwinseg<-seg[,c(1:5,startcol:endcol)]
  
  dwinsegechr<-dwinseg[dwinseg$Chr==j, ]
  head(dwinsegechr)
  for (n in dwinsegechr$Start){
    N<-dwinsegechr[dwinsegechr$Start==n,]$End
    b<-dwinsegechr[dwinsegechr$Start==n,][,11]
    h<-dwinsegechr[dwinsegechr$Start==n,][,12]
    a<-dwinsegechr[dwinsegechr$Start==n,][,10]
    
    #####draw line in the middle of the window,maybe need to use box later.
    #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry,chry+20*h),col="red")
    #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry+20*h,chry+20*h+20*b),col="blue")
    #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry+20*h+20*b,chry+20*h+20*b+20*a),col="dark green")
    
    rect(n/1000000,sy,N/1000000,sy+80*h,col="gray87",border = "gray87")
    rect(n/1000000,sy+80*h,N/1000000,sy+80*b+80*h,col="blue",border = "blue")
    h+a+b
    rect(n/1000000,sy+80*h+80*b,N/1000000,sy+80*h+80*b+80*a,col="limegreen",border = "limegreen")      
    axis(2, at =sy+40, labels = s ,pos=(2), cex.axis = 1,las=2,tick = F); # left side of axis
    lines(x=c(0, chrsize$Size[chrsize$ChrNum==j]/1000000),y=c(sy+40,sy+40),lty=2,lwd=0.5 )
    lines(x=c(-2,-2),y=c(sy,sy+80),lwd=0.5)
    lines(x=c(-2,-1),y=c(sy,sy),lwd=0.5)
    lines(x=c(-2,-1),y=c(sy+80,sy+80),lwd=0.5)
    #text(x =-5.5,y=sy,"0",cex=0.7)
    #text(x =-5.5,y=sy+80,"1",cex=0.7)
  } 
  
}


for (j in 2:5){
  #par(mar = c(5, 8, 3, 2),mfrow=c(1,1)) 
  #j<-9
  par(mar = c(2, 0, 3, 0))
  plot(NULL,NULL, ylim=c(0,270), xlim=c(0, chrsize$Size[chrsize$ChrNum==j]/1000000),
       bty="n", type="n",yaxt="n", cex.main=2.0, main=paste0("Chr",j), 
       xlab="", ylab="", cex.axis = 1, cex.lab = 1) 
  
  #legend(280,350,  c("A188","Hetero","B73"),pch=c(15,15,15), xjust=0,#yjust = 1,
  #      col=c("limegreen","gray87","blue"),horiz=TRUE,
  #     pt.cex = 1,cex = 0.5,box.lty = 0)
  #j<-1
  dp<-seg[seg$Chr==j,]
  p<-dp$pvalue
  bin<-dp$Pos
  #-log10(0.004)
  #-log10(p[which.min(p)])
  points((bin)/1000000, 200+80*(-log10(p)/1.5/1.5), pch=20, cex =0.4)
  lines(x=c(-2,-2),y=c(200,280),lwd=0.5)
  lines(x=c(-2,-1),y=c(200,200),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+80,200+80),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+2*80,200+2*80),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+3*80,200+3*80),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+4*80,200+4*80),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+5*80,200+5*80),lwd=0.5)
  
  lines(x=c(0,chrsize$Size[chrsize$ChrNum==j]/1000000),y=c(200+80*(-log10(0.05)),200+80*(-log10(0.05))), col="red", lwd=0.5,lty=2)
  
  text(x =-5.5,y=200+2,"0",cex=0.7)
  text(x =-5.5,y=200+80,"1",cex=0.7)
  text(x =-5.5,y=200+2*80,"2",cex=0.7)
  text(x =-5.5,y=200+3*80,"3",cex=0.7)
  text(x =-5.5,y=200+4*80,"4",cex=0.7)
  text(x =-5.5,y=200+5*80-2,"5",cex=0.7)
  #text(x=-15,y=200+80,"-log10(p)/1.5",cex=1,srt=90)
  #axis(2, at =200+2*80, labels = "-log10(p)/1.5" ,pos=(2), cex.axis = 0.8,las=2,tick = F)
  
  count<-0
  for (s in c("XT.II","XT.I")){####Plot by callus types
    #for (s in c("seg_I","seg_II")){####Plot by callus types
    # dwinseg<-s  
    count<-count+1
    sy<-200-100*count
    startcol<-6*count+2*(count-1)
    endcol<-6*count+2*(count-1)+7
    dwinseg<-seg[,c(1:5,startcol:endcol)]
    
    dwinsegechr<-dwinseg[dwinseg$Chr==j, ]
    head(dwinsegechr)
    for (n in dwinsegechr$Start){
      N<-dwinsegechr[dwinsegechr$Start==n,]$End
      b<-dwinsegechr[dwinsegechr$Start==n,][,11]
      h<-dwinsegechr[dwinsegechr$Start==n,][,12]
      a<-dwinsegechr[dwinsegechr$Start==n,][,10]
      #####draw line in the middle of the window,maybe need to use box later.
      #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry,chry+20*h),col="red")
      #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry+20*h,chry+20*h+20*b),col="blue")
      #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry+20*h+20*b,chry+20*h+20*b+20*a),col="dark green")
      
      rect(n/1000000,sy,N/1000000,sy+80*h,col="gray87",border = "gray87")
      rect(n/1000000,sy+80*h,N/1000000,sy+80*b+80*h,col="blue",border = "blue")
      rect(n/1000000,sy+80*h+80*b,N/1000000,sy+80*h+80*b+80*a,col="limegreen",border = "limegreen")      
      #axis(2, at =sy+40, labels = s ,pos=(2), cex.axis = 1,las=2,tick = F); # left side of axis
      lines(x=c(0, chrsize$Size[chrsize$ChrNum==j]/1000000),y=c(sy+40,sy+40),lty=2,lwd=0.5 )
      lines(x=c(-2,-2),y=c(sy,sy+80),lwd=0.5)
      lines(x=c(-2,-1),y=c(sy,sy),lwd=0.5)
      lines(x=c(-2,-1),y=c(sy+80,sy+80),lwd=0.5)
      #text(x =-5.5,y=sy,"0",cex=0.7)
      #text(x =-5.5,y=sy+80,"1",cex=0.7)
    } 
    
  }
  
}

j<-6

par(mar = c(5, 5, 6, 0))
plot(NULL,NULL, ylim=c(0,270), xlim=c(0, chrsize$Size[chrsize$ChrNum==j]/1000000),
     bty="n", type="n",yaxt="n", cex.main=2.0, main=paste0("Chr",j), 
     xlab="", ylab="", cex.axis = 1, cex.lab = 1) 

dp<-seg[seg$Chr==j,]

p<-dp$pvalue
bin<-dp$Pos
#-log10(0.004)
#-log10(p[which.min(p)])
points((bin)/1000000, 200+80*(-log10(p)/1.5), pch=20, cex =0.4)
lines(x=c(-2,-2),y=c(200,280),lwd=0.5)
lines(x=c(-2,-1),y=c(200,200),lwd=0.5)
lines(x=c(-2,-1),y=c(200+80,200+80),lwd=0.5)
lines(x=c(-2,-1),y=c(200+2*80,200+2*80),lwd=0.5)
lines(x=c(-2,-1),y=c(200+3*80,200+3*80),lwd=0.5)
lines(x=c(-2,-1),y=c(200+4*80,200+4*80),lwd=0.5)
lines(x=c(-2,-1),y=c(200+5*80,200+5*80),lwd=0.5)
lines(x=c(0,chrsize$Size[chrsize$ChrNum==j]/1000000),y=c(200+80*(-log10(0.05)),200+80*(-log10(0.05))), col="red", lwd=0.5,lty=2)

text(x =-5.5,y=200+2,"0",cex=0.7)
text(x =-5.5,y=200+80,"1",cex=0.7)
text(x =-5.5,y=200+2*80,"2",cex=0.7)
text(x =-5.5,y=200+3*80,"3",cex=0.7)
text(x =-5.5,y=200+4*80,"4",cex=0.7)
text(x =-5.5,y=200+5*80-2,"5",cex=0.7)
#text(x=-15,y=200+80,"-log10(p)/1.5",cex=1,srt=90)
axis(2, at =200+2*80, labels = "-log10(p)/1.5" ,pos=(2), cex.axis = 1,las=2,tick = F)

count<-0
for (s in c("XT.II","XT.I")){####Plot by callus types
  #for (s in c("seg_I","seg_II")){####Plot by callus types
  # dwinseg<-s  
  count<-count+1
  sy<-200-100*count
  startcol<-6*count+2*(count-1)
  endcol<-6*count+2*(count-1)+7
  dwinseg<-seg[,c(1:5,startcol:endcol)]
  
  dwinsegechr<-dwinseg[dwinseg$Chr==j, ]
  head(dwinsegechr)
  for (n in dwinsegechr$Start){
    N<-dwinsegechr[dwinsegechr$Start==n,]$End
    b<-dwinsegechr[dwinsegechr$Start==n,][,11]
    h<-dwinsegechr[dwinsegechr$Start==n,][,12]
    a<-dwinsegechr[dwinsegechr$Start==n,][,10]
    #####draw line in the middle of the window,maybe need to use box later.
    #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry,chry+20*h),col="red")
    #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry+20*h,chry+20*h+20*b),col="blue")
    #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry+20*h+20*b,chry+20*h+20*b+20*a),col="dark green")
    
    rect(n/1000000,sy,N/1000000,sy+80*h,col="gray87",border = "gray87")
    rect(n/1000000,sy+80*h,N/1000000,sy+80*b+80*h,col="blue",border = "blue")
    rect(n/1000000,sy+80*h+80*b,N/1000000,sy+80*h+80*b+80*a,col="limegreen",border = "limegreen")      
    axis(2, at =sy+40, labels = s ,pos=(2), cex.axis = 1,las=2,tick = F); # left side of axis
    lines(x=c(0, chrsize$Size[chrsize$ChrNum==j]/1000000),y=c(sy+40,sy+40),lty=2,lwd=0.5 )
    lines(x=c(-2,-2),y=c(sy,sy+80),lwd=0.5)
    lines(x=c(-2,-1),y=c(sy,sy),lwd=0.5)
    lines(x=c(-2,-1),y=c(sy+80,sy+80),lwd=0.5)
    #text(x =-5.5,y=sy,"0",cex=0.7)
    #text(x =-5.5,y=sy+80,"1",cex=0.7)
  } 
  
}

par(mar = c(5, 0, 6, 0))
for (j in 7){
  #par(mar = c(5, 8, 3, 2),mfrow=c(1,1)) 
  #j<-9
  
  plot(NULL,NULL, ylim=c(0,270), xlim=c(0, chrsize$Size[chrsize$ChrNum==j]/1000000),
       bty="n", type="n",yaxt="n", cex.main=2.0, main=paste0("Chr",j), 
       xlab="", ylab="", cex.axis = 1, cex.lab = 1) 
  
  #legend(280,350,  c("A188","Hetero","B73"),pch=c(15,15,15), xjust=0,#yjust = 1,
  #      col=c("limegreen","gray87","blue"),horiz=TRUE,
  #     pt.cex = 1,cex = 0.5,box.lty = 0)
  #j<-1
  dp<-seg[seg$Chr==j,]
  
  p<-dp$pvalue
  bin<-dp$Pos
  #-log10(0.004)
  #-log10(p[which.min(p)])
  points((bin)/1000000, 200+80*(-log10(p)/1.5), pch=20, cex =0.4)
  lines(x=c(-2,-2),y=c(200,280),lwd=0.5)
  lines(x=c(-2,-1),y=c(200,200),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+80,200+80),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+2*80,200+2*80),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+3*80,200+3*80),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+4*80,200+4*80),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+5*80,200+5*80),lwd=0.5)
  lines(x=c(0,chrsize$Size[chrsize$ChrNum==j]/1000000),y=c(200+80*(-log10(0.05)),200+80*(-log10(0.05))), col="red", lwd=0.5,lty=2)
  
  text(x =-5.5,y=200+2,"0",cex=0.7)
  text(x =-5.5,y=200+80,"1",cex=0.7)
  text(x =-5.5,y=200+2*80,"2",cex=0.7)
  text(x =-5.5,y=200+3*80,"3",cex=0.7)
  text(x =-5.5,y=200+4*80,"4",cex=0.7)
  text(x =-5.5,y=200+5*80-2,"5",cex=0.7)
  #text(x=-15,y=200+80,"-log10(p)/1.5",cex=1,srt=90)
  # axis(2, at =200+2*80, labels = "-log10(p)/1.5" ,pos=(2), cex.axis = 0.8,las=2,tick = F)
  
  count<-0
  for (s in c("XT.II","XT.I")){####Plot by callus types
    #for (s in c("seg_I","seg_II")){####Plot by callus types
    # dwinseg<-s  
    count<-count+1
    sy<-200-100*count
    startcol<-6*count+2*(count-1)
    endcol<-6*count+2*(count-1)+7
    dwinseg<-seg[,c(1:5,startcol:endcol)]
    
    dwinsegechr<-dwinseg[dwinseg$Chr==j, ]
    head(dwinsegechr)
    for (n in dwinsegechr$Start){
      N<-dwinsegechr[dwinsegechr$Start==n,]$End
      b<-dwinsegechr[dwinsegechr$Start==n,][,11]
      h<-dwinsegechr[dwinsegechr$Start==n,][,12]
      a<-dwinsegechr[dwinsegechr$Start==n,][,10]
      #####draw line in the middle of the window,maybe need to use box later.
      #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry,chry+20*h),col="red")
      #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry+20*h,chry+20*h+20*b),col="blue")
      #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry+20*h+20*b,chry+20*h+20*b+20*a),col="dark green")
      
      rect(n/1000000,sy,N/1000000,sy+80*h,col="gray87",border = "gray87")
      rect(n/1000000,sy+80*h,N/1000000,sy+80*b+80*h,col="blue",border = "blue")
      rect(n/1000000,sy+80*h+80*b,N/1000000,sy+80*h+80*b+80*a,col="limegreen",border = "limegreen")      
      # axis(2, at =sy+40, labels = s ,pos=(2), cex.axis = 1,las=2,tick = F); # left side of axis
      lines(x=c(0, chrsize$Size[chrsize$ChrNum==j]/1000000),y=c(sy+40,sy+40),lty=2,lwd=0.5 )
      lines(x=c(-2,-2),y=c(sy,sy+80),lwd=0.5)
      lines(x=c(-2,-1),y=c(sy,sy),lwd=0.5)
      lines(x=c(-2,-1),y=c(sy+80,sy+80),lwd=0.5)
      #text(x =-5.5,y=sy,"0",cex=0.7)
      #text(x =-5.5,y=sy+80,"1",cex=0.7)
    } 
    
  }
  
}


j<-8
plot(NULL,NULL, ylim=c(0,270), xlim=c(0, chrsize$Size[chrsize$ChrNum==j]/1000000),
     bty="n", type="n",yaxt="n", cex.main=2.0, main=paste0("Chr",j), 
     xlab="Physical coordinate (Mb)", ylab="", cex.axis = 1, cex.lab = 2) 

#legend(280,350,  c("A188","Hetero","B73"),pch=c(15,15,15), xjust=0,#yjust = 1,
#      col=c("limegreen","gray87","blue"),horiz=TRUE,
#     pt.cex = 1,cex = 0.5,box.lty = 0)
#j<-1
dp<-seg[seg$Chr==j,]

p<-dp$pvalue
bin<-dp$Pos
#-log10(0.004)
#-log10(p[which.min(p)])
points((bin)/1000000, 200+80*(-log10(p)/1.5), pch=20, cex =0.4)
lines(x=c(-2,-2),y=c(200,280),lwd=0.5)
lines(x=c(-2,-1),y=c(200,200),lwd=0.5)
lines(x=c(-2,-1),y=c(200+80,200+80),lwd=0.5)
lines(x=c(-2,-1),y=c(200+2*80,200+2*80),lwd=0.5)
lines(x=c(-2,-1),y=c(200+3*80,200+3*80),lwd=0.5)
lines(x=c(-2,-1),y=c(200+4*80,200+4*80),lwd=0.5)
lines(x=c(-2,-1),y=c(200+5*80,200+5*80),lwd=0.5)
lines(x=c(0,chrsize$Size[chrsize$ChrNum==j]/1000000),y=c(200+80*(-log10(0.05)),200+80*(-log10(0.05))), col="red", lwd=0.5,lty=2)

text(x =-5.5,y=200+2,"0",cex=0.7)
text(x =-5.5,y=200+80,"1",cex=0.7)
text(x =-5.5,y=200+2*80,"2",cex=0.7)
text(x =-5.5,y=200+3*80,"3",cex=0.7)
text(x =-5.5,y=200+4*80,"4",cex=0.7)
text(x =-5.5,y=200+5*80-2,"5",cex=0.7)
#text(x=-15,y=200+80,"-log10(p)/1.5",cex=1,srt=90)
#axis(2, at =200+2*80, labels = "-log10(p)/1.5" ,pos=(2), cex.axis = 0.8,las=2,tick = F)

count<-0
for (s in c("XT.II","XT.I")){####Plot by callus types
  #for (s in c("seg_I","seg_II")){####Plot by callus types
  # dwinseg<-s  
  count<-count+1
  sy<-200-100*count
  startcol<-6*count+2*(count-1)
  endcol<-6*count+2*(count-1)+7
  dwinseg<-seg[,c(1:5,startcol:endcol)]
  
  dwinsegechr<-dwinseg[dwinseg$Chr==j, ]
  head(dwinsegechr)
  for (n in dwinsegechr$Start){
    N<-dwinsegechr[dwinsegechr$Start==n,]$End
    b<-dwinsegechr[dwinsegechr$Start==n,][,11]
    h<-dwinsegechr[dwinsegechr$Start==n,][,12]
    a<-dwinsegechr[dwinsegechr$Start==n,][,10]
    #####draw line in the middle of the window,maybe need to use box later.
    #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry,chry+20*h),col="red")
    #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry+20*h,chry+20*h+20*b),col="blue")
    #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry+20*h+20*b,chry+20*h+20*b+20*a),col="dark green")
    
    rect(n/1000000,sy,N/1000000,sy+80*h,col="gray87",border = "gray87")
    rect(n/1000000,sy+80*h,N/1000000,sy+80*b+80*h,col="blue",border = "blue")
    rect(n/1000000,sy+80*h+80*b,N/1000000,sy+80*h+80*b+80*a,col="limegreen",border = "limegreen")      
    #   axis(2, at =sy+40, labels = s ,pos=(2), cex.axis = 1,las=2,tick = F); # left side of axis
    lines(x=c(0, chrsize$Size[chrsize$ChrNum==j]/1000000),y=c(sy+40,sy+40),lty=2,lwd=0.5 )
    lines(x=c(-2,-2),y=c(sy,sy+80),lwd=0.5)
    lines(x=c(-2,-1),y=c(sy,sy),lwd=0.5)
    lines(x=c(-2,-1),y=c(sy+80,sy+80),lwd=0.5)
    #text(x =-5.5,y=sy,"0",cex=0.7)
    #text(x =-5.5,y=sy+80,"1",cex=0.7)
  } 
  
}

for (j in 9:10){
  #par(mar = c(5, 8, 3, 2),mfrow=c(1,1)) 
  #j<-9
  
  plot(NULL,NULL, ylim=c(0,270), xlim=c(0, chrsize$Size[chrsize$ChrNum==j]/1000000),
       bty="n", type="n",yaxt="n", cex.main=2.0, main=paste0("Chr",j), 
       xlab="", ylab="", cex.axis = 1, cex.lab = 1) 
  
  #legend(280,350,  c("A188","Hetero","B73"),pch=c(15,15,15), xjust=0,#yjust = 1,
  #      col=c("limegreen","gray87","blue"),horiz=TRUE,
  #     pt.cex = 1,cex = 0.5,box.lty = 0)
  #j<-1
  dp<-seg[seg$Chr==j,]
  
  p<-dp$pvalue
  bin<-dp$Pos
  #-log10(0.004)
  #-log10(p[which.min(p)])
  points((bin)/1000000, 200+80*(-log10(p)/1.5), pch=20, cex =0.4)
  lines(x=c(-2,-2),y=c(200,280),lwd=0.5)
  lines(x=c(-2,-1),y=c(200,200),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+80,200+80),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+2*80,200+2*80),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+3*80,200+3*80),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+4*80,200+4*80),lwd=0.5)
  lines(x=c(-2,-1),y=c(200+5*80,200+5*80),lwd=0.5)
  lines(x=c(0,chrsize$Size[chrsize$ChrNum==j]/1000000),y=c(200+80*(-log10(0.05)),200+80*(-log10(0.05))), col="red", lwd=0.5,lty=2)
  
  text(x =-5.5,y=200+2,"0",cex=0.7)
  text(x =-5.5,y=200+80,"1",cex=0.7)
  text(x =-5.5,y=200+2*80,"2",cex=0.7)
  text(x =-5.5,y=200+3*80,"3",cex=0.7)
  text(x =-5.5,y=200+4*80,"4",cex=0.7)
  text(x =-5.5,y=200+5*80-2,"5",cex=0.7)
  #text(x=-15,y=200+80,"-log10(p)/1.5",cex=1,srt=90)
  #axis(2, at =200+2*80, labels = "-log10(p)/1.5" ,pos=(2), cex.axis = 0.8,las=2,tick = F)
  
  count<-0
  for (s in c("XT.II","XT.I")){####Plot by callus types
    #for (s in c("seg_I","seg_II")){####Plot by callus types
    # dwinseg<-s  
    count<-count+1
    sy<-200-100*count
    startcol<-6*count+2*(count-1)
    endcol<-6*count+2*(count-1)+7
    dwinseg<-seg[,c(1:5,startcol:endcol)]
    
    dwinsegechr<-dwinseg[dwinseg$Chr==j, ]
    head(dwinsegechr)
    for (n in dwinsegechr$Start){
      N<-dwinsegechr[dwinsegechr$Start==n,]$End
      b<-dwinsegechr[dwinsegechr$Start==n,][,11]
      h<-dwinsegechr[dwinsegechr$Start==n,][,12]
      a<-dwinsegechr[dwinsegechr$Start==n,][,10]
      #####draw line in the middle of the window,maybe need to use box later.
      #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry,chry+20*h),col="red")
      #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry+20*h,chry+20*h+20*b),col="blue")
      #lines(x = c((n+5000)/1000000,(n+5000)/1000000),y = c(chry+20*h+20*b,chry+20*h+20*b+20*a),col="dark green")
      
      rect(n/1000000,sy,N/1000000,sy+80*h,col="gray87",border = "gray87")
      rect(n/1000000,sy+80*h,N/1000000,sy+80*b+80*h,col="blue",border = "blue")
      rect(n/1000000,sy+80*h+80*b,N/1000000,sy+80*h+80*b+80*a,col="limegreen",border = "limegreen")      
      #axis(2, at =sy+40, labels = s ,pos=(2), cex.axis = 1,las=2,tick = F); # left side of axis
      lines(x=c(0, chrsize$Size[chrsize$ChrNum==j]/1000000),y=c(sy+40,sy+40),lty=2,lwd=0.5 )
      lines(x=c(-2,-2),y=c(sy,sy+80),lwd=0.5)
      lines(x=c(-2,-1),y=c(sy,sy),lwd=0.5)
      lines(x=c(-2,-1),y=c(sy+80,sy+80),lwd=0.5)
      #text(x =-5.5,y=sy,"0",cex=0.7)
      #text(x =-5.5,y=sy+80,"1",cex=0.7)
    } 
    
  }
  
}
dev.off()






