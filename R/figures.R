options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
path<-dirname(args[1])

calcyu1<-function(e){
  if (dim(e)[1]>1){
   return (e$V6 %*% diag(1/(apply(e[,5:10],1,max)))*30+40)
  }else{
   return ((e$V6/(apply(e[,5:10],1,max)))*30+40)
  }
}

calcyd1<-function(e){
  if (dim(e)[1]>1){
    return (e$V7 %*% diag(1/(apply(e[,5:10],1,max)))*30+40)
  }else{
    return ((e$V7/(apply(e[,5:10],1,max)))*30+40)
  }
}

calcyu2<-function(e){
  if (dim(e)[1]>1){
    return (e$V9 %*% diag(1/(apply(e[,5:10],1,max)))*30+40)
  }else{
    return ((e$V9/(apply(e[,5:10],1,max)))*30+40)
  }
}

calcyd2<-function(e){
  if (dim(e)[1]>1){
    return (e$V10 %*% diag(1/(apply(e[,5:10],1,max)))*30+40)
  }else{
    return ((e$V10/(apply(e[,5:10],1,max)))*30+40)
  }
}

calcymu1<-function(e){
  if (dim(e)[1]>1){
    return (e$V5 %*% diag(1/(apply(e[,5:10],1,max)))*30+40)
  }else{
    return ((e$V5/(apply(e[,5:10],1,max)))*30+40)
  }
}

calcymu2<-function(e){
  if (dim(e)[1]>1){
    return (e$V8 %*% diag(1/(apply(e[,5:10],1,max)))*30+40)
  }else{
    return ((e$V8/(apply(e[,5:10],1,max)))*30+40)
  }
}

#file<-read.table("~/crickshared/Bioinformatics/iRNASeq/geovadis/fusion/HsInv0379/svdg.rda",sep="\t")
file<-read.table(args[1],sep="\t")
genes<-unlist(unique(file$V1))

for (g in genes){
  print(g)
  t<-file[file$V1==g,]
  t<-t[order(t$V2),]
  e<-t[t$V3=="exon",]
  i<-t[t$V3!="exon",]
  ne<-length(e$V2)
  i<-i[1:(ne-1),]
  exonx1<-seq(1,length(e$V2)*10,10)+1
  exonx2<-seq(1,length(e$V2)*10,10)+6
  area<-exonx2-exonx1
  
  exonq1<-exonx1+1.5
  exonq2<-exonx1+3.5
  yu1<-calcyu1(e)
  yd1<-calcyd1(e)
  yu2<-calcyu2(e)
  yd2<-calcyd2(e)
  ymu1<-calcymu1(e)
  ymu2<-calcymu2(e)
  known<-exonx1+e$V4*area
  intronx1<-exonx2[1:(ne-1)]
  intronx2<-exonx1[2:ne]
  intronm<-intronx1+2.5
  pdf(paste0(args[2],"/",g,".pdf"))
  par(mar=c(1,1,3,1))
  plot(1,type='n',ylim=c(0,90),xlim=c(0,nrow(t)*5),
       axes=FALSE,xaxt = "n",yaxt="n",xlab="",ylab="",main=g)  
  #exons
  rect(exonx1,20,exonx2,30)  
  rect(exonx1,20,known,30,col="azure3") 
  if (ne>1){
   #introns
   cols<-sub("NO","black",i$V11)
   cols<-sub("YES","darkred",cols)
   segments(intronx1,25,intronx2,25,col=cols,lwd=1.5) 
   #intron quantity
   text(intronm,25,labels=round(i$V4,digits=3),cex=.7,pos=3)
  }
  #exons quantity
  segments(exonq1,yu1,exonq1,yd1,col="orange3",lwd=1.2)
  segments(exonq2,yu2,exonq2,yd2,col="green3",,lwd=1.2)
  points(exonq1,ymu1,pch=20,col="orange3")
  points(exonq2,ymu2,pch=20,col="green3")
  
  #legend
  text(1,3,labels="grey means % of overlap with known exons",cex=.8,pos=4)
  text(1,8,labels="green:control,orange:variation",cex=.8,pos=4)
  text(1,13,labels="red intron means intron in BP. Numbers above are average reads per sample in the variation group",cex=.8,pos=4)
  dev.off()
}
