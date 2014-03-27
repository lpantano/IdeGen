options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

library(DESeq2)
path<-dirname(args[1])
ref<-read.table(args[1],header=T,sep="\t")
ns<-ncol(ref)
counts<-ref[,7:ns]
ns<-ns-7+1
counts<-counts[apply(counts,1,sum)>ns,]
fakeconditions<-c(rep(0,round((ns-1)/2)),rep(1,round(ns/2)))
design<-data.frame(con=fakeconditions,row.names=names(counts))
dse<-DESeqDataSetFromMatrix(counts,design,design=~con)
dse<-estimateSizeFactors(dse)
write.table(sizeFactors(dse),paste0(path,"/sizefactors"))
#norm<-counts(dse,normalize=TRUE)
#norm.round<-round(norm)
