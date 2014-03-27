options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

library(DESeq2)
path<-dirname(args[1])
counts<-read.table(args[1],header=T,sep="\t")
row.names(counts)<-counts[,1]
gen<-read.table(args[2],header=T,sep="\t",row.names=1)
sizefactors<-read.table(args[3])
sorting<-read.table(args[4],row.names=1,sep="\t")

ns<-ncol(counts)
counts<-counts[,7:ns]
counts<-counts[apply(counts,1,sum)>ns,]
names(counts)<-as.character(row.names(sorting))

gen<-gen[names(counts),]
conditions<-as.character(gen[,args[5]])
design<-data.frame(con=conditions,row.names=row.names(gen))

dse<-DESeqDataSetFromMatrix(counts,design,design=~con)
dse<-estimateSizeFactors(dse)
sizeFactors(dse)<-sizefactors[,1]
norm<-round(counts(dse,normalize=TRUE))
colnames(norm)<-design[,1]

write.table(norm[,colnames(norm)=="INV",drop=FALSE],paste0(args[1],".norm1"),col.names=F,quote=F,sep="\t")
write.table(norm[,colnames(norm)=="STD",drop=FALSE],paste0(args[1],".norm2"),col.names=F,quote=F,sep="\t")
