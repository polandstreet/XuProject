
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome")
biocLite("BiocGenerics")
biocLite("BSgenome.Hsapiens.UCSC.hg38")

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd("C:/PKU_research/XuProject")

datafile<-read.table("SpliceRef.txt",
                               header=TRUE,sep="",fill=TRUE)


chr<-as.matrix(datafile[,"chrom"])

strand<-as.matrix(datafile[,"strand"]	)

startpos<-as.matrix(datafile[,"exonStarts"])

name<-as.matrix(datafile[,"name"])

name2<-as.matrix(datafile[,"name2"])

exonArray<-matrix("",ncol=5,nrow=sum(as.numeric(datafile[,"exonCount"])))

line<-1

for(i in 1:length(startpos)){
  
  exonStartPos<-unlist(strsplit(startpos[i],split=","))
  
  exonArray[line:(line+length(exonStartPos)-1),]<-
    cbind(name[i],name2[i],chr[i],strand[i],exonStartPos)
  
  line<-line+length(exonStartPos)
  
}

write.table(exonArray,
            "exonList.txt",
            append = FALSE, quote = FALSE, sep = ",",
            row.names = FALSE, col.names = TRUE)





seqList<-character(length=nrow(exonArray))



for(i in 1:length(seqList)){
  
  pos<-as.numeric(exonArray[i,5])
  
  seqList[i]<-getSeq(Hsapiens, exonArray[i,3],
                     start=pos-2,end=pos+2,
                     strand=exonArray[i,4],as.character=TRUE)
  
}				 


write.table(cbind(exonArray,seqList),
            "spliceSiteSeq.txt",
            append = FALSE, quote = FALSE, sep = ",",
            row.names = FALSE, col.names = TRUE)
