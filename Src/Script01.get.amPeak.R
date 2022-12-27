### Check for sangerseqR package.
if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")

if(!requireNamespace("sangerseqR", quietly = TRUE))
  BiocManager::install("sangerseqR")
  
if(!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")

if(!requireNamespace("this.path", quietly = TRUE))
  install.packages("this.path")
library(sangerseqR)

# Set working directory.
print(this.path::this.dir())


setwd(print(this.path::this.dir()))
setwd("../")

outpath="Amp/"

abi_files=list.files(path="ABI",pattern="*.ab1")


abi_files


for (num in 1:length(abi_files)){
	file=abi_files[num]
	seq <- readsangerseq(paste("ABI/",file,sep=""))
	seq <- makeBaseCalls(seq,ratio=0.2)
	peakAmp <- peakAmpMatrix(seq)
	peakAmp <- as.data.frame(peakAmp)
	colnames(peakAmp) <- c("A","C","G","T")
	peakAmp$ratio <- apply(peakAmp,1,function(x){a=sort(x,decreasing=T);a[2]/a[1]})
	peakAmp$sig <- ifelse(peakAmp$ratio>0.2,T,F)
	out_file=paste(outpath,file,".peakAmp.txt",sep="")
	write.table(peakAmp, file=out_file,row.names=F,quote=F,sep="\t")
}

