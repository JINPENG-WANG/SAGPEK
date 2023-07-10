# Class constructor.

setClass("sangerseq",
         representation(
           traceMatrix="matrix",
           peakPosMatrix="matrix",
           peakAmpMatrix="matrix"
         ))

setClass("abifHeader", 
         representation(
           abif="character",
           version="integer",
           name="character",
           number="integer",
           elementtype="integer",
           elementsize="integer",
           numelements="integer",
           dataoffset="integer",
           datahandle="integer"
         )
)

setClass("abifDirectory", 
         representation(
           name="character",
           tagnumber="integer",
           elementtype="integer",
           elementsize="integer",
           numelements="integer",
           datasize="integer",
           dataoffset="integer"
         )
)

setClass("abif", 
         representation(
           header="abifHeader",
           directory="abifDirectory",
           data="list"  #must be list because fields vary between file versions
         )
)

setGeneric("sangerseq", function(obj) standardGeneric("sangerseq"))


setGeneric("makeBaseCalls", 
           function(obj, ratio=.33) standardGeneric("makeBaseCalls"))


setGeneric("traceMatrix", 
           function(obj) standardGeneric("traceMatrix"))

setGeneric("traceMatrix<-", 
           function(obj, value) standardGeneric("traceMatrix<-"))

setGeneric("peakPosMatrix", 
           function(obj) standardGeneric("peakPosMatrix"))


setGeneric("peakPosMatrix<-", 
           function(obj, value) standardGeneric("peakPosMatrix<-"))

setGeneric("peakAmpMatrix", 
           function(obj) standardGeneric("peakAmpMatrix"))

setGeneric("peakAmpMatrix<-", 
           function(obj, value) standardGeneric("peakAmpMatrix<-"))



setMethod("sangerseq","abif",
          function(obj){
            res <- new("sangerseq")
            tracematrix <- matrix(c(obj@data$DATA.9,
                                    obj@data$DATA.10,
                                    obj@data$DATA.11,
                                    obj@data$DATA.12),
                                  ncol=4)
            orderedmatrix <- cbind(tracematrix[,regexpr("A",obj@data$FWO_.1)[1]],
                                   tracematrix[,regexpr("C",obj@data$FWO_.1)[1]],
                                   tracematrix[,regexpr("G",obj@data$FWO_.1)[1]],
                                   tracematrix[,regexpr("T",obj@data$FWO_.1)[1]]
            )
            
            basecalls1_old <- strsplit(obj@data$PBAS.2,"")[[1]]
            basecalls1_new <- basecalls1_old[basecalls1_old %in% DNA_ALPHABET]
            
            if (length(basecalls1_old) != length(basecalls1_new)) {
              #warning("Invalid characters removed from primary basecalls. This may result in basecalls being shifted. Please check chromatogram.")
            }
            basecalls1 <- paste0(basecalls1_new, collapse = "")
            
            # number of locations and calls do not always match
            basecallpositions1 <- obj@data$PLOC.2 + 1
            
            
            if(!is.null(obj@data$P2BA.1)) {
              basecalls2_old <- strsplit(obj@data$P2BA.1, "")[[1]]
              basecalls2_new <- basecalls2_old[basecalls2_old %in% DNA_ALPHABET]
              
              if (length(basecalls2_old) != length(basecalls2_new)) {
               # warning("Invalid characters removed from primary basecalls. This may result in basecalls being shifted. Please check chromatogram.")
              }
              basecalls2 <- paste0(basecalls2_new, collapse = "")
              
              # number of locations and calls do not always match
              #basecalls2 <- DNAString(substr(basecalls2,1,length(obj@data$PLOC.2)))
              basecallpositions2 <-obj@data$PLOC.2 + 1
            } else {
              #basecalls2 <- DNAString("")
              basecallpositions2 <- NA
            }
            
            if(!is.null(obj@data$P1AM.1)) {
              peakamps1 <- obj@data$P1AM.1
            } else {
              peakamps1 <- NA
            }
            if(!is.null(obj@data$P2AM.1)) {
              peakamps2 <- obj@data$P2AM.1
            } else {
              peakamps2 <- NA
            }
            
            res@traceMatrix <- orderedmatrix
            res@peakPosMatrix <- cbind(basecallpositions1, basecallpositions2, 
                                       NA, NA, deparse.level=0)
            res@peakAmpMatrix <- cbind(peakamps1, peakamps2, NA, NA, 
                                       deparse.level=0)
            return(res)
          }
)

setMethod("makeBaseCalls", "sangerseq",
          function(obj, ratio=.33) {
            #get peaks for each base
            Apeaks <- getpeaks(obj@traceMatrix[,1])
            Cpeaks <- getpeaks(obj@traceMatrix[,2])
            Gpeaks <- getpeaks(obj@traceMatrix[,3])
            Tpeaks <- getpeaks(obj@traceMatrix[,4])
            
            #get window around primary basecall peaks
            primarypeaks <- obj@peakPosMatrix[,1]
            diffs <- diff(c(0,primarypeaks))
            starts <- primarypeaks - 0.5*diffs
            stops <- c(primarypeaks[1:(length(primarypeaks)-1)] + 
                         0.5*diffs[2:length(diffs)], 
                       primarypeaks[length(diffs)] + 0.5*diffs[length(diffs)]
            ) 
            #hack for last peak. Just uses distance preceding peak 
            #as distance after peak
            
            #Now get max peak value for each channel in each peak window. 
            #If no peak return 0  
            primary <- NULL
            secondary <- NULL
            tempPosMatrix <- matrix(nrow=length(starts), ncol=4)
            tempAmpMatrix <- matrix(nrow=length(starts), ncol=4)
            
            ratio=0.33
            
            for(i in 1:length(starts)) {
              Apeak <- peakvalues(Apeaks, starts[i], stops[i])
              Cpeak <- peakvalues(Cpeaks, starts[i], stops[i])
              Gpeak <- peakvalues(Gpeaks, starts[i], stops[i])
              Tpeak <- peakvalues(Tpeaks, starts[i], stops[i])
              
              
              if(is.na(Apeak[2]) & 
                 is.na(Cpeak[2]) & 
                 is.na(Gpeak[2]) & 
                 is.na(Tpeak[2])) next #rare case where no peak found 
              signals <- c(Apeak[1], Cpeak[1], Gpeak[1], Tpeak[1])
              
              tempAmpMatrix[i,] <- signals
              positions <- c(Apeak[2], Cpeak[2], Gpeak[2], Tpeak[2])
              tempPosMatrix[i,] <- positions
              #print("signals")
              #print(signals)
              #print(positions)
              
              signalratios <- signals/max(signals, na.rm=TRUE)
              #print(signals)
              #print(max(signals,na.rm=TRUE))
              #print(signalratios)
              
              Bases <- c("A", "C", "G", "T")
              Bases[signalratios < ratio] <- NA
              #print(Bases)
              
              #sort by decreasing signal strength
              Bases <- Bases[order(signals, decreasing=TRUE)] 
              #print(Bases)
              positions <- positions[order(signals, decreasing=TRUE)]
              #print(positions)
              
            }
            
            tempPosMatrix[rowSums(!is.na(tempPosMatrix))>0,]
            tempAmpMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
            
            obj@peakPosMatrix <- tempPosMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
            obj@peakAmpMatrix <- tempAmpMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
            #obj@primarySeqID <- "sangerseq package primary basecalls"
            #obj@primarySeq <- DNAString(paste(primary, collapse=""))
            #obj@secondarySeqID <- "sangerseq package secondary basecalls"
            #obj@secondarySeq <- DNAString(paste(secondary, collapse=""))
            return(obj)
          }
)


setMethod("peakAmpMatrix", "sangerseq", function(obj) obj@peakAmpMatrix)


# Functions  used.


# Function-readBinaryData
readBinaryData <- function(filename) {
  fc <- file(filename, open = "rb")
  rawdata <- readBin(fc, what = "raw", n = 1.2*file.info(filename)$size)
  close(fc)
  return(rawdata)
}
###################################################



# Function-read.abif
read.abif <- function (filename) {
  
  #initialize object
  res <- new("abif")
  
  #get data
  rawdata <- readBinaryData(filename)
  
  res@header@abif <- RTC(rawdata[1:4])
  
  if (res@header@abif != "ABIF") 
    stop("file not in ABIF format")
  
  res@header@version <- SInt16(rawdata[5:6])
  res@header@name <- RTC(rawdata[7:10])
  res@header@number <- SInt32(rawdata[11:14])
  res@header@elementtype <- SInt16(rawdata[15:16])
  res@header@elementsize <- SInt16(rawdata[17:18])
  res@header@numelements <- SInt32(rawdata[19:22])
  res@header@dataoffset <- SInt32(rawdata[27:30])
  dataoffset <- res@header@dataoffset + 1
  res@header@datahandle <- SInt32(rawdata[31:34])
  #res@header@unused <- SInt16(rawdata[35:128], n = 47)
  #res@header@unused[1:length(res@header@unused)] <- 0
  
  #get directory
  for (i in seq_len(res@header@numelements)) {
    deb <- (i - 1) * res@header@elementsize + dataoffset
    direntry <- rawdata[deb:(deb + res@header@elementsize)]
    res@directory@name <- c(res@directory@name, RTC(direntry[1:4]))
    res@directory@tagnumber <- c(res@directory@tagnumber, 
                                 SInt32(direntry[5:8]))
    res@directory@elementtype <- c(res@directory@elementtype, 
                                   SInt16(direntry[9:10]))
    res@directory@elementsize <- c(res@directory@elementsize, 
                                   SInt16(direntry[11:12]))
    res@directory@numelements <- c(res@directory@numelements, 
                                   SInt32(direntry[13:16]))
    res@directory@datasize <- c(res@directory@datasize, 
                                SInt32(direntry[17:20]))
    res@directory@dataoffset <- c(res@directory@dataoffset, 
                                  SInt32(direntry[21:24]))
  }
  #fix for error in some .ab1 files that have the wrong data type for the 
  #PCON fields. Usually is 2 ("character") but should be 1 ("Uint8")
  res@directory@elementtype[res@directory@name == "PCON"] <- as.integer(1)
  
  #get data list
  res@data <- vector("list", length(res@directory@name))
  names(res@data) <- paste(res@directory@name, 
                           res@directory@tagnumber, 
                           sep = ".")
  for (i in seq_len(res@header@numelements)) {
    deb <- (i - 1) * res@header@elementsize + dataoffset
    if (res@directory@datasize[i] > 4) {
      debinraw <- res@directory@dataoffset[i] + 1
    }
    else {
      debinraw <- deb + 20
    }
    elementtype <- res@directory@elementtype[i]
    numelements <- res@directory@numelements[i]
    elementsize <- res@directory@elementsize[i]
    data <- rawdata[debinraw:(debinraw + numelements * elementsize)]
    if (elementtype == 1) 
      res@data[[i]] <- UInt8(data, n = numelements)
    if (elementtype == 2) 
      res@data[[i]] <- RTC(data)
    if (elementtype == 3) 
      res@data[[i]] <- UInt16(data, n = numelements)
    if (elementtype == 4) 
      res@data[[i]] <- SInt16(data, n = numelements)
    if (elementtype == 5) 
      res@data[[i]] <- SInt32(data, n = numelements)
    if (elementtype == 7) 
      res@data[[i]] <- f32(data, n = numelements)
    if (elementtype == 8) 
      res@data[[i]] <- f64(data, n = numelements)
    if (elementtype == 10) 
      res@data[[i]] <- list(year = SInt16(data, n = 1), 
                            month = UInt8(data[-(1:2)], n = 1), 
                            day = UInt8(data[-(1:3)], 
                                        n = 1)
      )
    if (elementtype == 11) 
      res@data[[i]] <- list(hour = UInt8(data, n = 1), 
                            minute = UInt8(data[-1], n = 1), 
                            second = UInt8(data[-(1:2)], n = 1), 
                            hsecond = UInt8(data[-(1:3)], n = 1)
      )
    if (elementtype == 18) {
      n <- SInt8(rawdata[debinraw])
      pString <- RTC(rawdata[(debinraw + 1):(debinraw + 
                                               n)])
      res@data[[i]] <- pString
    }
    if (elementtype == 19) 
      res@data[[i]] <- RTC(data[1:(length(data) - 1)])
    if (elementtype >= 1024) 
      res@data[[i]] <- data
    if (elementtype %in% c(12, 13)) 
      warning("unimplemented legacy type found in file")
    if (elementtype %in% c(6, 9, 14, 15, 16, 17, 20, 128, 
                           256, 384)) 
      warning("unsupported legacy type found in file")
  }
  return(res)
}


###Function-getpeaks.
getpeaks <- function(trace) {
  r <- rle(trace)
  indexes <- which(rep(diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2,
                       times = r$lengths))
  cbind(indexes, trace[indexes])
}


peakvalues <- function(x, pstart, pstop) {
  region <- x[x[,1] > pstart & x[,1] < pstop, ,drop=FALSE]
  if (length(region[,1]) == 0) return(c(0, NA))
  else return(c(max(region[,2], na.rm=TRUE), region[which.max(region[,2]),1]))
}

# Function-RTC
RTC <- function(x, multiple = TRUE, ...) {
  string <- suppressWarnings(rawToChar(x, multiple, ...))
  if(length(string) > 1) string <- paste(string, collapse="")
  #found that some ab1 files have unprinted characters at the end of the string
  #this is designed to remove them
  #string <- str_replace(string, "[^[:alnum:]]", "")
  return(string)
}
##########################################################

#Function-Seires Function
SInt32 <- function(f, n=length(f)/4) readBin(f, what = "integer", 
                                             signed = TRUE, endian = "big", 
                                             size = 4, n=n)
SInt16 <- function(f, n=length(f)/2) readBin(f, what = "integer", 
                                             signed = TRUE, endian = "big", 
                                             size = 2, n=n)
SInt8 <- function(f, n=length(f)) readBin(f, what = "integer", signed = TRUE, 
                                          endian = "big", size = 1, n=n)
UInt32 <- function(f, n=length(f)/4) readBin(f, what = "integer", 
                                             signed = FALSE, endian = "big", 
                                             size = 4, n=n)
UInt16 <- function(f, n=length(f)/2) readBin(f, what = "integer",
                                             signed = FALSE, endian = "big", 
                                             size = 2, n=n)
UInt8 <- function(f, n=length(f)) readBin(f, what = "integer", signed = FALSE, 
                                          endian = "big", size = 1, n=n)
f32 <- function(f, n=length(f)/4) readBin(f, what = "numeric", size = 4, n=n)
f64 <- function(f, n=length(f)/8) readBin(f, what = "numeric", size = 8, n=n)

#################################################


DNA_BASE_CODES <- c(A=1L, C=2L, G=4L, T=8L)
RNA_BASE_CODES <- c(A=1L, C=2L, G=4L, U=8L)


IUPAC_CODE_MAP <- c(
  A="A",
  C="C",
  G="G",
  T="T",
  M="AC",
  R="AG",
  W="AT",
  S="CG",
  Y="CT",
  K="GT",
  V="ACG",
  H="ACT",
  D="AGT",
  B="CGT",
  N="ACGT"
)

.IUPACcodes <- function(baseCodes)
{
  baseIsU <- names(baseCodes) == "U"
  if (any(baseIsU))
    names(baseCodes)[baseIsU] <- "T"
  code_list <- strsplit(IUPAC_CODE_MAP, "", fixed=TRUE)
  codes <- sapply(code_list, function(x) sum(baseCodes[x]))
  if (any(baseIsU))
    names(codes)[names(codes) == "T"] <- "U"
  codes
}

.DNAorRNAcodes <- function(baseCodes, baseOnly)
{
  #if (!isTRUEorFALSE(baseOnly))
  # stop("'baseOnly' must be TRUE or FALSE")
  codes <- .IUPACcodes(baseCodes)
  if (baseOnly) {
    codes[names(codes) %in% names(baseCodes)]
  } else {
    c(codes, `-`=16L, `+`=32L, `.`=64L)
  }
}

DNAcodes <- function(baseOnly) .DNAorRNAcodes(DNA_BASE_CODES, baseOnly)
RNAcodes <- function(baseOnly) .DNAorRNAcodes(RNA_BASE_CODES, baseOnly)

DNA_CODES <- DNAcodes(FALSE)
RNA_CODES <- RNAcodes(FALSE)

DNA_ALPHABET <- names(DNA_CODES)
RNA_ALPHABET <- names(RNA_CODES)

DNA_BASES <- names(DNAcodes(TRUE))
RNA_BASES <- names(RNAcodes(TRUE))


mergeIUPACLetters <- function(x)
{
  if (!is.character(x) || any(is.na(x)) || any(nchar(x) == 0))
    stop("'x' must be a vector of non-empty character strings")
  x <- CharacterList(strsplit(toupper(x), "", fixed=TRUE))
  yy <- unname(IUPAC_CODE_MAP[unlist(x, use.names=FALSE)])
  if (any(is.na(yy)))
    stop("some strings in 'x' contain non IUPAC letters")
  yy <- CharacterList(strsplit(yy, "", fixed=TRUE))
  y <- unstrsplit(sort(unique(IRanges:::regroupBySupergroup(yy, x))))
  names(IUPAC_CODE_MAP)[match(y, IUPAC_CODE_MAP)]
}


### Obtain arguments
args=commandArgs(T)
type=args[1]




#### Set working directories.

args=commandArgs()
scriptName=args[substr(args,1,7) == '--file=']

if (length(scriptName) == 0) {
  scriptName <- rstudioapi::getSourceEditorContext()$path
} else {
  scriptName <- substr(scriptName, 8, nchar(scriptName))
}

pathName = substr(
  scriptName, 
  1, 
  nchar(scriptName) - nchar(strsplit(scriptName, '.*[/|\\]')[[1]][2])
)


setwd(pathName)
setwd("../")
setwd("Amp")
allfile <- dir()
peakfile <- grep("*peakAmp.txt",allfile)
invisible(file.remove(allfile[peakfile]))
setwd("../")

outpath="Amp/"

abi_files=list.files(path="ABI",pattern="*.ab1")
if(type=="TEST"){
  abi_files=list.files(path="ABI/testdata",pattern="*.ab1")
}
message(" #SAGPEK: Step 1: Extracting signals from ABI-format files!")
N=length(abi_files)
width=options()$width

for (num in 1:length(abi_files)){
  file=abi_files[num]
  seq=""
  if(type=="TEST"){
    abifile <- read.abif(paste("ABI/testdata/",file,sep=""))
    seq <- sangerseq(abifile)
  }else{
    abifile <- read.abif(paste("ABI/",file,sep=""))
    seq <- sangerseq(abifile)
  }
  
  seq <- makeBaseCalls(seq,ratio=0.2)
  peakAmp <- peakAmpMatrix(seq)
  peakAmp <- as.data.frame(peakAmp)
  colnames(peakAmp) <- c("A","C","G","T")
  peakAmp$ratio <- apply(peakAmp,1,function(x){a=sort(x,decreasing=T);a[2]/a[1]})
  peakAmp$sig <- ifelse(peakAmp$ratio>0.2,T,F)
  out_file=paste(outpath,file,".peakAmp.txt",seq="")
  write.table(peakAmp,file=out_file,row.names=F,quote=F,sep="\t")
  
  # progress bar
  
  cat('[', paste0(rep('#', num/N*width), collapse=''),
      paste0(rep('-', width - num/N*width), collapse=''),
      ']',
      round(num/N*100),'%')
  Sys.sleep(0.05)
  if(num==N) cat ('   DONE!\n')
  else cat('\r')
}
  
    
  













