#########################
args <- commandArgs(TRUE)
INPUTDir<- args[1]	
OutDir <- args[2]	
bin_size <- as.numeric(args[3])
#library
print("-------------load library")
library(QDNAseq)
library(BiocGenerics)
library(Rsamtools)
library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg19)
library(QDNAseq.hg19)
library(DNAcopy)
library(CGHregions)     

chrs <- c(1:22,"X","Y")
bins <- getBinAnnotations(binSize=bin_size)
print(" loading done ! < -----------")

##################################################################
#STEP1: calculate Mitochondrial average coverage without correction
# input files (bam)
bamFile <- list.files(bamDir,"\\.bam$") 
bamID <- head(unlist(strsplit(tail(unlist(strsplit(bamFile, "/")),1),"[.]")),1)
INPUTBAM<-paste(bamDir,bamFile,sep="/")
loop_samples(bins,INPUTBAM,bin_size,bamID)
avecovFile<-list.files(INPUTDir,"*avgCov.txt$")
INPUTcov<-paste(INPUTDir,avecovFile,sep="/")
MTcov<-MTcovcal(INPUTcov)
#STEP2: investigate segmental aneuploidy
#library(QDNAseq)
#library(CGHbase)
#input files (copynumbersegmentedfiles)

copynumbersegFile <- list.files(INPUTDir,"\\.rds$")
segID <- head(unlist(strsplit(tail(unlist(strsplit(copynumberseg, "/")),1),"[.]")),1)
INPUTcopynumberseg<-paste(INPUTDir,copynumberseg,sep="/")
process(INPUTcopynumberseg)

#sTEP3:calculate correctionfactor

tableFile <- list.files(OutDir,"\\.txt$")
tableID <- head(unlist(strsplit(tail(unlist(strsplit(tableFile, "/")),1),"[.]")),1)
INPUTtable<-paste(OutDir,tableFile,sep="/")

F<-correctionfactor(INPUTtable,method="log2ratio")

#STEP4: calculate Mitochondrial average coverage after correction
MTcovcorrected<-MTcov*F

resultcov<-data.frame(MTcov,F,MTcovcorrected)
OUTFILE<-paste0(OutDir,"/",tableID,"_","covFc",".txt")
write.table(resultcov,OUTFILE)

######################################
################### Define functions
#bamcoverage is to calculate averagecoverage for each chromosome.
bamcoverage <- function (bamfile) {
	  # read in the bam file
	  bam <- scanBam(bamfile)[[1]] # the result comes in nested lists
	  # filter reads without match position
	  ind <- ! is.na(bam$pos)
	  ## remove non-matches, they are not relevant to us
	  bam <- lapply(bam, function(x) x[ind])
	  ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
	  ## names of the bam data frame:
	  ## "qname"  "flag"   "rname"  "strand" "pos"    "qwidth"
	  ## "mapq"   "cigar"  "mrnm"   "mpos"   "isize"  "seq"    "qual"
	  ## construc: genomic ranges object containing all reads
	  ranges <- GRanges(seqnames=Rle(bam$rname), ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname )
	  ## returns a coverage for each reference sequence (aka. chromosome) in the bam file
	  return (mean(coverage(ranges)))      
}  
#loog_samples is to run copy number analysis genome-wide 
loop_samples <- function(bins,InputFile,BIN_S,bamID) {

print(paste0("============= InputFile:", InputFile))
	OutFile <- paste0(OutDir,"/",bamID,"_","avgCov",".txt")
	write.table(bamcoverage(InputFile),OutFile,sep="\t",col.names=F)
print(paste0("============= binReadCounts"))
	readCounts <- binReadCounts(bins, bamfiles=InputFile)
	#write.table(readCounts,"~/tmp1.csv",sep="\t",col.names=F)
print(paste0("============= applyFilters"))
	readCounts <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)
	#write.table(readCounts,"~/tmp2.csv",sep="\t",col.names=F)
print(paste0("============= estimateCorrection"))
	readCounts <- estimateCorrection(readCounts)
	#write.table(readCounts,"~/tmp3.csv",sep="\t",col.names=F)

print(paste0("============= plot noise"))
	OutFile <- paste0(OutDir,"/",bamID,"_","noisePlot",BIN_S,".pdf")
	pdf(OutFile)
	noisePlot(readCounts)
	dev.off()
print(paste0("============= applyFilters correctBins normaizeBins smoothOutlierBins"))
	readCounts <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE, chromosomes=NA)	# includes both X & Y
	copyNumbers <- correctBins(readCounts)
	copyNumbersNormalized <- normalizeBins(copyNumbers)
	copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
print(paste0("============= segmentBins & normalizeSegmentedBins"))
	copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
	copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

print(paste0("============= exportBins"))
     	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"_seg.bed")
     	exportBins(copyNumbersSegmented,file = OutFile, format="bed", filter=F,type=c("segments"))

print(paste0("============= output copyNumbersSegmented into rds file"))
    	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"before.rds")
	saveRDS(copyNumbersSegmented, file=OutFile)
##copyNumbersSegmented <- readRDS("")


print(paste0("============= plot genome-wide ============="))
	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,".pdf")
	print(OutFile)
	pdf(OutFile,paper='a4r')
	#pdf(OutFile)
	plot(copyNumbersSegmented, plot.type="w",main=paste0("bin=",BIN_S))	# paste0(unlist(strsplit(bamID, "-"))[2],",bin=",BIN_S))
	dev.off()


print(paste0("============= plot per_chr ============="))
	f.data <- as(featureData(copyNumbersSegmented), "data.frame")
	OutFile <- paste0(OutDir,"/",bamID,"_","byChr",BIN_S,".pdf")
	pdf(OutFile)
#	for (chromosome in c(1:22,"X","Y")){
	for (chromosome in unique(f.data$chromosome)) {
		select <- which(f.data$chromosome == chromosome)
		print(paste("Plotting chromosome:", chromosome, sep = " "))
		## Looping over samples, make a plot per sample, per chromosome
		for (c in 1:ncol(copyNumbersSegmented)){
			sample.name <- sampleNames(phenoData(copyNumbersSegmented))[c] # You can use this, for plot title or plot name
			sample.name <- bamID
			plot(copyNumbersSegmented[select, c], ylim=c(-2,2),main=paste0("Chr",chromosome,",bin=",BIN_S))         
		}
	}
	dev.off()

   
print(paste0("============= set cutoff & plot ============="))
	copyNumbersCalled <- try(callBins(copyNumbersSegmented),silent=T)
	if (inherits(copyNumbersCalled,'try-error')) {
		copyNumbersCalled <- try(callBins(copyNumbersSegmented,method="cutoff",cutoffs=c(-0.1, 0.1)),silent=T)
	}  	

	cgh <- makeCgh(copyNumbersCalled) 
	regions <- CGHregions(cgh)
	OutFile <- paste0(OutDir,"/",bamID,"_",BIN_S,"_after_regions.bed")
	exportBins(regions,file = OutFile, format="bed", filter=T,type=c("calls"))

	OutFile <- paste0(OutDir,"/",bamID,"_","Calls",BIN_S,".pdf")
	pdf(OutFile)
	plot(copyNumbersCalled)
	dev.off()

}
##MTcalculate is to calculate MTcoverage compared to autosomecoverage, input is avecov.txt caculate from bamcoverage()
#output is a single numeric vector
MTcovcal<-function(object){
  data<-read.csv(object,sep="\t",header = FALSE)
  MTaveragecoverage<-awk$V2[agrep("MT$", data$V1)]
  autosomeaeragecoverage<-data$V2[1:22]
  MTcopynumber<-mean(MTaveragecoverage/autosomeaeragecoverage*2)
  MTcopynumber
}

## segmentation function is to calculate the size of each segment, process is to excute segmentation function with 
#preprocessing input (copynumbersegmented.rds) object file. The output is a dataframe (loadpackage QDNAseq,CGHbase)
process<-function(object){
  object<-readRDS(object)
  chromosome <- fData(object)$chromosome
  start <- fData(object)$start
  end <- fData(object)$end
  log2ratio<-log2(assayDataElement(object, "segmented")[,1])
  copyNumbersCalled <- try(callBins(object),silent=T)
  if (inherits(copyNumbersCalled,'try-error')) {
    copyNumbersCalled <- try(callBins(object,method="cutoff",cutoffs=c(-0.1, 0.1)),silent=T)
  }  	
  loss<-probloss(copyNumbersCalled)[,1]
  gain<-probgain(copyNumbersCalled)[,1]
  data<-data.frame(chromosome,start,end,log2ratio,loss,gain)
  ind <- sapply(data$log2ratio, function(x) all(is.na(x)))
  data<- data[ !ind, ]
  data$chr<-factor(data$chr,levels = c(seq(1,22),"X","Y"))
  OutFile <- paste0(OutDir,"/",bedID,"_",".txt")
  write.table(segmentation(data),OutFile)
  segmentation(data)
}

segmentation<-function(data){
  segment<-NULL 
  for (k in c(seq(1,22),"X","Y")) {
    datachr<-data[data$chromosome==k,]
    runs<-rle(datachr$log2ratio)
    chromosome<-rep(k,length(runs$values))
    segstart<-NULL
    segend<-NULL
    probloss<-NULL
    probgain<-NULL
    for (i in c(1:length(runs$values))) {
      runs.lengths.cumsum = cumsum(runs$lengths)
      ends = runs.lengths.cumsum[i]
      newindex = ifelse(i>1, i-1, 0)
      starts= runs.lengths.cumsum[newindex] + 1
      if (0 %in% newindex){
        starts = c(1)
      }
      segstart[i]<-datachr$start[starts]
      segend[i]<-datachr$end[ends]
      probloss[i]<-datachr$loss[starts]
      probgain[i]<-datachr$gain[starts]
    }
    log2ratio<-runs$values
    runlength<-runs$lengths
    seglength<-segend-segstart
    df<-data.frame(chromosome,segstart,segend,seglength,log2ratio,probloss,probgain)
    segment[[k]]<-df
  }
  segment<- do.call("rbind", segment)
  rownames(segment)<-NULL
  return(segment)
}
################# correctionfactor function is to calculate the correction factor, based on genome size of each sample. This
#function allow custom-defined log2ratio threshold, or probability threshold. Either method has to be specified when using the
#function, default value for log2ratio is log2(c(loss=0.7,gain=1.6)),default value for ploss=0.9,pgain=0.9.
correctionfactor<-function (object,method=c("log2ratio","probability"),cutoff=log2(c(loss=0.7,gain=1.6)),ploss=0.9,pgain=0.9) {
  method <- match.arg(method)
  table<-read.csv2(object,header=TRUE,sep = "",stringsAsFactors = FALSE)
  if (method == "log2ratio")  {
    table$log2ratio<-as.numeric(table$log2ratio)
    losslength<-sum(as.numeric(table[table$log2ratio<as.numeric(cutoff["loss"]),"seglength"]))
    gainlength<-sum(as.numeric(table[table$log2ratio>as.numeric(cutoff["gain"]),"seglength"]))
  }
  if (method == "probability")  {
    losslength<-sum(table[table$probloss>ploss,"seglength"])
    gainlength<-sum(table[table$probgain>pgain,"seglength"])
  }
  REF<-6072607692 
 REF<-6072607692 
  F<-(REF-losslength+gainlength)/REF #QDNAseq automatically detected Y chromosome as one copy loss.
}


print("-------------define function done")
