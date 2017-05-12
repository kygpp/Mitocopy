#    create MT bins for QDNAseq analysis (Adapted from QDNAseq package)
print("-------------xiao> load library")
library(QDNAseq)
library(BiocGenerics)
library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg19)
library(QDNAseq.hg19)
#define binsize
binSize<-1000
print("-------------> load library done")


bsgenome<-BSgenome.Hsapiens.UCSC.hg19

chrs <- GenomeInfoDb::seqnames(bsgenome)
        info <- GenomeInfoDb::genomeStyles(GenomeInfoDb::organism(bsgenome))
        style <- GenomeInfoDb::seqlevelsStyle(bsgenome)[1]
        chrs <- info[, style]
       selectedMT <- grep("^(chr)?M(T)?$", chrs)
      chrs <- chrs[selectedMT]
lengths <- GenomeInfoDb::seqlengths(bsgenome)[chrs]
binWidth <- as.integer(binSize * 1L)
chrdatafunction<-function(chr) {
        chr.size <- lengths
        chr.starts <- seq(from=1L, to=chr.size, by=binWidth)
        chr.ends <- chr.starts + binWidth - 1L
        chr.ends[length(chr.ends)] <- chr.size
        chr.seq <- BSgenome::getSeq(bsgenome, chr, as.character=TRUE)
        bin.seq <- substring(chr.seq, first=chr.starts, last=chr.ends)
        acgt <- gsub("[^ACGT]", "", bin.seq)
        at<-gsub("[^AT]", "", acgt)  ####2GC, AT along with CT, AG, AC and GT
        ag<-gsub("[^AG]", "", acgt)
        gt<-gsub("[^GT]", "", acgt)
        ac<-gsub("[^AC]", "", acgt)
        ct<-gsub("[^CT]", "", acgt)
        cg <- gsub("[^CG]", "", acgt)
        chr.bases <- nchar(acgt) / (binWidth) * 100
        chr.gc <- nchar(cg) / nchar(acgt) * 100
       list(start=chr.starts, end=chr.ends, bases=chr.bases, gc=chr.gc)
    }
chrData<-chrdatafunction(chrs)
chromosome<-rep(chrs, ceiling(lengths/binWidth))
bins <- data.frame(chromosome, chrData, stringsAsFactors=FALSE)
colnames(bins)<-c("chromosome","start","end","bases","gc")
bins$chromosome <- sub("^chr", "", bins$chromosome)
rownames(bins)<-sprintf("%s:%i-%i", bins$chromosome, bins$start, bins$end)
bins$mappability <- calculateMappability(bins, bigWigFile='wgEncodeCrgMapabilityAlign50mer.bigWig', bigWigAverageOverBed='bigWigAverageOverBed')
bins$use <- bins$bases > 0
#saveRDS(bins,"hg19_MTcreated_binsize1000bp.rds")
#readRDS(bins,"hg19_MTcreated_binsize1000bp.rds")
##Use custom control files for calculating residuals
control<-c(c1.bam,c2.bam,c3.bam,c4.bam)
bins$residual <- iterateResiduals(binReadCounts(bins,bamfiles=control))
#For now, there's no blacklist available for MT genome.
saveRDS(bins,"hg19_MTcreated_binsize1000bp.rds")
#now,you have binsize of 1 kb for MT genome.




 

