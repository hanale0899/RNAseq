library(Rsubread)
##getfiles
rseq.files <- list.files(path = './data', pattern = ".fq.gz$", full.names = TRUE)
rseq.files

##build_index
buildindex(basename="GRCh38",reference="GRCh38_latest_genomic.fq.gz") #Download the Index from NCBI database

##align
A<-align(index="data/GRCh38",
      readfile1=rseq.files,
      unique = FALSE,
      nBestLocations = 6,
      output_format="BAM",
      output_file=out.multi.bam)
args(A)

##Counting
bam.files <- list.files(path = "./", pattern = ".bam$", full.names = TRUE)
fc <- featureCounts(bam.files, annot.ext="GCF_000001405.40_GRCh38.p14_genomic.gtf",isGTFAnnotationFile=TRUE) #Download gtf file from NCBI database
rawcount <- fc$counts
colnames(rawcount)<-c("C01","C02","C03","C04","C05","C06","C07","C08","C09","C10","T01","T02","T03","T04","T05","T06","T07","T08","T09","T10")
write.table(rawcount,file="Control.txt",sep='\t')

