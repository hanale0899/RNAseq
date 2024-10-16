library(limma)
library(edgeR)
library(gplots)

##inputdata
seqdata <- read.delim("Control.txt", stringsAsFactors = FALSE)

##Preprocess the raw counts
S <- colSums(seqdata)
stupidcpm <- cpm(seqdata) #normalize into cpm
head(stupidcpm)
thresh <- stupidcpm > 0.5 #keep records with cpm > 0.5
head(thresh)
rowSums(head(thresh))
keep <- rowSums(thresh) >= 5 #keep the genes with more than 5 records meet the cpm threshold
counts.keep <- seqdata[keep,] #create a new data frame with only kept genes
gene.list1 <- row.names(counts.keep) #export the gene list
write.table(gene.list1,file="controlgenes.txt",sep="\t")