library(tidyverse)
library(ggplot2)

meta_file_path = "metadata.tsv"
data_file_path = "GSE45291_SLE.tsv"

df <- read.table(data_file_path, header = TRUE, sep = "\t")
meta <- read.table(meta_file_path, header = FALSE, sep = "\t")

# drop every columns in meta dataframe except V1 (Sample) and V5 (Condition)
meta <- meta[, c(1, 5)]

df1 <- as.data.frame(t(df))
df1 <- df1 %>% rownames_to_column(var = "gene")

# make dataframes less messy
colnames(df1) <- df1[1,]
colnames(meta) <- meta[1,]
df1 <- df1[-1,]
meta <- meta[-1,]

#merge 2 dataframe into 1
data <- merge(meta, df1, by.x = "Sample", by.y = 'gene', all = TRUE)

#convert cell values to numeric, exclude the first and second column
data[,3:ncol(data)] <- sapply(data[,3:ncol(data)], as.numeric)
rownames(data) <- data[,1]
data <- data[-1]

data_pca <- data[-1]

# run PCA on the data frame data
pca_result <- prcomp(data_pca)

# Step 2: Extract PCA Results
# Extract scores (transformed data) and loadings
pca_scores <- as.data.frame(pca_result$x)  # Principal component scores
pca_loadings <- as.data.frame(pca_result$rotation)  # Loadings

#assign color for conditions
colors <- as.list(ifelse(data[1]=='Healthy', 'red', 'blue'))

# Plot PC1 against PC2
pdf("PCAplot.pdf")
ggplot(pca_scores, aes(x = PC1, y = PC2, color=colors )) +
  geom_point() +
  labs(x = paste0("PC1 (", round(pca_result$sdev[1] * 100 / sum(pca_result$sdev), 2), "%)"),
       y = paste0("PC2 (", round(pca_result$sdev[2] * 100 / sum(pca_result$sdev), 2), "%)")) +
  ggtitle("PCA Plot")
dev.off()

#Run Differential Expression analysis
##BiocManager::install("DESeq2") #install DESeq2
##BiocManager::install("edgeR")
library(DESeq2)
library(edgeR)
library(ggplot2)
library(limma)
library(gplots)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(Glimma)
library(RColorBrewer)
dge <- DGEList(df,group=meta$Condition)
dds <- DESeqDataSetFromMatrix(countData=round(dge$counts), 
                              colData=meta, 
                              design=~Condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Condition", "Healthy", "SLE"))
res <- data.frame(res)

#label the top 20 significant genes
library(ggrepel)
res$deg_state = ifelse(res$pvalue<0.05, ifelse(res$log2FoldChange>0, 'under', 'over'), "none")
res$rank <- rank(res$pvalue)
res$gene_name <- df[1]
label<-c()
for (i in (1:length(df[,1]))){
  if (res$rank[i]<=20) {
    label <- c(label,df[i,1])
  }
  else{
    label <- c(label, NA)
  }
}
res$label<-label

#Draw volcano plot
pdf("volcano_plot.pdf", width=8, height=6)
p<-ggplot(res, aes(x=log2FoldChange, y=-log10(pvalue), col=deg_state, label=label)) + 
  geom_point() +
  ggtitle("Volcano plot of the edgeR results") +
  theme_minimal() +
  ylab("-log10(P)") +
  xlab("log2 fold change") +
  geom_vline(xintercept=0) +
  geom_text_repel(color="black", max.overlaps = 15) +
  scale_color_manual(values=c("none"="grey", "over"="red", "under"="blue"))
print(p)
dev.off()

#Segment RNA count-dataframe to collect 50 significantly expressed genes
library(dplyr)
sigs <- res[res$pvalue<0.05,]
name<-as.list(sigs$gene_name)
df_new<-data[,colnames(data) %in% name$gene]
df_new<-cbind(data[,1],df_new)
colnames(df_new)[1]="condition"
df_new<-data.frame(df_new)
df_new<-df_new[order(df_new$condition),]
kp<-df_new[,-1]

#Run Kmeans clustering for the DE-analyzed genes
library(cluster)
library(fpc)
clust <- kmeans(kp, centers = 4)
pdf("kmeans.pdf")
k<-clusplot(kp, clust$cluster, color=TRUE, shade=TRUE, 
            labels=0.5, lines=0)
print(k)
dev.off()

#Calculate the mean count per condition group
Healthy <- df_new[df_new$condition == "Healthy",]
geneH<-c()
valH<-c()
for(i in (2:51)){
  geneH<-c(geneH,colnames(Healthy)[i])
  valH<-c(valH,mean(Healthy[,i]))
}

SLE <- df_new[df_new$condition == "SLE",]
geneS<-c()
valS<-c()
for (i in (2:51)){
  geneS<-c(geneS,colnames(SLE)[i])
  valS<-c(valS,mean(SLE[,i]))
}

#assemble the mean value of each group into a dataset
hm<-rbind(valH,valS)
hm<-data.frame(hm)
colnames(hm)<-geneH
rownames(hm)<-c("Healthy","SLE")

#Draw boxplot for the gene with biggest difference
pdf("boxplot.pdf", width = 3, height = 5)
box<-boxplot(Healthy$RPS4Y1, SLE$RPS4Y1, col = "pink", ylim = c(0, 15))
print(box)
dev.off()

#plot heatmap to visualize the expression of significant DE genes in 2 conditions
library(ggplot2)
library(RColorBrewer)
library(gplots)
hm<-as.matrix(hm)
pdf("heatmap.pdf")
heat<-heatmap.2(hm, Rowv = NA, col = colorRampPalette(c("blue","red"))(25),
          key=TRUE,keysize=1.25,trace="none",
          cexRow = 1,
          cexCol = 0.5)
print(heat)
dev.off()

#Run correlation amongst significantly DE genes
correlation=cor(kp)

#Plot the heatmap for the correlation coefficient
pdf("correlation.pdf")
c<-heatmap.2(correlation, col = colorRampPalette(c("blue","red"))(25),
             key=TRUE,trace="none",cexRow = 0.5, cexCol = 0.5)
print(c)
dev.off()

#Build a dataframe that contains pairs of genes with strong correlation (r>0.7)
gene1 <- c()
gene2 <- c()
value <- c()
correlation <- data.frame(correlation)
for(i in (1:49)){
  for(j in ((i+1):50)){
    if (correlation[i,j]>0.7){
      gene1 = c(gene1, colnames(correlation[i,])[j])
      gene2 = c(gene2, colnames(correlation[j,])[i])
      value = c(value, correlation[i,j])
    }
  }
}
cor_important <- data.frame(gene1, gene2, value)
cor_sup_important <- cor_important[cor_important$value>0.95,]

#plot the strong linear correlation (>0.95)
for (i in (1:length(cor_sup_important[,1]))){
  gene1 = cor_sup_important[i,1]
  gene2 = cor_sup_important[i,2]
  dfn = kp[c(gene1, gene2)]
  lm<-lm(data=dfn)
  print(lm$coefficients[1])
  print(lm$coefficients[2])
  pdf(paste("plot", i, ".pdf", sep=""), width=4, height=4)
  p <- plot(dfn, pch = 16, cex = 1, col = "blue") + abline(lm$coefficients[1], lm$coefficients[2])
  print(p)
  dev.off()
}

