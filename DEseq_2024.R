#Load the Libraries
#install.packages("DESeq2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)
install.packages("ggrepel")
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)


#set working directory
setwd("C:/Users/YU/OneDrive - Yeshiva University/courses_2023/bioinformatics_2022/Labs/Labs_2023/Lab12-13-MachineLearning/Clustering")
#Load the count data
mat<-read.csv("airway_scaledcounts.csv",header=TRUE,row.names=1,
              check.names=FALSE)
colnames(mat)
head(mat)
#Load the sample information
metadata<-read.csv("airway_metadata.csv",header=TRUE,row.names=1,
                   check.names=FALSE)
colnames(metadata)
head(metadata)


#set factor levels
metadata$dex <- factor(metadata$dex)



#Create a DESeq object and import the count data and sample information
dds <- DESeqDataSetFromMatrix(countData = mat, colData = metadata, design = ~dex)
#Set the reference for the Treatment factor
dds$dex <- factor(dds$dex, levels = c("control","treated"))
#
# Filter the genes
keep <- rowSums(counts(dds)) >=5
dds <- dds[keep,]
dds
#
# Perform the statistical test to identify differentially expressed genes
dds <- DESeq(dds)
deseq_result <- results(dds)
deseq_result
#
# Change DESeq results to R dataframe
deseq_result <- as.data.frame(deseq_result)
class(deseq_result)
head(deseq_result)


#Order the result table by increasing p-value
colnames(deseq_result)
deseq_result_ordered <- deseq_result[order(deseq_result$pvalue),]
head(deseq_result_ordered)


#Make some queries
deseq_result["ENSG00000152583",]


#Extract the most differentially expressed genes log2fold <-1 or >1
#Step1 Filter based on p-value

filtered <- deseq_result %>% filter (deseq_result$padj <0.05)
#Step 2 -- Filter based on Fold change
filtered <- filtered %>% filter(abs(filtered$log2FoldChange) >1)

dim(deseq_result)
dim(filtered)
write.csv(deseq_result,'de_result.all.csv')
write.csv(filtered,'de_result.filtered.csv')


#Save the normalized count
normalized_counts <- counts(dds,normalized=TRUE)
head(normalized_counts)
write.csv(normalized_counts,'normalized_counts.csv')
#Visualization
#Dispersion Plot
plotDispEsts(dds)
#PCA 
#Dimensionality Reduction Technique
#variance stabilizing transformation
vsd <- vst(dds,blind=FALSE)
plotPCA(vsd,intgroup=c("dex"))
#HeatMaps
#We will use the top 20 Genes
#sample-to-sample distance (with clustering)
#Generate the distance matrix
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix)
#set color scheme
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
#generate the heatmap
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, col=colors)
#heatmap of log2normalized counts of top 20 genes
top_20_hits <- deseq_result[order(deseq_result$padj), ][1:20,]
top_20_hits <- row.names(top_20_hits)
top_20_hits
#
rld <- rlog(dds,blind=FALSE)
pheatmap(assay(rld)[top_20_hits,], cluster_rows=FALSE,show_rownames=TRUE,cluster_cols=FALSE)
#with clustering
pheatmap(assay(rld)[top_20_hits,])
#add annotation
annot_info <- as.data.frame(colData(dds)[,c('dex','sizeFactor')])
colData(dds)
annot_info
pheatmap(assay(rld)[top_20_hits,],cluster_rows=FALSE,show_rownames=TRUE,cluster_cols=FALSE,
          annotation_col=annot_info)
#heatmap of z-scores

cal_z_score <- function(x) {(x-mean(x))/sd(x)}
zscore_all <- t(apply(normalized_counts, 1, cal_z_score))
zscore_subset <- zscore_all[top_20_hits,]
pheatmap(zscore_subset)
