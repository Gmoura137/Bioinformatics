# Visualiation of clusters with heatmaps
install.packages("pheatmap")
install.packages("tidyverse")
install.packages("ggplotify")
install.packages("heatmaply")
library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
library(ggplotify) ## to convert pheatmap to ggplot2
library(heatmaply) ## for constructing interactive heatmap
setwd("C:/Users/YU/OneDrive - Yeshiva University/courses_2023/bioinformatics_2022/Labs/Labs_2023/Lab12-13-MachineLearning/Clustering")
mat<-read.csv("airway_scaledcounts.csv",header=TRUE,row.names=1,
              check.names=FALSE)
head(mat)
metadata<-read.csv("airway_metadata.csv",header=TRUE,row.names=1,
              check.names=FALSE)
head(metadata)
# Remove all the NA rows from mat and assign it to a new dataframe.
mat_na <- na.omit(mat)
head(mat_na)
mat_na$controlmean <- (mat_na$SRR1039508+mat_na$SRR1039512+mat_na$SRR1039516+mat_na$SRR1039520)/4
mat_na$treatedmean <- (mat_na$SRR1039509+mat_na$SRR1039513+mat_na$SRR1039517+mat_na$SRR1039521)/4
head(mat_na)
mat_zero <- mat_na[apply(mat_na, 1, function(row) all(row !=0 )), ]  # Remove zero-rows
head(mat_zero)                                                         # Print updated datafilter(mat_na, controlmean > 0 & treatedmean > 0)
mat_zero$normalized <- log2(mat_zero$treatedmean/mat_zero$controlmean)
head(mat_zero)
ordered_data <- mat_zero[order(-mat_zero$normalized),]
head(ordered_data)
top_20 <- head(ordered_data,20)
top_20
rownames(top_20)
pheatmap(top_20, scale='column')
pheatmap(top_20,scale="row",
         cutree_rows=2,cutree_cols=2,main="Expression and Clustering of Top 20 DE genes",color=colorRampPalette(c("navy", "white", "red"))(50))
jpeg("heatmap_de.jpg", res=300, width=7, height=4.5, unit="in")
pheatmap(top_20,scale="row",color=colorRampPalette(c("navy", "white", "red"))(50),
         cutree_cols=2, cutree_rows=2,
         main="Expression and clustering of top DE genes",
         fontsize=11, cellwidth=35, cellheight=10.25)
rna_seq<-read.csv("RNAseq_mat_top20.csv",header=TRUE,row.names=1,
                  check.names=FALSE)
head(rna_seq)
pheatmap(rna_seq, scale="row")
# color customization
pheatmap(rna_seq,scale="row",
         color=colorRampPalette(c("navy", "white", "red"))(50))
#create data frame for annotations
library(dplyr)
library(tibble)
df <- data.frame(sample=as.character(colnames(rna_seq)), dex="Treatment") %>% 
  column_to_rownames("sample")
head(df)
df$dex<-ifelse(rownames(df) %in% c("508","512","516","520"),
               "untrt","trt")
df
pheatmap(rna_seq,scale="row", annotation_col = df,
         annotation_colors=list(dex=c(trt="orange",untrt="black")),
         color=colorRampPalette(c("navy", "white", "red"))(50))
pheatmap(rna_seq,scale="row", annotation_col = df,
         annotation_colors =list(dex=c(trt="orange",untrt="black")),
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cutree_cols=2, cutree_rows=2)
# add a title; adjust the font size etc.
pheatmap(rna_seq,scale="row", annotation_col = df,
         annotation_colors=list(dex=c(trt="orange",untrt="black")),
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cutree_cols=2, cutree_rows=2, 
         main="Expression and clustering of top 20 DE genes",
         fontsize=11, cellwidth=35, cellheight=10.25)
# save a jpeg file for printing
jpeg("rna_seq.jpg", res=300, width=7, height=4.5, unit="in")
pheatmap(rna_seq,scale="row", annotation_col = df,
         annotation_colors=list(dex=c(trt="orange",untrt="black")),
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cutree_cols=2, cutree_rows=2,
         main="Expression and clustering of top DE genes",
         fontsize=11, cellwidth=35, cellheight=10.25)
dev.off()
dev.off()
