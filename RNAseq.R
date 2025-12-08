# Visualiation of clusters with heatmaps
# https://btep.ccr.cancer.gov/docs/data-visualization-with-r/Lesson5_intro_to_ggplot/#load-the-libraries
install.packages("pheatmap")
install.packages("tidyverse")
install.packages("ggplotify")
install.packages("heatmaply")
library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
library(ggplotify) ## to convert pheatmap to ggplot2
library(heatmaply) ## for constructing interactive heatmap
#setwd("C:/Users/raji/OneDrive - Yeshiva University/courses/bioinformatics_2022/Labs/Labs_2023/Lab11-MachineLearning/Clustering")
setwd("C:/Users/YU/OneDrive - Yeshiva University/courses_2023/bioinformatics_2022/Labs/Labs_2023/Lab11-MachineLearning/Clustering")
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