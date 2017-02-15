# change yo directory 

countdata <- read.table("~/Downloads/GSE85567_RNASeq_normalizedcounts.txt", row.names = 1, check.names = FALSE)
samples <- read.table("~/Downloads/GSE85566_processed_meth_data_covarates.txt", sep=",", header= TRUE, row.names=1)

library(DESeq2)
library("ggplot2")
library(tidyverse)

samplenames <- as.character(colnames(countdata))
samples <- samples %>% filter(ID %in% samplenames)

samplenames1 <- as.character(samples$ID)
countdata <- countdata[,samplenames1]


#putting things in order
samples$ID <- as.character(samples$ID)
samples <- samples[order(samples$ID),]
countdata <- countdata[,order(names(countdata))]

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = samples, design = ~ Status)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

data <- plotPCA(vsd, intgroup=c("ID"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=ID)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))


plo