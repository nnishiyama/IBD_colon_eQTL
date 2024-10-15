#!/usr/bin/R

# Calculate RUVseq factors

args <- commandArgs(TRUE)

library(data.table)
library(tidyverse)
library(DESeq2)
library(RUVSeq)
library(vsn)

# load raw counts
cts <- fread(args[1], header=TRUE)
# load meta data
coldata <- fread(args[2])
# arrange sample IDs to match
idx <- intersect(colnames(cts), coldata$sample)
coldata <- subset(coldata, sample %in% idx) %>% arrange(match(sample, idx))
# For QTLtools BED-formatted files:
cts <- cts %>% select(colnames(cts)[1:6], all_of(idx))
# THIS MUST RETURN TRUE!
all(coldata$ID == colnames(cts)[-(1:6)])
# Filter out regions with low coverage
cts <- cts[apply(cts[,-(1:6)], 1, function(x) length(which(x>=5))/length(x) >0.25),]
# Pull out gene info
genes <- cts[,1:6]
# Subset just count matrix and convert to integers
cts <- round(cts[,-(1:6)])
# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~1) # Use 1 for intercept only model (no contrast between groups)
# Normalize based on sequencing depth
vsd <- vst(dds, blind = TRUE)
meanSdPlot(assay(vsd))
plotPCA(vsd,intgroup='batch') + theme(plot.title = element_text(hjust = 0.5)) + ggtitle('PCA on uncorrected gene counts') + labs(color = 'Batch')
plotPCA(vsd,intgroup='disease') + theme(plot.title = element_text(hjust = 0.5)) + ggtitle('PCA on uncorrected gene counts') + labs(color = 'Disease')
plotPCA(vsd, intgroup='sex') + theme(plot.title = element_text(hjust = 0.5)) + ggtitle('PCA on uncorrected gene counts') + labs(color = 'Sex')
# Create model matrix to run limma removeBatchEffect
mm <- model.matrix(~1, colData(vsd))
sex <- as.numeric(as.factor(coldata$sex))
batch <- as.numeric(as.factor(coldata$batch))
pc1 <- coldata$genotypePC1
pc2 <- coldata$genotypePC2
pc3 <- coldata$genotypePC3
pc4 <- coldata$genotypePC4
# Remove batch effect
assay(vsd) <- removeBatchEffect(assay(vsd), batch=batch, covariates=cbind(sex, pc1, pc2, pc3, pc4), design=mm)
# corrected PCAs
plotPCA(vsd,'batch') + theme(plot.title = element_text(hjust = 0.5)) + ggtitle('PCA removing batch & sex effects') + labs(color = 'Batch')
plotPCA(vsd,"disease") + theme(plot.title = element_text(hjust = 0.5)) + ggtitle('PCA removing batch & sex effects') + labs(color = 'Disease')
plotPCA(vsd, 'sex') + theme(plot.title = element_text(hjust = 0.5)) + ggtitle('PCA removing batch & sex effects') + labs(color = 'Sex')
################################################################################
# RUVseq
# Pull out the normalized matrix
nc <- as.data.frame(assay(vsd))
row.names(nc) <- genes$gene
# Calculate the average gene expression of each genes and take the top 5000 highest expressed
nc$avg_exp <- rowMeans(nc,na.rm=TRUE)
nc <- nc %>% arrange(desc(avg_exp))
nc <- nc[1:5000,]
nc <- nc %>% select(-avg_exp)
# Calculate the variance of each genes, and choose the lowest 1000 genes as the negative control gene
nc$cv <- rowSds(as.matrix(nc))/rowMeans(nc,na.rm=TRUE)
nc <- nc %>% arrange(cv)
nc <- nc[1:1000,]
nc <- nc %>% select(-cv)
# Create a new Expression Set and perform an upper quartile normalization
nc <- round(nc*1000)
pre <- newSeqExpressionSet(as.matrix(nc),phenoData = data.frame(coldata, row.names = coldata$sample))
pre <- betweenLaneNormalization(pre, which="upper")
plotRLE(pre, outline=FALSE, ylim=c(-.1, .1), col=as.numeric(as.factor(coldata$disease)))
plotPCA(pre, col=as.numeric(as.factor(coldata$disease)), k=2, cex=1.2)
# set k
k <- round(0.25*length(idx))
# Run RUVseq
post <- RUVg(pre, rownames(nc), k=k) 
# Create an RLE plot for review
plotRLE(post, outline=FALSE, ylim=c(-.1, .1), col=as.numeric(as.factor(coldata$disease)))
# Creates a PCA plot of the negative control genes. Should not see much of a pattern.
plotPCA(post, col=as.numeric(as.factor(coldata$disease)), k=2, cex=1.2)
# Pull out the metadata and RUV factors
coldata_new<-pData(post)
ruv_factors <- coldata_new[,-(1:ncol(coldata))]
################################################################################
# Save RUV factors to file
# need to transpose matrix before saving
ruv_factors <- t(ruv_factors)
write.table(ruv_factors,  sep = '\t', quote = FALSE, row.names = FALSE,
            file = 'data/UNC/RUV_factors.txt')
################################################################################
# Set up new DESeq object to visualize RUV correction
dds1 <- DESeqDataSetFromMatrix(countData=cts,colData=coldata_new,design = ~1)
vsd <- vst(dds1)
# Remove batch effects
assay(vsd) <- removeBatchEffect(assay(vsd), batch = batch, covariate = cbind(sex,ruv_factors[,1:19]),design=mm)
# Corrected PCA with RUV factors
plotPCA(vsd,"disease") + theme(plot.title = element_text(hjust = 0.5)) + ggtitle('PCA removing batch, sex, RUV factors') + labs(color = 'Disease')
