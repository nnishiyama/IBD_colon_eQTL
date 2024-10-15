#!/usr/bin/R

# Normalize counts and inverse normal transform gene expression data

args <- commandArgs(TRUE)

library(data.table)
library(edgeR)
library(tidyverse)

# read in raw count data: m genes (rows) x annotations(col 1:6) + n samples (columns)
df <- fread(args[1], header=TRUE)
# read in sample list
coldata <- fread(args[2])
# find shared set of samples 
idx <- intersect(colnames(df), coldata$sample)
# keep only samples in list
df <- df %>% select(colnames(df)[1:6], all_of(idx))
coldata <- subset(coldata, sample %in% idx) %>% arrange(match(sample, idx))
# filter df rows with at least 5 reads per sample across 25% of samples
df <- df[apply(df[,-(1:6)], 1, function(x) length(which(x>=5))/length(x) >0.25),]
# subset count matrix
dat <- df[,-(1:6)]
# create group vector
group <- coldata$disease
# create DGElist
dge <- DGEList(counts=dat, group=group)
# calc TMM normalization factors
tmm <- calcNormFactors(dge,method='TMM')
# convert to TMM normalized counts
norm <- cpm(tmm)
# inverse normal transform expression data
int <- matrix(NA, nrow = nrow(norm), ncol = ncol(norm))
colnames(int) <- colnames(norm)
for (i in 1:nrow(norm)) {
  x = norm[i,]
  y = qnorm((rank(x, na.last='keep')-0.5)/sum(!is.na(x)))
  int[i,] = y
}
# add gene annotations back in
int.out <- cbind(df[,1:6],int)
# export INT data: m genes (rows) x annotations(col 1:6) + n samples (columns)
    write.table(int.out, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
