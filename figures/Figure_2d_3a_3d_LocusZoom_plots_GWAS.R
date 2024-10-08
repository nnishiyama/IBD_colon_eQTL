#Figures 2d, 3a, 3d. LocusZoom plots: ELMO1, ABO, TNFRSF14
library(locuszoomr)
library(LDlinkR)
library(AnnotationHub)
library(tidyverse)
library(data.table)
library(memoise)

# load coloc SNPs
load('data/UNC/RData/locusZoom_matrices_GWAS-eQTL_coloc_H4.RData')
# reformat
df.barcuva <- df.barcuva %>% unite('id', chr:pos, sep = ':', remove = FALSE)
df.gtex <- df.gtex %>% unite('id', chr:pos, sep = ':', remove = FALSE)
df.gwas <- df.gwas %>% unite('id', chr:pos, sep = ':', remove = FALSE)
df.ibd <- df.ibd %>% unite('id', chr:pos, sep = ':', remove = FALSE)
df.barcuva$chr <- gsub('chr', '', df.barcuva$chr)
df.gtex$chr <- gsub('chr', '', df.gtex$chr)
df.gwas$chr <- gsub('chr', '', df.gwas$chr)
df.ibd$chr <- gsub('chr', '', df.ibd$chr)
################################################################################
# set up for LocusZoom plots
# LDlink token
t <- '' # requires personal access token
# load Ensembl database
ah <- AnnotationHub()
ensembl.v105 <- ah[['AH98047']]
# create functions to calculate LD
mem_LDmatrix <- memoise(LDlinkR::LDmatrix)
get_LD <- function(l, coord = 'id', index_snp = '', pop = 'EUR', genome_build = 'grch38', token = '', ...){
  labs <- l$labs
  index_snp <- index_snp
  r <- l$data[which(l$data[, labs] == index_snp), coord]
  rslist <- l$data[, coord]
  int <- 1000 #1000
  if (length(rslist > int)) {
    rslist <- rslist[order(l$data$logP, decreasing = TRUE)[seq_len(int)]]
  }
  rslist <- unique(c(rslist, r))
  ldm <- mem_LDmatrix(rslist, pop = pop, genome_build = genome_build, token = token, ...)
  ld <- ldm[, index_snp]
  l$data$ld <- ld[match(l$data[, labs], ldm$RS_number)]
  l
}
get_LD_small <- function(l, coord = 'id', index_snp = '', pop = 'EUR', genome_build = 'grch38', token = '', ...){
  labs <- l$labs
  index_snp <- index_snp
  r <- l$data[which(l$data[, labs] == index_snp), coord]
  rslist <- l$data[, coord]
  # no filtering
  # int <- 1000 #1000
  # if (length(rslist > int)) {
  #   rslist <- rslist[order(l$data$logP, decreasing = TRUE)[seq_len(int)]]
  # }
  rslist <- unique(c(rslist, r))
  ldm <- mem_LDmatrix(rslist, pop = pop, genome_build = genome_build, token = token, ...)
  ld <- ldm[, index_snp]
  l$data$ld <- ld[match(l$data[, labs], ldm$RS_number)]
  l
}
################################################################################
# ELMO1
gene <-'ELMO1'
# subset data
barcuva.gene <- df.barcuva[df.barcuva$gene == gene,]
gtex.gene <- df.gtex[df.gtex$gene == gene,]
gwas.gene <- df.gwas[df.gwas$gene == gene,]
ibd.gene <- df.ibd[df.ibd$gene == gene,]
# set locus
loc.barcuva <- locus(barcuva.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'pval_nominal', 
                     gene = gene, fix_window = 8e5)
loc.gtex <- locus(gtex.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'pval_nominal', 
                  gene = gene, fix_window = 8e5)
loc.gwas <- locus(gwas.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'P-value', 
                  gene = gene, fix_window = 8e5)
loc.ibd <- locus(ibd.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'nom_pval', 
                 gene = gene, fix_window = 8e5)
# find GWAS index variant
index_snp <- loc.gwas$index_snp
g <- loc.gwas$data[which(loc.gwas$data$rsID == index_snp),1:6]
# add index variant, if not exists (for calculating LD)
if (! index_snp %in% loc.barcuva$data$rsID){
  loc.barcuva$data <- rbind(loc.barcuva$data, unname(unlist(cbind(g, 1, NA, NA, NA, 0)))) 
  loc.barcuva$data$pos <- as.integer(loc.barcuva$data$pos)
  loc.barcuva$data$pval_nominal <- as.double(loc.barcuva$data$pval_nominal)
  loc.barcuva$data$logP <- as.double(loc.barcuva$data$logP)
}

if (! index_snp %in% loc.gtex$data$rsID){
  loc.gtex$data <- rbind(loc.gtex$data, unname(unlist(cbind(g, 1, NA, NA, NA, 0)))) 
  loc.gtex$data$pos <- as.integer(loc.gtex$data$pos)
  loc.gtex$data$pval_nominal <- as.double(loc.gtex$data$pval_nominal)
  loc.gtex$data$logP <- as.double(loc.gtex$data$logP)
}

if (! index_snp %in% loc.ibd$data$rsID){
  loc.ibd$data <- rbind(loc.ibd$data, unname(unlist(cbind(g, 1, NA, NA, NA, 0)))) 
  loc.ibd$data$pos <- as.integer(loc.ibd$data$pos)
  loc.ibd$data$nom_pval <- as.double(loc.ibd$data$nom_pval)
  loc.ibd$data$logP <- as.double(loc.ibd$data$logP)
}
loc.barcuva$data$pos <- as.integer(loc.barcuva$data$pos)
loc.gtex$data$pos  <- as.integer(loc.gtex$data$pos)
loc.ibd$data$pos <- as.integer(loc.ibd$data$pos)
# add LD info
loc.barcuva <- get_LD(loc.barcuva, token = t, index_snp = index_snp)
loc.gtex <- get_LD(loc.gtex, token = t, index_snp = index_snp)
loc.gwas <- get_LD(loc.gwas, token = t, index_snp = index_snp)
loc.ibd <- get_LD(loc.ibd, token = t, index_snp = index_snp)
# plot
file = paste('plots/Figure_2d_LocusZoom_plot_', gene, '_GWAS.png', sep = '')
png(filename = file, res = 300, units = 'in', height = 12, width = 6)
oldpar <- set_layers(4)
scatter_plot(loc.gwas, pcutoff = NULL, index_snp = index_snp, labels = 'index',
             border = TRUE) 
title(gene, adj = 0, line = 2.4, cex.main = 2)
title('GWAS', adj = 0.1, line = 0.6, cex.main = 2)
scatter_plot(loc.ibd, pcutoff = NULL, index_snp = index_snp, labels = 'index',
             border = TRUE, legend_pos = NULL)
title('UNC', adj = 0 , line = 0.6, cex.main = 2)
scatter_plot(loc.gtex, pcutoff = NULL, index_snp = index_snp, labels = 'index', 
             border = TRUE, legend_pos = NULL)
title('GTEx', adj = 0, line = 0.6, cex.main = 2)
scatter_plot(loc.barcuva, pcutoff = NULL, index_snp = index_snp, labels = 'index',
             border = TRUE, legend_pos = NULL)
title('BarcUVa', adj = 0, line = 0.6, cex.main = 2)
genetracks(loc.gwas)
dev.off()
################################################################################
# ABO
gene <-'ABO'
# subset data
barcuva.gene <- df.barcuva[df.barcuva$gene == gene,]
gtex.gene <- df.gtex[df.gtex$gene == gene,]
gwas.gene <- df.gwas[df.gwas$gene == gene,]
ibd.gene <- df.ibd[df.ibd$gene == gene,]
# set locus
loc.barcuva <- locus(barcuva.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'pval_nominal', 
                     gene = gene, fix_window = 1e5)
loc.gtex <- locus(gtex.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'pval_nominal', 
                  gene = gene, fix_window = 1e5)
loc.gwas <- locus(gwas.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'P-value', 
                  gene = gene, fix_window = 1e5)
loc.ibd <- locus(ibd.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'nom_pval', 
                 gene = gene, fix_window = 1e5)
# find GWAS index variant
index_snp <- loc.gwas$index_snp
g <- loc.gwas$data[which(loc.gwas$data$rsID == index_snp),1:6]
# add index variant, if not exists (for calculating LD)
if (! index_snp %in% loc.barcuva$data$rsID){
  loc.barcuva$data <- rbind(loc.barcuva$data, unname(unlist(cbind(g, 1, NA, NA, NA, 0)))) 
  loc.barcuva$data$pos <- as.integer(loc.barcuva$data$pos)
  loc.barcuva$data$pval_nominal <- as.double(loc.barcuva$data$pval_nominal)
  loc.barcuva$data$logP <- as.double(loc.barcuva$data$logP)
}

if (! index_snp %in% loc.gtex$data$rsID){
  loc.gtex$data <- rbind(loc.gtex$data, unname(unlist(cbind(g, 1, NA, NA, NA, 0)))) 
  loc.gtex$data$pos <- as.integer(loc.gtex$data$pos)
  loc.gtex$data$pval_nominal <- as.double(loc.gtex$data$pval_nominal)
  loc.gtex$data$logP <- as.double(loc.gtex$data$logP)
}

if (! index_snp %in% loc.ibd$data$rsID){
  loc.ibd$data <- rbind(loc.ibd$data, unname(unlist(cbind(g, 1, NA, NA, NA, 0)))) 
  loc.ibd$data$pos <- as.integer(loc.ibd$data$pos)
  loc.ibd$data$nom_pval <- as.double(loc.ibd$data$nom_pval)
  loc.ibd$data$logP <- as.double(loc.ibd$data$logP)
}
loc.barcuva$data$pos <- as.integer(loc.barcuva$data$pos)
loc.gtex$data$pos  <- as.integer(loc.gtex$data$pos)
loc.ibd$data$pos <- as.integer(loc.ibd$data$pos)
# add LD info
loc.barcuva <- get_LD_small(loc.barcuva, token = t, index_snp = index_snp)
loc.gtex <- get_LD_small(loc.gtex, token = t, index_snp = index_snp)
loc.gwas <- get_LD(loc.gwas, token = t, index_snp = index_snp)
loc.ibd <- get_LD_small(loc.ibd, token = t, index_snp = index_snp)
# remove index SNP (not in original data)
loc.barcuva$data <- subset(loc.barcuva$data, ! rsID == index_snp)
# plot
file = paste('plots/Figure_3a_LocusZoom_plot_',
             gene, '_GWAS.png', sep = '')
png(filename = file, res = 300, units = 'in', height = 12, width = 6)
oldpar <- set_layers(4)
scatter_plot(loc.gwas, pcutoff = NULL, index_snp = index_snp, labels = 'index',
             border = TRUE, ylim = c(0,10)) 
title(gene, adj = 0, line = 2.4, cex.main = 2)
title('GWAS', adj = 0.1, line = 0.6, cex.main = 2)
scatter_plot(loc.ibd, pcutoff = NULL, index_snp = index_snp, labels = 'index',
             border = TRUE, legend_pos = NULL, ylim = c(0,45))
title('UNC', adj = 0 , line = 0.6, cex.main = 2)
scatter_plot(loc.gtex, pcutoff = NULL, index_snp = index_snp, labels = 'index', 
             border = TRUE, legend_pos = NULL, ylim = c(0,45))
title('GTEx', adj = 0, line = 0.6, cex.main = 2)
scatter_plot(loc.barcuva, pcutoff = NULL, index_snp = index_snp,
             border = TRUE, legend_pos = NULL, ylim = c(0,45))
title('BarcUVa', adj = 0, line = 0.6, cex.main = 2)
genetracks(loc.gwas)
dev.off()
################################################################################
# TNFRSF14
gene <-'TNFRSF14'
# subset data
barcuva.gene <- df.barcuva[df.barcuva$gene == gene,]
gtex.gene <- df.gtex[df.gtex$gene == gene,]
gwas.gene <- df.gwas[df.gwas$gene == gene,]
ibd.gene <- df.ibd[df.ibd$gene == gene,]
# set locus
loc.barcuva <- locus(barcuva.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'pval_nominal', 
                     gene = 'MMEL1', fix_window = 5e5)
loc.gtex <- locus(gtex.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'pval_nominal', 
                  gene = 'MMEL1', fix_window = 5e5)
loc.gwas <- locus(gwas.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'P-value', 
                  gene = 'MMEL1', fix_window = 5e5)
loc.ibd <- locus(ibd.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'nom_pval', 
                 gene = 'MMEL1', fix_window = 5e5)
# find GWAS index variant
index_snp <- loc.gwas$index_snp
g <- loc.gwas$data[which(loc.gwas$data$rsID == index_snp),1:6]
# add index variant, if not exists (for calculating LD)
if (! index_snp %in% loc.barcuva$data$rsID){
  loc.barcuva$data <- rbind(loc.barcuva$data, unname(unlist(cbind(g, 1, NA, NA, NA, 0)))) 
  loc.barcuva$data$pos <- as.integer(loc.barcuva$data$pos)
  loc.barcuva$data$pval_nominal <- as.double(loc.barcuva$data$pval_nominal)
  loc.barcuva$data$logP <- as.double(loc.barcuva$data$logP)
}

if (! index_snp %in% loc.gtex$data$rsID){
  loc.gtex$data <- rbind(loc.gtex$data, unname(unlist(cbind(g, 1, NA, NA, NA, 0)))) 
  loc.gtex$data$pos <- as.integer(loc.gtex$data$pos)
  loc.gtex$data$pval_nominal <- as.double(loc.gtex$data$pval_nominal)
  loc.gtex$data$logP <- as.double(loc.gtex$data$logP)
}

if (! index_snp %in% loc.ibd$data$rsID){
  loc.ibd$data <- rbind(loc.ibd$data, unname(unlist(cbind(g, 1, NA, NA, NA, 0)))) 
  loc.ibd$data$pos <- as.integer(loc.ibd$data$pos)
  loc.ibd$data$nom_pval <- as.double(loc.ibd$data$nom_pval)
  loc.ibd$data$logP <- as.double(loc.ibd$data$logP)
}
loc.barcuva$data$pos <- as.integer(loc.barcuva$data$pos)
loc.gtex$data$pos  <- as.integer(loc.gtex$data$pos)
loc.ibd$data$pos <- as.integer(loc.ibd$data$pos)
# add LD info
loc.barcuva <- get_LD_small(loc.barcuva, token = t, index_snp = index_snp)
loc.gtex <- get_LD(loc.gtex, token = t, index_snp = index_snp)
loc.gwas <- get_LD(loc.gwas, token = t, index_snp = index_snp)
loc.ibd <- get_LD(loc.ibd, token = t, index_snp = index_snp)
# plot
file = paste('plots/Figure_3d_LocusZoom_plot_',
             gene, '_GWAS.png', sep = '')
png(filename = file, res = 300, units = 'in', height = 12, width = 6)
oldpar <- set_layers(4)
scatter_plot(loc.gwas, pcutoff = NULL, index_snp = index_snp, labels = 'index',
             border = TRUE, ylim = c(0,17)) 
title(gene, adj = 0, line = 2.4, cex.main = 2)
title('GWAS', adj = 0.1, line = 0.6, cex.main = 2)
scatter_plot(loc.ibd, pcutoff = NULL, index_snp = index_snp, labels = 'index',
             border = TRUE, legend_pos = NULL, ylim = c(0,14), label_y = 13, label_x = 1)
title('UNC', adj = 0 , line = 0.6, cex.main = 2)
scatter_plot(loc.gtex, pcutoff = NULL, index_snp = index_snp, labels = 'index', 
             border = TRUE, legend_pos = NULL, ylim = c(0,14), label_y = 10, label_x = 1)
title('GTEx', adj = 0, line = 0.6, cex.main = 2)
scatter_plot(loc.barcuva, pcutoff = NULL, index_snp = index_snp, labels = 'index',
             border = TRUE, legend_pos = NULL, ylim = c(0,14), label_y = 10, label_x = 1)
title('BarcUVa', adj = 0, line = 0.6, cex.main = 2)
genetracks(loc.gwas)
dev.off()
