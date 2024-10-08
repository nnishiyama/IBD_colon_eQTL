# Figure 1c. LocusZoom plot H3 eQTL-eQTL colocalization example
library(tidyverse)
library(data.table)
library(locuszoomr)
library(LDlinkR)
library(AnnotationHub)
library(memoise)

# load data
load('data/UNC/RData/locusZoom_matrices_eQTL-eQTL_H3.RData')
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
  int <- 1800
  if (length(rslist > int)) {
    rslist <- rslist[order(l$data$logP, decreasing = TRUE)[seq_len(int)]]
  }
  rslist <- unique(c(rslist, r))
  ldm <- mem_LDmatrix(rslist, pop = pop, genome_build = genome_build, token = token, ...)
  ld <- ldm[, index_snp]
  l$data$ld <- ld[match(l$data[, labs], ldm$RS_number)]
  l
}
################################################################################
# SLC45A4
gene <-'SLC45A4'
# subset data
barcuva.gene <- df.barcuva[df.barcuva$gene == gene,]
gtex.gene <- df.gtex[df.gtex$gene == gene,]
ibd.gene <- df.ibd[df.ibd$gene == gene,]
# set locus
loc.barcuva <- locus(barcuva.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'pval_nominal', 
                     gene = 'PTP4A3', fix_window = 7e5)
loc.gtex <- locus(gtex.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'pval_nominal', 
                  gene = 'PTP4A3', fix_window = 7e5)
loc.ibd <- locus(ibd.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'nom_pval', 
                 gene = 'PTP4A3', fix_window = 7e5)
# find IBD eQTL index variant
index_snp <- loc.ibd$index_snp
g <- loc.ibd$data[which(loc.ibd$data$rsID == index_snp),1:6]
# add index variant, if not exists (for calculating LD)
if (! index_snp %in% loc.barcuva$data$rsID) {
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
loc.barcuva$data$pos <- as.integer(loc.barcuva$data$pos)
loc.gtex$data$pos  <- as.integer(loc.gtex$data$pos)
loc.ibd$data$pos <- as.integer(loc.ibd$data$pos)
# add LD info
loc.barcuva <- get_LD(loc.barcuva, token = t, index_snp = index_snp)
loc.gtex <- get_LD(loc.gtex, token = t, index_snp = index_snp)
loc.ibd <- get_LD(loc.ibd, token = t, index_snp = index_snp)
# plot
file = paste('plot/Figure_1c_LocusZoom_plot_', gene, '_eQTL.png', sep = '')
png(filename = file, res = 300, units = 'in', height = 12, width = 6)
oldpar <- set_layers(3)
try(scatter_plot(loc.ibd, pcutoff = NULL, index_snp = index_snp, labels = 'index',
                 border = TRUE, label_x = 2, label_y = 0.5, ylim = c(0,6.1)), silent = TRUE)
title(gene, adj = 0, line = 2.4, cex.main = 2)
title('UNC', adj = 1 , line = 0.6, cex.main = 2)
try(scatter_plot(loc.gtex, pcutoff = NULL, index_snp = index_snp, labels = 'index',
                 border = TRUE, ylim = c(0,6.1), legend_pos = NULL), silent = TRUE)
title('GTEx', adj = 1, line = 0.6, cex.main = 2)
try(scatter_plot(loc.barcuva, pcutoff = NULL, index_snp = index_snp, labels = 'index',
                 border = TRUE, ylim = c(0,6.1), legend_pos = NULL), silent = TRUE)
title('BarcUVa', adj = 1, line = 0.6, cex.main = 2)
genetracks(loc.ibd)
dev.off()
