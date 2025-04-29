# Figure 5b-h.
library(data.table)
library(tidyverse)
library(mashr)
library(qvalue)
library(ggplot2)
library(ggsignif)
library(locuszoomr)
library(LDlinkR)
library(AnnotationHub)
library(memoise)
library(ggsignif)
library(UpSetR)
library(ComplexHeatmap)

################################################################################
# Figure 5b. violin plots of unadjusted eQTL-eQTL colocalization effect sizes

# load full summary stats
load('data/UNC/RData/coloc_all_data.RData')
remove(cd, ibd, uc)
# remove duplicates
eqtl.barcuva <- eqtl.barcuva %>% distinct(gene_id, snp, .keep_all = TRUE)
eqtl.gtex <- eqtl.gtex %>% distinct(gene_id, snp, .keep_all = TRUE)
eqtl.ibd <- eqtl.ibd %>% distinct(phe_id, snp, .keep_all = TRUE)
# calculate cv
eqtl.barcuva <- eqtl.barcuva %>% mutate(cv = abs(slope/(slope_se*sqrt(445))))
eqtl.gtex <- eqtl.gtex %>% mutate(cv = abs(slope/(slope_se*sqrt(368))))
eqtl.ibd <- eqtl.ibd %>% mutate(cv = abs(slope/(slope_se*sqrt(252))))
# filter by cv
eqtl.barcuva <- eqtl.barcuva[eqtl.barcuva$cv < 1,]
eqtl.gtex <- eqtl.gtex[eqtl.gtex$cv < 1,]
eqtl.ibd <- eqtl.ibd[eqtl.ibd$cv < 1,]
# load eQTL-eQTL colocalization results
coloc.barcuva.unc <- fread('data/UNC/Supplementary_Table_5_coloc_IBD_eQTL_BarcUVa_eQTL.txt')
coloc.gtex.unc <- fread('data/UNC/Supplementary_Table_4_coloc_IBD_eQTL_GTEx_eQTL.txt')
coloc.nibd <- fread('data/GTEx/Supplementary_Table_6_coloc_GTEx_eQTL_BarcUVa_eQTL.txt')
# separate out ENSEMBL IDs
coloc.barcuva.unc <- coloc.barcuva.unc %>% separate(gene_id_IBD, c('gene', 'ver'), remove = FALSE)
coloc.gtex.unc <- coloc.gtex.unc %>% separate(gene_id_IBD, c('gene', 'ver'), remove = FALSE)
coloc.nibd <- coloc.nibd %>% separate(gene_id_BarcUVa, c('gene', 'ver'), remove = FALSE)
# limit results to genes with a significant eQTL in IBD tissue
coloc.barcuva.unc <- coloc.barcuva.unc[coloc.barcuva.unc$qval_IBD < 0.05,]
coloc.gtex.unc <- coloc.gtex.unc[coloc.gtex.unc$qval_IBD < 0.05,]
coloc.nibd <- coloc.nibd[coloc.nibd$qval_BarcUVa < 0.05 | coloc.nibd$qval_GTEx < 0.05,]
# create a gene id key
key1 <- coloc.barcuva.unc %>% dplyr::select(gene, gene_id_BarcUVa, gene_id_IBD, nom_pval_IBD, nom_pval_BarcUVa)
key2 <- coloc.gtex.unc %>% dplyr::select(gene, gene_id_GTEx, gene_id_IBD, nom_pval_IBD, nom_pval_GTEx)
key3 <- coloc.nibd %>% dplyr::select(gene, gene_id_BarcUVa, gene_id_GTEx, nom_pval_BarcUVa, nom_pval_GTEx)
gene.key <- full_join(key1, key2)
gene.key <- full_join(gene.key, key3)
# set PP thresholds
h1.b <- coloc.barcuva.unc[coloc.barcuva.unc$H1_hypothesis_IBD > 0.5,]
h2.b <- coloc.barcuva.unc[coloc.barcuva.unc$H2_hypothesis_BarcUVa > 0.5,]
h3.b <- coloc.barcuva.unc[coloc.barcuva.unc$H3_hypothesis_different > 0.5,]
h4.b <- coloc.barcuva.unc[coloc.barcuva.unc$H4_hypothesis_shared > 0.5,]
h1.g <- coloc.gtex.unc[coloc.gtex.unc$H1_hypothesis_IBD > 0.5,]
h2.g <- coloc.gtex.unc[coloc.gtex.unc$H2_hypothesis_GTEx > 0.5,]
h3.g <- coloc.gtex.unc[coloc.gtex.unc$H3_hypothesis_different > 0.5,]
h4.g <- coloc.gtex.unc[coloc.gtex.unc$H4_hypothesis_shared > 0.5,]
h1.n <- coloc.nibd[coloc.nibd$H1_hypothesis_GTEx > 0.5,]
h2.n <- coloc.nibd[coloc.nibd$H2_hypothesis_BarcUVa > 0.5,]
h3.n <- coloc.nibd[coloc.nibd$H3_hypothesis_different > 0.5,]
h4.n <- coloc.nibd[coloc.nibd$H4_hypothesis_shared > 0.5,]
# pull out genes
barcuva.unique.genes <- subset(gene.key, gene %in% setdiff(c(h2.b$gene, h3.b$gene, h2.n$gene, h3.n$gene),
                                c(h1.n$gene, h4.n$gene, h1.g$gene, h2.g$gene, h3.g$gene, h4.g$gene)))$gene_id_BarcUVa
gtex.unique.genes <- subset(gene.key, gene %in% setdiff(c(h2.g$gene, h3.g$gene, h1.n$gene, h3.n$gene), 
                             c(h2.n$gene, h4.n$gene, h1.b$gene, h2.b$gene, h3.b$gene, h4.b$gene)))$gene_id_GTEx
ibd.unique.genes <- setdiff(c(h1.b$gene_id_IBD, h3.b$gene_id_IBD, h1.g$gene_id_IBD, h3.g$gene_id_IBD),
                            c(h2.b$gene_id_IBD, h4.b$gene_id_IBD, h2.g$gene_id_IBD, h4.g$gene_id_IBD))
barcuva.shared.genes <- unique(c(h4.b$gene_id_BarcUVa, h4.n$gene_id_BarcUVa))
gtex.shared.genes <- unique(c(h4.g$gene_id_GTEx, h4.n$gene_id_GTEx))
ibd.shared.genes <- unique(c(h4.b$gene_id_IBD, h4.g$gene_id_IBD))
# pull out full signals
barcuva.eqtl.unique <- subset(eqtl.barcuva, gene_id %in% barcuva.unique.genes)
gtex.eqtl.unique <- subset(eqtl.gtex, gene_id %in% gtex.unique.genes)
ibd.eqtl.unique <- subset(eqtl.ibd, phe_id %in% ibd.unique.genes)
barcuva.eqtl.shared <- subset(eqtl.barcuva, gene_id %in% barcuva.shared.genes)
gtex.eqtl.shared <- subset(eqtl.gtex, gene_id %in% gtex.shared.genes)
ibd.eqtl.shared <- subset(eqtl.ibd, phe_id %in% ibd.shared.genes)
# pull out nominal p-value theshold
ibd.thresh <- gene.key %>% select(gene_id_IBD, nom_pval_IBD) %>% 
  rename(phe_id = gene_id_IBD)
gtex.thresh <- gene.key %>% select(gene_id_GTEx, nom_pval_GTEx) %>%
  rename(gene_id = gene_id_GTEx)
barcuva.thresh <- gene.key %>% select(gene_id_BarcUVa, nom_pval_BarcUVa) %>%
  rename(gene_id = gene_id_BarcUVa)
# filter by nominal p-value threshold
barcuva.eqtl.shared <- barcuva.eqtl.shared %>% left_join(barcuva.thresh, by = 'gene_id') %>%
  filter(pval_nominal <= nom_pval_BarcUVa) %>% select(-c(nom_pval_BarcUVa))
barcuva.eqtl.unique <- barcuva.eqtl.unique %>% left_join(barcuva.thresh, by = 'gene_id') %>%
  filter(pval_nominal <= nom_pval_BarcUVa) %>% select(-c(nom_pval_BarcUVa))
gtex.eqtl.shared <- gtex.eqtl.shared %>% left_join(gtex.thresh, by = 'gene_id') %>%
  filter(pval_nominal <= nom_pval_GTEx) %>% select(-c(nom_pval_GTEx))
gtex.eqtl.unique <- gtex.eqtl.unique %>% left_join(gtex.thresh, by = 'gene_id') %>%
  filter(pval_nominal <= nom_pval_GTEx) %>% select(-c(nom_pval_GTEx))
ibd.eqtl.shared <- ibd.eqtl.shared %>% left_join(ibd.thresh, by = 'phe_id') %>% 
  filter(nom_pval <= nom_pval_IBD) %>% select(-c(nom_pval_IBD))
ibd.eqtl.unique <- ibd.eqtl.unique %>% left_join(ibd.thresh, by = 'phe_id') %>% 
  filter(nom_pval <= nom_pval_IBD) %>% select(-c(nom_pval_IBD))
# format for plotting
ibd.eqtl.shared$group <- 'UNC'
ibd.eqtl.unique$group <- 'UNC'
gtex.eqtl.shared$group <- 'GTEx'
gtex.eqtl.unique$group <- 'GTEx'
barcuva.eqtl.shared$group <- 'BarcUVa'
barcuva.eqtl.unique$group <- 'BarcUVa'
unique.eqtl <- full_join(barcuva.eqtl.unique, gtex.eqtl.unique)
unique.eqtl <- full_join(unique.eqtl, ibd.eqtl.unique)
unique.eqtl$eGene <- 'unique'
shared.eqtl <- full_join(barcuva.eqtl.shared, gtex.eqtl.shared)
shared.eqtl <- full_join(shared.eqtl, ibd.eqtl.shared)
shared.eqtl$eGene <- 'shared in two or more'
p.eqtl <- full_join(unique.eqtl, shared.eqtl)
# convert to log10 scale 
p.eqtl <- p.eqtl %>% mutate(log10.abs = log10(abs(slope)))
# summary stats
ss.eqtl <- p.eqtl %>% group_by(group, eGene) %>% summarize(n=n(), med_es=round(median(log10.abs),3),
                                                           yn = -1.35, ym = -1.3)
# stats
wilcox.test(abs(barcuva.eqtl.shared$slope), abs(barcuva.eqtl.unique$slope))
wilcox.test(abs(gtex.eqtl.shared$slope), abs(gtex.eqtl.unique$slope))
wi <- wilcox.test(abs(ibd.eqtl.shared$slope), abs(ibd.eqtl.unique$slope))
p.stats <- c(rep('p < 2.2e-16', 2), paste('p =', round(wi$p.value, 3)))
# plot effect size
png(res = 300, units = 'in', height = 8, width = 10,
    filename = 'plots/Figure_5b_violin_plot_unadjusted_effect_size_eQTL_colocalizing_full_signal_log10.png')
ggplot(p.eqtl, aes(group, log10.abs, fill = eGene)) + geom_violin(lwd=1) + 
  geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) +
  labs(x='', y='log10(absoulte unadjusted efffect size)', title = 'eQTL-eQTL colocalizations: unadjusted effect sizes') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 16)) + 
  geom_text(data = ss.eqtl, aes(group, y = ym, label = paste0('med = ', med_es)), 
            position = position_dodge(width = 0.9), size = 4) +
  geom_text(dat = ss.eqtl, aes(group, y = yn, label = paste0('n = ', n)), 
            position = position_dodge(width = 0.9), size = 4) +
  geom_signif(y_position=c(0.4, 0.4, 0.25), xmin=c(0.77, 1.77, 2.77), xmax=c(1.22, 2.22, 3.22),
              annotation=p.stats, tip_length=0)
dev.off()
################################################################################
# Figure 5c. violin plots of mash-adjusted eQTL-eQTL colocalization effect sizes

# load mash results
load('data/UNC/RData/mashr_all_fit_full_signal.RData')
# reformat
all.sd$eqtl <- all.beta$eqtl
# subset sd's
barcuva.sd <- subset(all.sd, eqtl %in% ma.b$eqtl) %>% arrange(match(eqtl, ma.b$eqtl)) %>%
  mutate(se = barcuva/sqrt(445))
gtex.sd <- subset(all.sd, eqtl %in% ma.g$eqtl) %>% arrange(match(eqtl, ma.g$eqtl)) %>%
  mutate(se = gtex/sqrt(368))
ibd.sd <- subset(all.sd, eqtl %in% ma.i$eqtl) %>% arrange(match(eqtl, ma.i$eqtl)) %>% 
  mutate(se = ibd/sqrt(252))
# calculate CVs
ma.b <- ma.b %>% mutate(cv = abs(barcuva.sd$barcuva/barcuva))
ma.g <- ma.g %>% mutate(cv = abs(gtex.sd$gtex/gtex))
ma.i <- ma.i %>% mutate(cv = abs(ibd.sd$ibd/ibd))
# filter out results with high CV (>1)
ma.b <- ma.b[ma.b$cv < 1,]
ma.g <- ma.g[ma.g$cv < 1,]
ma.i <- ma.i[ma.i$cv < 1,]
# pull out genes
mash.barcuva.unique.genes <- setdiff(c(h2.b$gene, h3.b$gene, h2.n$gene, h3.n$gene),
                                                           c(h1.n$gene, h4.n$gene, h1.g$gene, h2.g$gene, h3.g$gene, h4.g$gene))
mash.gtex.unique.genes <- setdiff(c(h2.g$gene, h3.g$gene, h1.n$gene, h3.n$gene), 
                                                        c(h2.n$gene, h4.n$gene, h1.b$gene, h2.b$gene, h3.b$gene, h4.b$gene))
mash.ibd.unique.genes <- setdiff(c(h1.b$gene, h3.b$gene, h1.g$gene, h3.g$gene),
                            c(h2.b$gene, h4.b$gene, h2.g$gene, h4.g$gene))
mash.barcuva.shared.genes <- unique(c(h4.b$gene, h4.n$gene))
mash.gtex.shared.genes <- unique(c(h4.g$gene, h4.n$gene))
mash.ibd.shared.genes <- unique(c(h4.b$gene, h4.g$gene))
# pull out full signals
mash.barcuva.eqtl.unique <- subset(ma.b, gene %in% mash.barcuva.unique.genes) %>% dplyr::select(-c(gtex, ibd))
mash.gtex.eqtl.unique <- subset(ma.g, gene %in% mash.gtex.unique.genes) %>% dplyr::select(-c(barcuva, ibd))
mash.ibd.eqtl.unique <- subset(ma.i, gene %in% mash.ibd.unique.genes) %>% dplyr::select(-c(gtex, barcuva))
mash.barcuva.eqtl.shared <- subset(ma.b, gene %in% mash.barcuva.shared.genes) %>% dplyr::select(-c(gtex, ibd))
mash.gtex.eqtl.shared <- subset(ma.g, gene %in% mash.gtex.shared.genes) %>% dplyr::select(-c(barcuva, ibd))
mash.ibd.eqtl.shared <- subset(ma.i, gene %in% mash.ibd.shared.genes) %>% dplyr::select(-c(gtex, barcuva))
# format for plotting
mash.barcuva.eqtl.shared$group <- 'BarcUVa'
mash.barcuva.eqtl.unique$group <- 'BarcUVa'
mash.gtex.eqtl.shared$group <- 'GTEx'
mash.gtex.eqtl.unique$group <- 'GTEx'
mash.ibd.eqtl.shared$group <- 'UNC'
mash.ibd.eqtl.unique$group <- 'UNC'
colnames(mash.barcuva.eqtl.shared)[1] <- 'slope'
colnames(mash.barcuva.eqtl.unique)[1] <- 'slope'
colnames(mash.gtex.eqtl.shared)[1] <- 'slope'
colnames(mash.gtex.eqtl.unique)[1] <- 'slope'
colnames(mash.ibd.eqtl.shared)[1] <- 'slope'
colnames(mash.ibd.eqtl.unique)[1] <- 'slope'
mash.eqtl.shared <- full_join(mash.barcuva.eqtl.shared, mash.gtex.eqtl.shared)
mash.eqtl.shared <- full_join(mash.eqtl.shared, mash.ibd.eqtl.shared)
mash.eqtl.shared$eGene <- 'shared in two or more'
mash.eqtl.unique <- full_join(mash.barcuva.eqtl.unique, mash.gtex.eqtl.unique)
mash.eqtl.unique <- full_join(mash.eqtl.unique, mash.ibd.eqtl.unique)
mash.eqtl.unique$eGene <- 'unique'
p.mash.eqtl <- full_join(mash.eqtl.shared, mash.eqtl.unique)
# convert to log10 scale
p.mash.eqtl <- p.mash.eqtl %>% mutate(log10.abs = log10(abs(slope)))
# summary stats
ss.mash.eqtl <- p.mash.eqtl %>% group_by(group, eGene) %>% summarize(n=n(), med_es=round(median(log10.abs),3),
                                                           yn = -1.7, ym = -1.65)
# stats
wilcox.test(log10(abs(mash.barcuva.eqtl.shared$slope)), log10(abs(mash.barcuva.eqtl.unique$slope)))
wilcox.test(abs(mash.gtex.eqtl.shared$slope), abs(mash.gtex.eqtl.unique$slope))
wilcox.test(abs(mash.ibd.eqtl.shared$slope), abs(mash.ibd.eqtl.unique$slope))
p.stats <- rep('p < 2.2e-16', 3)
# plot effect size
png(res = 300, units = 'in', height = 8, width = 10,
    filename = 'plots/Figure_5c_violin_plot_mash_effect_size_eQTL_colocalizing_full_signal_log10.png')
ggplot(p.mash.eqtl, aes(group, log10.abs, fill = eGene)) + geom_violin(lwd=1) +
  geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) +
  labs(x='', y='log10(absoulte mash-adjusted efffect size)', title = 'eQTL-eQTL colocalizations: mash-adjusted effect sizes') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 16)) + 
  geom_text(data = ss.mash.eqtl, aes(group, y = ym, label = paste0('med = ', med_es)), 
            position = position_dodge(width = 0.9), size = 4) +
  geom_text(dat = ss.mash.eqtl, aes(group, y = yn, label = paste0('n = ', n)), 
            position = position_dodge(width = 0.9), size = 4) +
  geom_signif(y_position=c(0.4, 0.4, 0.35), xmin=c(0.77, 1.77, 2.77), xmax=c(1.22, 2.22, 3.22),
              annotation=p.stats, tip_length=0) 
dev.off()
################################################################################
# Figure 5d. violin plots of unadjusted effect sizes: full signal GWAS colocalizations

# load GWAS results
coloc.barcuva <- fread('data/BarcUVa/Supplementary_Table_9_coloc_H4_BarcUVa_eQTL_GWAS.txt')
coloc.gtex <- fread('data/GTEx/Supplementary_Table_8_coloc_H4_GTEx_eQTL_colon_transverse_GWAS.txt')
coloc.ibd <- fread('data/UNC/Supplementary_Table_7_coloc_H4_IBD_eQTL_GWAS.txt')
# reformat gene ids
coloc.ibd <- coloc.ibd %>% separate(gene, c('gene', 'ver'))
coloc.barcuva <- coloc.barcuva %>% separate(gene, c('gene', 'ver'))
coloc.gtex <- coloc.gtex %>% separate(gene, c('gene', 'ver'))
# UpSet plot for colocalizing eGenes
set.seed(12345)
# create a list of colocalizing eGenes
set.egenes <- list(UNC = coloc.ibd$gene, 
                   GTEx = coloc.gtex$gene,
                   BarcUVa = coloc.barcuva$gene)
# create combination matrix
mat.egenes <- make_comb_mat(set.egenes, mode = 'distinct')
# pull out colocalizing eGenes unique to IBD
unique.egenes <- extract_comb(mat.egenes, '100')
barcuva.unique <- subset(coloc.barcuva, gene %in% extract_comb(mat.egenes, '001'))
gtex.unique <- subset(coloc.gtex, gene %in% extract_comb(mat.egenes, '010'))
ibd.unique <- subset(coloc.ibd, gene %in% unique.egenes)
# combine in shared in at least one
ibd.comb <- subset(coloc.ibd, ! gene %in% ibd.unique$gene)
gtex.comb <- subset(coloc.gtex, ! gene %in% gtex.unique$gene)
barcuva.comb <- subset(coloc.barcuva, ! gene %in% barcuva.unique$gene)
# reformat ENSEMBL IDs
ibd.comb <- ibd.comb %>% unite('phe_id', c(gene, ver), sep = '.', remove = FALSE)
gtex.comb <- gtex.comb %>% unite('gene_id', c(gene, ver), sep = '.', remove = FALSE)
barcuva.comb <- barcuva.comb %>% unite('gene_id', c(gene, ver), sep = '.', remove = FALSE)
ibd.unique <- ibd.unique %>% unite('phe_id', c(gene, ver), sep = '.', remove = FALSE)
gtex.unique <- gtex.unique %>% unite('gene_id', c(gene, ver), sep = '.', remove = FALSE)
barcuva.unique <- barcuva.unique %>% unite('gene_id', c(gene, ver), sep = '.', remove = FALSE)
# subset full summary statistics
ibd.full.comb <- subset(eqtl.ibd, phe_id %in% ibd.comb$phe_id)
gtex.full.comb <- subset(eqtl.gtex, gene_id %in% gtex.comb$gene_id)
barcuva.full.comb <- subset(eqtl.barcuva, gene_id %in% barcuva.comb$gene_id)
ibd.full.unique <- subset(eqtl.ibd, phe_id %in% ibd.unique$phe_id)
gtex.full.unique <- subset(eqtl.gtex, gene_id %in% gtex.unique$gene_id)
barcuva.full.unique <- subset(eqtl.barcuva, gene_id %in% barcuva.unique$gene_id)
# pull out nominal p-value thresholds
ibd.thresh <- coloc.ibd %>% unite('phe_id', c(gene, ver), sep = '.', remove = FALSE) %>% 
  select(phe_id, nom_pval_thresh)
barcuva.thresh <- coloc.barcuva %>% unite('gene_id', c(gene, ver), sep = '.', remove = FALSE) %>% 
  select(gene_id, pval_nominal_threshold)
gtex.thresh <- coloc.gtex %>% unite('gene_id', c(gene, ver), sep = '.', remove = FALSE) %>% 
  select(gene_id, pval_nominal_threshold)
# filter based on nominal p-value thresholds
ibd.full.comb <- ibd.full.comb %>% left_join(ibd.thresh, by = 'phe_id') %>% 
  filter(nom_pval <= nom_pval_thresh) %>% select(-c(nom_pval_thresh))
gtex.full.comb <- gtex.full.comb %>% left_join(gtex.thresh, by = 'gene_id') %>% 
  filter(pval_nominal <= pval_nominal_threshold) %>% select(-c(pval_nominal_threshold))
barcuva.full.comb <- barcuva.full.comb %>% left_join(barcuva.thresh, by = 'gene_id') %>% 
  filter(pval_nominal <= pval_nominal_threshold) %>% select(-c(pval_nominal_threshold))
ibd.full.unique <- ibd.full.unique %>% left_join(ibd.thresh, by = 'phe_id') %>% 
  filter(nom_pval <= nom_pval_thresh) %>% select(-c(nom_pval_thresh))
gtex.full.unique <- gtex.full.unique %>% left_join(gtex.thresh, by = 'gene_id') %>% 
  filter(pval_nominal <= pval_nominal_threshold) %>% select(-c(pval_nominal_threshold))
barcuva.full.unique <- barcuva.full.unique %>% left_join(barcuva.thresh, by = 'gene_id') %>% 
  filter(pval_nominal <= pval_nominal_threshold) %>% select(-c(pval_nominal_threshold))
# format for plotting
ibd.full.comb$group <- 'UNC'
ibd.full.unique$group <- 'UNC'
gtex.full.comb$group <- 'GTEx'
gtex.full.unique$group <- 'GTEx'
barcuva.full.comb$group <- 'BarcUVa'
barcuva.full.unique$group <- 'BarcUVa'
unique.full <- full_join(barcuva.full.unique, gtex.full.unique)
unique.full <- full_join(unique.full, ibd.full.unique)
unique.full$eGene <- 'unique'
comb.full <- full_join(barcuva.full.comb, gtex.full.comb)
comb.full <- full_join(comb.full, ibd.full.comb)
comb.full$eGene <- 'shared in two or more'
p.full <- full_join(unique.full, comb.full)
# convert to log10 scale
p.full <- p.full %>% mutate(log10.abs = log10(abs(slope)))
# summary stats
ss.full <- p.full %>% group_by(group, eGene) %>% summarize(n=n(), med_es=round(median(log10.abs),3),
            yn = -1.2, ym = -1.15)
# stats
wb <- wilcox.test(abs(barcuva.full.comb$slope), abs(barcuva.full.unique$slope))
wg <- wilcox.test(abs(gtex.full.comb$slope), abs(gtex.full.unique$slope))
wi <- wilcox.test(abs(ibd.full.comb$slope), abs(ibd.full.unique$slope))
#p.stats <- c(wb$p.value, wg$p.value, wi$p.value)
#p.stats <- gsub('^', 'p = ', p.stats)
#p.stats <- rep('p < 2.2e-16', 3)
p.stats <- c('p < 2.2e-16', paste('p =', round(wg$p.value, 3)), 'p < 2.2e-16')
# plot effect size
png(res = 300, units = 'in', height = 8, width = 10,
    filename = 'plots/Figure_5d_violin_plot_unadjusted_effect_size_GWAS_colocalizing_full_signal_log10.png')
ggplot(p.full, aes(group, log10.abs, fill = eGene)) + geom_violin(lwd=1) +  
  geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) +
  labs(x='', y='log10(absoulte unadjusted efffect size)', title = 'GWAS-eQTL colocalizations: unadjusted effect sizes') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 16)) + 
  geom_text(data = ss.full, aes(group, y = ym, label = paste0('med = ', med_es)), 
            position = position_dodge(width = 0.9), size = 4) +
  geom_text(dat = ss.full, aes(group, y = yn, label = paste0('n = ', n)), 
            position = position_dodge(width = 0.9), size = 4) +
  geom_signif(y_position=c(0.25, 0.25, 0.15), xmin=c(0.77, 1.77, 2.77), xmax=c(1.22, 2.22, 3.22),
              annotation=p.stats, tip_length=0)
dev.off()
################################################################################
# Figure 5e. violin plots of mash-adjusted effect sizes: full signal GWAS colocalizations

# load mash results
load(file = 'data/UNC/RData/mashr_coloc_fit_full_signal.RData')
# format to pull standard errors
coloc.sd$rsID <- df.coloc$rsID
coloc.sd$eqtl <- df.coloc$eqtl
coloc.sd <- coloc.sd %>% separate(eqtl, c('gene', 'snp'), sep = '-', remove = FALSE) %>%
  distinct(eqtl, .keep_all = TRUE)
# keep only distinct results
mc.b <- mc.b %>% distinct(eqtl, .keep_all = TRUE)
mc.g <- mc.g %>% distinct(eqtl, .keep_all = TRUE)
mc.i <- mc.i %>% distinct(eqtl, .keep_all = TRUE)
# subset sd's
barcuva.sd <- subset(coloc.sd, eqtl %in% mc.b$eqtl) %>% arrange(match(eqtl, mc.b$eqtl)) %>%
  mutate(se = barcuva/sqrt(445))
gtex.sd <- subset(coloc.sd, eqtl %in% mc.g$eqtl) %>% arrange(match(eqtl, mc.g$eqtl)) %>%
  mutate(se = gtex/sqrt(368))
ibd.sd <- subset(coloc.sd, eqtl %in% mc.i$eqtl) %>% arrange(match(eqtl, mc.i$eqtl)) %>% 
  mutate(se = ibd/sqrt(252))
# calculate CVs
mc.b <- mc.b %>% mutate(cv = abs(barcuva.sd$barcuva/barcuva))
mc.g <- mc.g %>% mutate(cv = abs(gtex.sd$gtex/gtex))
mc.i <- mc.i %>% mutate(cv = abs(ibd.sd$ibd/ibd))
# filter out results with high CV (>1)
mc.b <- mc.b[mc.b$cv < 1,]
mc.g <- mc.g[mc.g$cv < 1,]
mc.i <- mc.i[mc.i$cv < 1,]
# subset mash results
mash.barcuva.comb <- subset(mc.b, gene %in% barcuva.comb$gene) %>% dplyr::select(-c(gtex, ibd))
mash.gtex.comb <- subset(mc.g, gene %in% gtex.comb$gene) %>% dplyr::select(-c(barcuva, ibd))
mash.ibd.comb <- subset(mc.i, gene %in% ibd.comb$gene) %>% dplyr::select(-c(gtex, barcuva))
mash.barcuva.unique <- subset(mc.b, gene %in% barcuva.unique$gene) %>% dplyr::select(-c(gtex, ibd))
mash.gtex.unique <- subset(mc.g, gene %in% gtex.unique$gene) %>% dplyr::select(-c(barcuva, ibd))
mash.ibd.unique <- subset(mc.i, gene %in% ibd.unique$gene) %>% dplyr::select(-c(gtex, barcuva))
# format for plotting
mash.barcuva.comb$group <- 'BarcUVa'
mash.barcuva.unique$group <- 'BarcUVa'
mash.gtex.comb$group <- 'GTEx'
mash.gtex.unique$group <- 'GTEx'
mash.ibd.comb$group <- 'UNC'
mash.ibd.unique$group <- 'UNC'
colnames(mash.barcuva.comb)[1] <- 'slope'
colnames(mash.barcuva.unique)[1] <- 'slope'
colnames(mash.gtex.comb)[1] <- 'slope'
colnames(mash.gtex.unique)[1] <- 'slope'
colnames(mash.ibd.comb)[1] <- 'slope'
colnames(mash.ibd.unique)[1] <- 'slope'
mash.comb <- full_join(mash.barcuva.comb, mash.gtex.comb)
mash.comb <- full_join(mash.comb, mash.ibd.comb)
mash.comb$eGene <- 'shared in two or more'
mash.unique <- full_join(mash.barcuva.unique, mash.gtex.unique)
mash.unique <- full_join(mash.unique, mash.ibd.unique)
mash.unique$eGene <- 'unique'
p.mash <- full_join(mash.comb, mash.unique)
# convert to log10 scale
p.mash <- p.mash %>% mutate(log10.abs = log10(abs(slope)))
# summary stats
ss.mash <- p.mash %>% group_by(group, eGene) %>% summarize(n=n(), med_es=round(median(log10.abs),3),
          yn = -1.65, ym = -1.6)
# stats
wb <- wilcox.test(abs(mash.barcuva.comb$slope), abs(mash.barcuva.unique$slope))
wg <- wilcox.test(abs(mash.gtex.comb$slope), abs(mash.gtex.unique$slope))
wi <- wilcox.test(abs(mash.ibd.comb$slope), abs(mash.ibd.unique$slope))
wilcox.test(abs(mash.gtex.comb$slope), abs(mash.ibd.comb$slope))
p.stats <- c('p < 2.2e-16', gsub('^', 'p = ', round(wg$p.value, 3)), 'p < 2.2e-16')
# plot effect size
png(res = 300, units = 'in', height = 8, width = 10,
    filename = 'plots/Figure_5e_violin_plot_mash_effect_size_GWAS_colocalizing_full_signal_log10.png')
ggplot(p.mash, aes(group, log10.abs, fill = eGene)) + geom_violin(lwd=1) + 
  geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) +
  labs(x='', y='log10(absoulte mash-adjusted efffect size)', title = 'GWAS-eQTL colocalizations: mash-adjusted effect sizes') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 16)) + 
geom_text(data = ss.mash, aes(group, y = ym, label = paste0('med = ', med_es)), 
          position = position_dodge(width = 0.9), size = 4) +
  geom_text(dat = ss.mash, aes(group, y = yn, label = paste0('n = ', n)), 
            position = position_dodge(width = 0.9), size = 4) +
  geom_signif(y_position=c(0.2, 0.08, 0.1), xmin=c(0.77, 1.77, 2.77), xmax=c(1.22, 2.22, 3.22),
              annotation=p.stats, tip_length=0) 
dev.off()
################################################################################
# Figure 5h. LocusZoom plot of mash effect sizes

# LDlink token
t <- '' # requires personal access token
# load Ensembl database
ah <- AnnotationHub()
ensembl.v105 <- ah[['AH98047']]
# create LD lookup function
mem_LDmatrix <- memoise(LDlinkR::LDmatrix)
get_LD <- function(l, coord = 'id', index_snp = '', pop = 'EUR', genome_build = 'grch38', token = '', ...){
  labs <- l$labs
  index_snp <- index_snp
  r <- l$data[which(l$data[, labs] == index_snp), coord]
  rslist <- l$data[, coord]
  rslist <- unique(c(rslist, r))
  ldm <- mem_LDmatrix(rslist, pop = pop, genome_build = genome_build, token = token, ...)
  ld <- ldm[, index_snp]
  l$data$ld <- ld[match(l$data[, labs], ldm$RS_number)]
  l
}
get_LD_long <- function(l, coord = 'id', index_snp = '', pop = 'EUR', genome_build = 'grch38', token = '', ...){
  labs <- l$labs
  index_snp <- index_snp
  r <- l$data[which(l$data[, labs] == index_snp), coord]
  rslist <- l$data[, coord]
  if (length(rslist > 1000)) {
    rslist <- rslist[order(l$data$logP, decreasing = TRUE)[seq_len(1000)]]
  }
  rslist <- unique(c(rslist, r))
  ldm <- mem_LDmatrix(rslist, pop = pop, genome_build = genome_build, token = token, ...)
  ld <- ldm[, index_snp]
  l$data$ld <- ld[match(l$data[, labs], ldm$RS_number)]
  l
}
# load coloc SNPs
df.gwas <- load('data/UNC/RData/locusZoom_matrices_GWAS-eQTL_coloc_H4.RData')
# reformat
df.gwas <- df.gwas %>% unite('id', chr:pos, sep = ':', remove = FALSE)
df.gwas$chr <- gsub('chr', '', df.gwas$chr)
# select candidate genes
df.candidates <- data.frame(gene_symb = c('ABO', 'APOBR', 'CTSS', 'FLRT3', 'FUT2', 'NXPE1', 'TMEM170A', 'TNFRSF14', 'TUFM', 'ZFP90'),
                            gene_id = c('ENSG00000175164','ENSG00000184730','ENSG00000163131','ENSG00000125848','ENSG00000176920','ENSG00000095110','ENSG00000166822','ENSG00000157873','ENSG00000178952','ENSG00000184939'),
                            win_size = c(5e4, 1e5, 9e4, 7e5, 1e5, 1e5, 3e5, 1e5, 1e6, 2e5),
                            center = c('ABO', 'APOBR', 'CTSS', 'MACROD2-IT1', 'FUT2', 'NXPE1', 'CFDP1', 'TNFRSF14', 'TUFM', 'ZFP90'))
# pull out colocalizing gene info
cb <- coloc.barcuva %>% select(gene, gene_name, GWAS_index_snp)
cg <- coloc.gtex %>% select(gene, gene_name, GWAS_index_snp)
ci <- coloc.ibd %>% select(gene, name, GWAS_index_snp)
# merge GWAS results
genes <- merge(ci, cg, by.x = colnames(ci), by.y = colnames(cg), suffixes = c('_IBD', '_GTEx'), all = TRUE) %>%
  distinct(gene, .keep_all = TRUE)
genes <- merge(genes, cb, by.x = colnames(genes), by.y = colnames(cb), all = TRUE) %>%
  distinct(gene, .keep_all = TRUE)
# loop through each gene
for (i in 1:nrow(df.candidates)) {
  # select gene
  gene.id <- df.candidates$gene_id[i]
  gene.sym <- df.candidates$gene_symb[i]
  center <- df.candidates$center[i]
  print(c(i, gene.sym))
  # set window size
  win.size <- df.candidates$win_size[i]
  # subset mash results
  m.g <- subset(mc.g, gene == gene.id) %>% distinct(eqtl, .keep_all = TRUE)
  m.i <- subset(mc.i, gene == gene.id) %>% distinct(eqtl, .keep_all = TRUE)
  snps <- intersect(m.g$snp, m.i$snp)
  if (length(snps) == 0) {
    next
  }
  # subset variants for locus
  df.locus <- subset(m.g, snp %in% snps) %>% arrange(match(snp, snps))
  df.locus <- df.locus %>% separate(snp, c('chr', 'pos', 'ref', 'alt'), remove = FALSE) %>% unite('id', chr:pos, sep = ':', remove = FALSE)
  df.locus$chr <- gsub('chr', '', df.locus$chr)
  df.locus$abs.gtex <- abs(df.locus$gtex)
  df.locus$abs.ibd <- abs(df.locus$ibd)
  # add standard errors
  sd <- subset(coloc.sd, eqtl %in% df.locus$eqtl) %>% distinct(eqtl, .keep_all = TRUE) %>% arrange(match(eqtl, df.locus$eqtl))
  # subset GWAS
  gwas.gene <- df.gwas[df.gwas$gene == gene.sym,]
  if (nrow(gwas.gene) == 0) {
    next
  }
  # set locus
  loc.gwas <- locus(gwas.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'P-value', 
                    gene = center, fix_window = win.size)
  loc.gtex <- locus(df.locus, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', yvar = 'abs.gtex', 
                    gene = center, fix_window = win.size)
  loc.ibd <- locus(df.locus, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', yvar = 'abs.ibd', 
                   gene = center, fix_window = win.size)
  # set index snp
  index_snp <- df.locus[df.locus$snp == subset(genes, gene == gene.id)$GWAS_index_snp,]$rsID
  x <- loc.ibd$data[which(loc.ibd$data$rsID == index_snp),1:6]
  # add LD info
  if (nrow(loc.gwas$data) > 999) {
    loc.gwas <- get_LD_long(loc.gwas, token = t, index_snp = index_snp)
  } else {
    loc.gwas <- get_LD(loc.gwas, token = t, index_snp = index_snp)
  }
  
  if (nrow(loc.gtex$data) > 999) {
    loc.gtex <- get_LD_long(loc.gtex, token = t, index_snp = index_snp)
  } else {
    loc.gtex <- get_LD(loc.gtex, token = t, index_snp = index_snp)
  }
  
  if (nrow(loc.ibd$data) > 999) {
    loc.ibd <- get_LD_long(loc.ibd, token = t, index_snp = index_snp)
  } else {
    loc.ibd <- get_LD(loc.ibd, token = t, index_snp = index_snp)
  }

  # plot effect sizes on same panel
  file = paste('plots/LocusZoom_plot_mash_effect_sizes_', gene.sym, '_eQTL.png', sep = '')
  png(filename = file, res = 300, units = 'in', height = 12, width = 6)
  oldpar <- set_layers(2)
  try(scatter_plot(loc.gwas, index_snp = index_snp, labels = 'index',
                   border = TRUE, pcutoff = NULL, cex = 3), silent = TRUE)
  title(gene.sym, adj = 0, line = 2.55)
  title('GWAS', adj = 0.2, line = 1)
  try(scatter_plot(loc.ibd, index_snp = index_snp, labels = 'index',
                   border = TRUE, legend_pos = NULL, ylab = 'Absolute mash-adjusted effect size', 
                   cex = 3, ylim = c(0, round(max(max(df.locus$gtex, df.locus$ibd)), 1)+0.1)), 
      silent = TRUE)
  title('mash-adjusted effect sizes', adj = 0, line = 1)
  try(scatter_plot(loc.gtex, index_snp = index_snp, labels = 'index', 
                   border = TRUE, legend_pos = NULL, ylab = 'Absolute mash-adjusted effect size', 
                   cex = 3, ylim = c(0, round(max(max(df.locus$gtex, df.locus$ibd)), 1)+0.1),
                   pch = 24, add = TRUE), silent = TRUE)
  legend(x = 'topright', legend=c("UNC", "GTEx"), pch = c(16,17), cex = 1.5, pt.cex = 3, 
         col = c('black', 'black'), bty = 'n')
  genetracks(loc.ibd)
  dev.off()
}

################################################################################
