# Figure 4a. violin plot for colocalizing eGene TSS distances
library(data.table)
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)
library(qvalue)
library(ggplot2)
library(ggsignif)

coloc.ibd <- fread('data/UNC/Supplementary_Table_7_coloc_h4_IBD_eQTL_GWAS.txt')
coloc.gtex <- fread('data/GTEx/Supplementary_Table_8_coloc_H4_GTEx_eQTL_colon_transverse_GWAS.txt')
coloc.barcuva <- fread('data/BarcUVa/Supplementary_Table_9_coloc_H4_BarcUVa_eQTL_GWAS.txt')
# reformat
coloc.ibd <- coloc.ibd %>% separate(gene, c('gene', 'ver'))
coloc.barcuva <- coloc.barcuva %>% separate(gene, c('gene', 'ver'))
coloc.gtex <- coloc.gtex %>% separate(gene, c('gene', 'ver'))
################################################################################
# UpSet plot for colocalizing eGenes
set.seed(12345)
# create a list of colocalizing eGenes
set.egenes <- list(UNC = coloc.ibd$gene, 
                   GTEx = coloc.gtex$gene,
                   BarcUVa = coloc.barcuva$gene)
# create combination matrix
mat.egenes <- make_comb_mat(set.egenes, mode = 'distinct')
################################################################################
# subset unique egenes
barcuva.unique <- subset(coloc.barcuva, gene %in% extract_comb(mat.egenes, '001'))
gtex.unique <- subset(coloc.gtex, gene %in% extract_comb(mat.egenes, '010'))
ibd.unique <- subset(coloc.ibd, gene %in% extract_comb(mat.egenes, '100'))
# subset shared egenes
shared.egenes <- extract_comb(mat.egenes, '111')
ibd.shared <- subset(coloc.ibd, gene %in% shared.egenes)
gtex.shared <- subset(coloc.gtex, gene %in% shared.egenes)
barcuva.shared <- subset(coloc.barcuva, gene %in% shared.egenes)
# calculate median distance to TSS
median(abs(barcuva.shared$tss_distance))
median(abs(gtex.shared$tss_distance))
median(abs(ibd.shared$dist_phe_var))
median(c(abs(barcuva.shared$tss_distance),abs(gtex.shared$tss_distance),abs(ibd.shared$dist_phe_var)))
median(c(abs(gtex.shared$tss_distance),abs(ibd.shared$dist_phe_var)))
median(abs(barcuva.unique$tss_distance))
median(abs(gtex.unique$tss_distance))
median(abs(ibd.unique$dist_phe_var))
# test for difference
wb <- wilcox.test(abs(barcuva.shared$tss_distance), abs(barcuva.unique$tss_distance))
wg <- wilcox.test(abs(gtex.shared$tss_distance), abs(gtex.unique$tss_distance))
wi <- wilcox.test(abs(ibd.shared$dist_phe_var), abs(ibd.unique$dist_phe_var))
p.stats <- c(round(wb$p.value, 4), round(wg$p.value, 3), round(wi$p.value, 4))
p.stats <- gsub('^', 'p = ', p.stats)
# format for plotting
b.s <- barcuva.shared %>% select(eqtl, tss_distance)
g.s <- gtex.shared %>% select(eqtl, tss_distance)
i.s <- ibd.shared %>% select(eqtl, dist_phe_var)
colnames(i.s) <- c('eqtl', 'tss_distance')
b.s$group <- 'BarcUVa'
g.s$group <- 'GTEx'
i.s$group <- 'UNC'
shared <- full_join(b.s, g.s)
shared <- full_join(shared, i.s)
shared$eGene <- 'shared'
b.u <- barcuva.unique %>% select(eqtl, tss_distance)
g.u <- gtex.unique %>% select(eqtl, tss_distance)
i.u <- ibd.unique %>% select(eqtl, dist_phe_var)
colnames(i.u) <- c('eqtl', 'tss_distance')
b.u$group <- 'BarcUVa'
g.u$group <- 'GTEx'
i.u$group <- 'UNC'
uniq <- full_join(b.u, g.u)
uniq <- full_join(uniq, i.u)
uniq$eGene <- 'unique'
p <- full_join(shared, uniq)
# summary stats
ss <- p %>% group_by(group, eGene) %>% summarize(n=n(), med=median(abs(tss_distance)), 
      yn = -25000, ym = -10000)
ss <- ss %>% mutate(med = round(med/1000, 2))
# plot
png(res = 300, units = 'in', height = 8, width = 8,
    filename = 'plots/Figure_4a_Violin_plot_eGene_TSS_distance.png')
ggplot(p, aes(group, abs(tss_distance), fill = eGene)) + geom_violin(lwd=1) +
  geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) +
  labs(x='', y='Absoulte Distance to eGene TSS') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  geom_signif(y_position=c(1e6, 8.5e5, 8.5e5), xmin=c(0.77, 1.77, 2.77), xmax=c(1.22, 2.22, 3.22),
              annotation=p.stats, tip_length=0) + 
  geom_text(data = ss, aes(group, y = ym, label = paste('med =', med, 'kb')), 
            position = position_dodge(width = 0.9), size = 3) +
  geom_text(dat = ss, aes(group, y = yn, label = paste0('n = ', n)), 
            position = position_dodge(width = 0.9), size = 3)
dev.off()
