# Figure 4c. mash effect size distribution for colocalizing eGenes
library(data.table)
library(tidyverse)
library(mashr)
library(ggplot2)
library(ggsignif)

# load mash data
load(file = 'data/UNC/RData/mashr_coloc_fit_full_signal.RData')
# drop alternate rsIDs
mc.b <- mc.b %>% distinct(eqtl, .keep_all = TRUE)
mc.g <- mc.g %>% distinct(eqtl, .keep_all = TRUE)
mc.i <- mc.i %>% distinct(eqtl, .keep_all = TRUE)
# create gene sets
shared.genes <- Reduce(intersect, list(coloc.barcuva$gene, coloc.gtex$gene, coloc.ibd$gene))
gtex.genes <- setdiff(coloc.gtex$gene, c(coloc.barcuva$gene, coloc.ibd$gene))
ibd.genes <- setdiff(coloc.ibd$gene, c(coloc.barcuva$gene, coloc.gtex$gene))
# subset mash results
gtex.all.shared <- subset(mc.g, gene %in% shared.genes)
gtex.all.unique <- subset(mc.g, gene %in% gtex.genes)
ibd.all.shared <- subset(mc.i, gene %in% shared.genes)
ibd.all.unique <- subset(mc.i, gene %in% ibd.genes)
# stats
ws <- wilcox.test(abs(gtex.all.shared$gtex), abs(ibd.all.shared$ibd))
wu <- wilcox.test(abs(gtex.all.unique$gtex), abs(ibd.all.unique$ibd))
median(abs(gtex.all.shared$gtex))
median(abs(ibd.all.shared$ibd))
median(abs(gtex.all.unique$gtex))
median(abs(ibd.all.unique$ibd))
wilcox.test(abs(gtex.all.shared$gtex), abs(gtex.all.unique$gtex))
wilcox.test(abs(ibd.all.shared$ibd), abs(ibd.all.unique$ibd))
p.stats <- c(paste('p =', round(ws$p.value, 3)), 'p < 2.2e-16')
# reformat for plotting
gtex.all.shared$eGene <- 'Shared eGene'
gtex.all.unique$eGene <- 'Unique eGene'
gtex.all.shared$slope <- gtex.all.shared$gtex
gtex.all.unique$slope <- gtex.all.unique$gtex
g <- full_join(gtex.all.shared, gtex.all.unique)
g$Group <- 'GTEx'
ibd.all.shared$eGene <- 'Shared eGene'
ibd.all.unique$eGene <- 'Unique eGene'
ibd.all.shared$slope <- ibd.all.shared$ibd
ibd.all.unique$slope <- ibd.all.unique$ibd
i <- full_join(ibd.all.shared, ibd.all.unique)
i$Group <- 'UNC'
p <- full_join(g,i)
# convert to log10 scale
p <- p %>% mutate(log10.abs = log10(abs(slope)))
# summary stats
ss <- p %>% group_by(eGene, Group) %>% summarize(n=n(), med=round(median(abs(slope)), 3),
      yn=-0.05, ym=-0.02)
ss <- p %>% group_by(eGene, Group) %>% summarize(n=n(), med=round(median(log10.abs), 3),
                                                 yn=-1.67, ym=-1.6)
# violin plot
png(filename = 'plots/Figure_4c_violin_plot_mash_shared_vs_unique_eGenes_GTEx_UNC_all_eQTL.png',
    res = 300, units = 'in', height = 5, width = 6)
  ggplot(p, aes(eGene, abs(slope), fill = Group)) + geom_violin(lwd=1) + 
    geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) + 
    labs(x='', y='mash-Adjusted Absoulte eQTL Effect Size', title='eQTL that colocalize with GWAS') +
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
          plot.title = element_text(hjust = 0.5, size = 16)) + 
  geom_signif(y_position=c(1.07, 1.16), xmin=c(0.78, 1.78), xmax=c(1.22, 2.22),
              annotation=p.stats, tip_length=0) +
    geom_text(data = ss, aes(eGene, y=ym, label = paste0('med = ', med)),
              position = position_dodge(width = 0.9), size = 3) +
    geom_text(data = ss, aes(eGene, y=yn, label = paste0('n = ', n)), 
              position = position_dodge(width = 0.9), size = 3)
dev.off()

png(filename = 'plots/Figure_4c_violin_plot_mash_shared_vs_unique_eGenes_GTEx_UNC_all_eQTL_log10.png',
    res = 300, units = 'in', height = 5, width = 6)
ggplot(p, aes(eGene, log10.abs, fill = Group)) + geom_violin(lwd=1) + 
  geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) + 
  labs(x='', y='log10(absolute mash-adjusted eQTL effect size)', title='eQTL that colocalize with GWAS') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 16)) + 
  geom_signif(y_position=c(0.1, 0.1), xmin=c(0.78, 1.78), xmax=c(1.22, 2.22),
              annotation=p.stats, tip_length=0) +
  geom_text(data = ss, aes(eGene, y=ym, label = paste0('med = ', med)),
            position = position_dodge(width = 0.9), size = 4) +
  geom_text(data = ss, aes(eGene, y=yn, label = paste0('n = ', n)), 
            position = position_dodge(width = 0.9), size = 4)
dev.off()
