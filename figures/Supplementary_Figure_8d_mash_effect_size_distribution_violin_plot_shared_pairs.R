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
# subset mash results
gtex.all.shared <- subset(mc.g, gene %in% shared.genes)
ibd.all.shared <- subset(mc.i, gene %in% shared.genes)
# subset to only shared eSNP-eGene pairs for a sanity check
shared.eqtl <- intersect(gtex.all.shared$eqtl, ibd.all.shared$eqtl)
shared.eqtl <- subset(gtex.all.shared, eqtl %in% shared.eqtl)
median(abs(shared.eqtl$gtex))
median(abs(shared.eqtl$ibd))
wilcox.test(abs(shared.eqtl$gtex),abs(shared.eqtl$ibd))
# plot
si <- shared.eqtl %>% select(eqtl, ibd)
colnames(si) <- c('eqtl', 'slope')
shared.eqtl$Group <- 'GTEx'
si$Group <- 'UNC'
shared.eqtl <- full_join(shared.eqtl, si)
# plot
png(filename = 'plots/Supplementary_Figure_8d_violin_plot_mash_shared_eSNP-eGene_pairs_GTEx_UNC.png',
    res = 300, units = 'in', height = 5, width = 6)
ggplot(shared.eqtl, aes(Group, abs(slope), fill = Group)) + geom_violin(lwd=1) + 
  geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) + 
  labs(x='', y='mash-Adjusted Absoulte eQTL Effect Size', title='Shared eVariant-eGene pairs that colocalize with GWAS') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 16)) + 
  geom_signif(test = 'wilcox.test', comparisons = list(c('GTEx', 'UNC'))) 
dev.off()
