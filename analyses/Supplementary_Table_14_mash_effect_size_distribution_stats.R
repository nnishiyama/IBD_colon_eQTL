# Supplementary Table 14. Compare mash effect size distributions for colocalizing eQTL
library(data.table)
library(tidyverse)
library(mashr)
library(qvalue)
library(ggplot2)
library(ggsignif)

# load mash data
load(file = '/work/users/n/n/nnishi/mashr/mashr_coloc_fit_full_signal.RData')
# pull out colocalizing gene info
genes <- c(shared.genes, gtex.genes, ibd.genes)
cg <- subset(coloc.gtex, gene %in% genes) %>% select(gene, gene_name, GWAS_index_snp)
ci <- subset(coloc.ibd, gene %in% genes) %>% select(gene, name, GWAS_index_snp)
# merge
genes <- merge(ci, cg, by.x = colnames(ci), by.y = colnames(cg), suffixes = c('_IBD', '_GTEx'), all = TRUE) %>%
  distinct(gene, .keep_all = TRUE)
# loop through and compare distributions for each colocalizing eGene
df.all <- data.frame()
df.overlap <- data.frame()
for (x in 1:nrow(genes)) {
  print(x)
  ## first: using all passing SNPs
  egene <- genes$gene[x]
  g <- subset(mc.g, gene == egene) %>% distinct(eqtl, .keep_all = TRUE)
  i <- subset(mc.i, gene == egene) %>% distinct(eqtl, .keep_all = TRUE)
  # annotate number of snps
  num.g <- length(unique(g$snp))
  num.i <- length(unique(i$snp))
  # calculate median mash-adjusted effect size
  med.g <- median(abs(g$gtex))
  med.i <- median(abs(i$ibd))
  # annotate data set with larger median effect size
  if (is.na(med.g) | is.na(med.i)) {
    dist <- NA
  } else if (med.g > med.i) {
    dist <- 'GTEx'
  } else dist <- 'UNC'
  # check for no significant results
  if (num.g == 0 | num.i == 0) {
    print('No significant results for at least one data set')
    # append & move on
    out <- c(egene, genes$name[x], bin, genes$GWAS_index_snp[x], num.g, num.i, med.g, med.i, dist, NA, NA)
    df.all <- rbind(df.all, out)
    out <- c(egene, genes$name[x], bin, genes$GWAS_index_snp[x], 0, NA, NA, dist, NA, NA, NA, NA)
    df.overlap <- rbind(df.overlap, out)
    next
  }
  # statistic to compare distribution
  stat <- wilcox.test(abs(g$gtex), abs(i$ibd))
  w <- stat$statistic
  w.p <- stat$p.value
  # bin gene
  if (egene %in% g.genes) {
    bin <- 'GTEx only'
  } else if (egene %in% i.genes) {
    bin <- 'UNC only'
  } else bin <- 'shared'
  # append
  out <- c(egene, genes$name[x], bin, genes$GWAS_index_snp[x], num.g, num.i, med.g, med.i, dist, w, w.p)
  df.all <- rbind(df.all, out)
  ## second: using only shared SNPs
  snps <- intersect(g$snp, i$snp)
  df <- subset(g, snp %in% snps) %>% arrange(match(snp, snps))
  # annotate number of overlapping snps
  num.snps <- length(snps)
  # calculate median mash-adjusted effect size
  med.g <- median(abs(df$gtex))
  med.i <- median(abs(df$ibd))
  # annotate data set with larger median effect size
  if (is.na(med.g) | is.na(med.i)) {
    dist <- NA
  } else if (med.g > med.i) {
    dist <- 'GTEx'
  } else dist <- 'UNC'
  # statistic to compare distribution
  stat <- wilcox.test(abs(df$gtex), abs(df$ibd))
  w <- stat$statistic
  w.p <- stat$p.value
  # check number of snps
  if (num.snps < 3) {
    # append & move on
    print('Not enough snps')
    out <- c(egene, genes$name[x], bin, genes$GWAS_index_snp[x], num.snps, med.g, med.i, dist, w, w.p, NA, NA)
    df.overlap <- rbind(df.overlap, out)
    next
  }
  # correlation
  corr <- cor.test(df$gtex, df$ibd)
  r <- corr$estimate
  r.p <- corr$p.value
  # append
  out <- c(egene, genes$name[x], bin, genes$GWAS_index_snp[x], num.snps, med.g, med.i, dist, w, w.p, r, r.p)
  df.overlap <- rbind(df.overlap, out)
  # plot
  gg <- data.frame(eqtl = c(df$eqtl, df$eqtl), slope = c(df$gtex, df$ibd), group = c(rep('GTEx', nrow(df)), rep('UNC', nrow(df))))
  # summary stats for plotting
  ss <- gg %>% group_by(group) %>% summarize(n=n(), med=round(median(abs(slope)), 3),
        yn=0.025, ym=0.03)
  # violin plot
  file <- paste('plots/mashr_violin_plot_coloc_full_signal_', genes$name[x], '.png', sep = '')
  p <- ggplot(gg, aes(group, abs(slope), fill = group)) + geom_violin(lwd=1) +
    geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) +
    labs(x='', y='mash-Adjusted Absoulte eQTL Effect Size', title=paste(genes$name[x], 'eQTL signal')) +
    theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
          plot.title = element_text(hjust = 0.5, size = 16), legend.position = 'none') +
    geom_signif(test = 'wilcox.test', comparisons = list(c('GTEx', 'UNC'))) +
    geom_text(data = ss, aes(group, y=ym, label = paste0('med = ', med)), size = 3) +
    geom_text(data = ss, aes(group, y=yn, label = paste0('n = ', n)), size = 3)
  ggsave(p, file=file, width = 6, height = 6, units = "in", dpi=300)
  # correlation plot
  file <- paste('plots/mashr_correlation_plot_coloc_full_signal_', genes$name[x], '.png', sep = '')
  png(filename = file, res = 300, units = 'in', height = 6, width = 6)
  plot(df$ibd, df$gtex, xlim = c(-1,1), ylim = c(-1,1),
        xlab = 'UNC mash-adjusted effect sizes', ylab = 'GTEx mash-adjusted effect sizes',
        main = paste(genes$name[x], 'eQTL signal'))
  abline(0,1, col='red')
  abline(v=0,h=0)
  dev.off()
}
# label columns
colnames(df.all) <- c('gene', 'name', 'gene_sharing_bin', 'GWAS_index_snp',
                      'number_variants_GTEx', 'number_variants_IBD', 
                      'abs_median_effect_size_GTEx', 'abs_median_effect_size_UNC',
                      'larger_distribution','wilcoxon_statistic', 'wilcoxon_pvalue')
colnames(df.overlap) <- c('gene', 'name', 'gene_sharing_bin', 'GWAS_index_snp', 'number_variants_overlapping',
                      'abs_median_effect_size_GTEx', 'abs_median_effect_size_UNC', 'larger_distribution',
                      'wilcoxon_statistic', 'wilcoxon_pvalue', 'pearson_r', 'pearson_pvalue')
# convert to numeric
df.all[names(df.all)[c(5:8,10:11)]] <- lapply(df.all[names(df.all)[c(5:8,10:11)]], as.numeric)
df.overlap[names(df.overlap)[c(5:7,9:12)]] <- lapply(df.overlap[names(df.overlap)[c(5:7,9:12)]], as.numeric)
# multiple hypothesis testing correction
q <- qvalue(df.all$wilcoxon_pvalue, pi0 = 1)
df.all <- df.all %>% mutate(wilcoxon_padj = q$qvalues, .after = wilcoxon_pvalue)
q <- qvalue(df.overlap$wilcoxon_pvalue, pi0 = 1)
df.overlap <- df.overlap %>% mutate(wilcoxon_padj = q$qvalues, .after = wilcoxon_pvalue)
q <- qvalue(df.overlap$pearson_pvalue, pi0 = 1)
df.overlap <- df.overlap %>% mutate(pearson_padj = q$qvalues, .after = pearson_pvalue)
# sort
df.overlap <- df.overlap %>% arrange(desc(larger_distribution), desc(gene_sharing_bin))
################################################################################
# write out results
write.table(df.overlap, quote = FALSE, row.names = FALSE, sep = '\t',
            file = 'data/UNC/Supplementary_Table_14_mashr_coloc_fit_full_signal_effect_size_distribution_stats.txt')
################################################################################
# filter
df.filtered <- df.overlap[df.overlap$wilcoxon_padj < 0.05,] %>% drop_na(wilcoxon_padj)
table(df.filtered$gene_sharing_bin, df.filtered$larger_distribution)
