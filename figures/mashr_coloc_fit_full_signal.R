# mash analysis: pull out whole signal to compare distributions
library(data.table)
library(tidyverse)
library(mashr)
library(qvalue)
library(ggplot2)
library(ggsignif)
################################################################################
# load summary statistics
load('/work/users/n/n/nnishi/eqtl/freeze/final/RData/full_summary_stats_rsIDs.RData')
# drop GWAS data
remove(cd, ibd, uc)
# load coloc results
coloc.barcuva <- fread('/work/users/n/n/nnishi/gwas/liu/coloc_H4_BarcUVa_eQTL_colon_GWAS_Liu_all_pop_20240111.txt')
coloc.gtex <- fread('/work/users/n/n/nnishi/gwas/liu/coloc_H4_GTEx_eQTL_colon_transverse_GWAS_Liu_all_pop_20240111.txt')
coloc.ibd <- fread('/work/users/n/n/nnishi/gwas/liu/coloc_h4_IBD_MAF002_GWAS_Liu_all_pop_20240111.txt')
# drop NAs
coloc.barcuva <- coloc.barcuva %>% drop_na(GWAS_index_snp)
# reformat
coloc.barcuva <- coloc.barcuva %>% separate(gene, c('gene', 'ver'))
coloc.gtex <- coloc.gtex %>% separate(gene, c('gene', 'ver'))
coloc.ibd <- coloc.ibd %>% separate(gene, c('gene', 'ver'))
# create a union list of colocalizing eGenes
genes <- unique(c(coloc.barcuva$gene, coloc.gtex$gene, coloc.ibd$gene))
# drop BarcUVa-unique eGenes
genes <- genes[! genes %in% setdiff(coloc.barcuva$gene, c(coloc.gtex$gene, coloc.ibd$gene))]
# pull out colocalizing gene info
cg <- subset(coloc.gtex, gene %in% genes) %>% select(gene, gene_name, GWAS_index_snp)
ci <- subset(coloc.ibd, gene %in% genes) %>% select(gene, name, GWAS_index_snp)
# merge
genes <- merge(ci, cg, by.x = colnames(ci), by.y = colnames(cg), suffixes = c('_IBD', '_GTEx'), all = TRUE) %>%
  distinct(gene, .keep_all = TRUE)
# create categorized gene lists
shared.genes <- unique(c(intersect(coloc.ibd$gene, coloc.barcuva$gene), intersect(coloc.ibd$gene, coloc.gtex$gene)))
g.genes <- setdiff(coloc.gtex$gene, shared.genes) # for now keeping genes only shared between gtex & barcuva
i.genes <- setdiff(coloc.ibd$gene, shared.genes)
# reformat ENSEMBL IDs
eqtl.barcuva <- eqtl.barcuva %>% separate(gene_id, c('gene', 'ver'))
eqtl.gtex <- eqtl.gtex %>% separate(gene_id, c('gene', 'ver'))
eqtl.ibd <- eqtl.ibd %>% separate(phe_id, c('gene', 'ver'))
# annotate eqtl
eqtl.barcuva$eqtl <- paste(eqtl.barcuva$gene, eqtl.barcuva$snp, sep = '-')
eqtl.gtex$eqtl <- paste(eqtl.gtex$gene, eqtl.gtex$snp, sep = '-')
eqtl.ibd$eqtl <- paste(eqtl.ibd$gene, eqtl.ibd$snp, sep = '-')
# subset data to gene list
eb <- subset(eqtl.barcuva, gene %in% genes$gene) %>% select(eqtl, slope, slope_se, rsID)
eg <- subset(eqtl.gtex, gene %in% genes$gene) %>% select(eqtl, slope, slope_se, rsID)
ei <- subset(eqtl.ibd, gene %in% genes$gene) %>% select(eqtl, slope, slope_se, rsID)
colnames(ei) <- c('eqtl', 'slope_IBD', 'slope_se_IBD', 'rsID')
# drop unnecessary data
remove(eqtl.barcuva, eqtl.gtex, eqtl.ibd, cg, ci)
# merge
df.coloc <- merge(eb, eg, by = c('eqtl', 'rsID'), all = TRUE, suffixes = c('_BarcUVa', '_GTEx'))
df.coloc <- merge(df.coloc, ei, by = c('eqtl', 'rsID'), all = TRUE)
# reformat
bhat <- cbind(df.coloc$slope_BarcUVa, df.coloc$slope_GTEx, df.coloc$slope_IBD)
colnames(bhat) <- c('barcuva', 'gtex', 'ibd')
row.names(bhat) <- df.coloc$eqtl
sehat <- cbind(df.coloc$slope_se_BarcUVa, df.coloc$slope_se_GTEx, df.coloc$slope_se_IBD)
colnames(sehat) <- c('barcuva', 'gtex', 'ibd')
row.names(sehat) <- df.coloc$eqtl
dof <- cbind(rep(443, nrow(df.coloc)), rep(366, nrow(df.coloc)), rep(250, nrow(df.coloc)))
colnames(dof) <- c('barcuva', 'gtex', 'ibd')
row.names(dof) <- df.coloc$eqtl
# load mash data
load('/work/users/n/n/nnishi/mashr/mashr_colon_eqtl_random_null_corr.RData')
# create a mash object
data.coloc <- mash_set_data(Bhat = bhat, Shat = sehat, df = dof, V = data.strong$V)
# apply fit
m.coloc <- mash(data.coloc, g = get_fitted_g(m), fixg = TRUE)
# total number of effects that mash considers significant
ic <- get_significant_results(m.coloc)
# significant results by condition
length(get_significant_results(m.coloc, conditions = 1))
length(get_significant_results(m.coloc, conditions = 2))
length(get_significant_results(m.coloc, conditions = 3))
# calculate the proportion of pairwise sharing between studies
get_pairwise_sharing(m.coloc)
# extract posterior matrices & subset to mash significant eQTL
coloc.beta <- as.data.frame(m.coloc$result$PosteriorMean)
coloc.sd <- as.data.frame(m.coloc$result$PosteriorSD)
coloc.lfsr <- as.data.frame(m.coloc$result$lfsr)
# reformat
coloc.beta$rsID <- df.coloc$rsID
coloc.beta$eqtl <- df.coloc$eqtl
coloc.beta <- coloc.beta %>% separate(eqtl, c('gene', 'snp'), sep = '-', remove = FALSE)
# pull out significant results for each condition
mc.b <- coloc.beta[get_significant_results(m.coloc, conditions =1),]
mc.g <- coloc.beta[get_significant_results(m.coloc, conditions =2),]
mc.i <- coloc.beta[get_significant_results(m.coloc, conditions =3),]
# save workspace
save.image(file = '/work/users/n/n/nnishi/mashr/mashr_coloc_fit_full_signal.RData')
# format for writing output
colnames(coloc.lfsr) <- paste0(colnames(coloc.lfsr), '_lfsr')
coloc.lfsr$eqtl <- coloc.beta$eqtl
coloc.sd$eqtl <- coloc.beta$eqtl
out <- merge(coloc.beta, coloc.sd, by = 'eqtl', suffixes = c('_beta', '_sd')) %>% select(-c(gene:snp)) #%>% relocate(gene:snp, .after = eqtl)
out <- merge(out, coloc.lfsr, by = 'eqtl') %>% distinct()
# write out
write.table(out, sep = '\t', quote = FALSE, row.names = FALSE, 
            file = '/work/users/n/n/nnishi/mashr/mashr_coloc_fit_full_signal_colocalizing_eQTL.txt')
################################################################################
# focus only on GTEx & UNC eGenes
load(file = '/work/users/n/n/nnishi/mashr/mashr_coloc_fit_full_signal.RData')
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
# summary stats
ss <- p %>% group_by(eGene, Group) %>% summarize(n=n(), med=round(median(abs(slope)), 3),
      yn=-0.05, ym=-0.02)
# violin plot
png(filename = '/work/users/n/n/nnishi/eqtl/freeze/final/figures/Figure4b_violin_plot_mash_shared_vs_unique_eGenes_GTEx_UNC_all_eQTL.png',
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
# subset to only shared eSNP-eGene pairs for a sanity check
shared.eqtl <- intersect(gtex.all.shared$eqtl, ibd.all.shared$eqtl)
shared.eqtl <- subset(gtex.all.shared, eqtl %in% shared.eqtl)
median(abs(shared.eqtl$gtex))
median(abs(shared.eqtl$ibd))
wilcox.test(abs(shared.eqtl$gtex),abs(shared.eqtl$ibd))
# calculate difference
shared.eqtl <- shared.eqtl %>% mutate(diff = abs(ibd) - abs(gtex))
median(abs(shared.eqtl$diff))
hist(abs(shared.eqtl$diff))
hist(shared.eqtl$diff)
range(abs(shared.eqtl$diff))
cor.test(abs(shared.eqtl$gtex), abs(shared.eqtl$ibd))
plot(shared.eqtl$gtex, shared.eqtl$ibd)
abline(v=0, h=0)
abline(0,1, col = 'red')
# annotate direction
shared.eqtl <- shared.eqtl %>% mutate(direction = ifelse(gtex > 0 & ibd > 0, 'same',
                                                  ifelse(gtex < 0 & ibd < 0, 'same', 'opposite')))
# plot
si <- shared.eqtl %>% select(eqtl, ibd)
colnames(si) <- c('eqtl', 'slope')
shared.eqtl$Group <- 'GTEx'
si$Group <- 'UNC'
shared.eqtl <- full_join(shared.eqtl, si)
# plot
png(filename = '/work/users/n/n/nnishi/eqtl/freeze/final/figures/violin_plot_mash_shared_eSNP-eGene_pairs_GTEx_UNC.png',
    res = 300, units = 'in', height = 5, width = 6)
ggplot(shared.eqtl, aes(Group, abs(slope), fill = Group)) + geom_violin(lwd=1) + 
  geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) + 
  labs(x='', y='mash-Adjusted Absoulte eQTL Effect Size', title='Shared eVariant-eGene pairs that colocalize with GWAS') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 16)) + 
  geom_signif(test = 'wilcox.test', comparisons = list(c('GTEx', 'UNC'))) 
dev.off()
################################################################################
# top variants only

# subset mash results
gtex.top.shared <- subset(mc.g, eqtl %in% subset(coloc.gtex, gene %in% shared.genes)$eqtl)
gtex.top.unique <- subset(mc.g, eqtl %in% subset(coloc.gtex, gene %in% gtex.genes)$eqtl)
ibd.top.shared <- subset(mc.i, eqtl %in% subset(coloc.ibd, gene %in% shared.genes)$eqtl)
ibd.top.unique <- subset(mc.i, eqtl %in% subset(coloc.ibd, gene %in% ibd.genes)$eqtl)
# stats 
ws.top <- wilcox.test(abs(gtex.top.shared$gtex), abs(ibd.top.shared$ibd))
wu.top <- wilcox.test(abs(gtex.top.unique$gtex), abs(ibd.top.unique$ibd))
median(abs(gtex.top.shared$gtex))
median(abs(ibd.top.shared$ibd))
median(abs(gtex.top.unique$gtex))
median(abs(ibd.top.unique$ibd))
p.stats <- c(round(ws.top$p.value, 3), round(wu.top$p.value, 3))
p.stats <- gsub('^', 'p = ', p.stats)
# reformat for plotting
gtex.top.shared$eGene <- 'Shared eGene'
gtex.top.unique$eGene <- 'Unique eGene'
gtex.top.shared$slope <- gtex.top.shared$gtex
gtex.top.unique$slope <- gtex.top.unique$gtex
g.top <- full_join(gtex.top.shared, gtex.top.unique)
g.top$Group <- 'GTEx'
ibd.top.shared$eGene <- 'Shared eGene'
ibd.top.unique$eGene <- 'Unique eGene'
ibd.top.shared$slope <- ibd.top.shared$ibd
ibd.top.unique$slope <- ibd.top.unique$ibd
i.top <- full_join(ibd.top.shared, ibd.top.unique)
i.top$Group <- 'UNC'
p.top <- full_join(g.top,i.top)
# violin plot
png(filename = '/work/users/n/n/nnishi/eqtl/freeze/final/figures/violin_plot_mash_shared_vs_unique_eGenes_GTEx_UNC_top_eQTL.png',
    res = 300, units = 'in', height = 5, width = 6)
ggplot(p.top, aes(eGene, abs(slope), fill = Group)) + geom_violin(lwd=1) + 
  geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) + 
  labs(x='', y='mash-Adjusted Absoulte eQTL Effect Size', title='Top eQTL that colocalize with GWAS') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 16)) + 
  geom_signif(y_position=c(1.07, 1.06), xmin=c(0.78, 1.78), xmax=c(1.22, 2.22),
              annotation=p.stats, tip_length=0) 
dev.off()

################################################################################
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
  file <- paste('/work/users/n/n/nnishi/mashr/plots/mashr_violin_plot_coloc_full_signal_',
               genes$name[x], '.png', sep = '')
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
  # file <- paste('/work/users/n/n/nnishi/mashr/plots/mashr_correlation_plot_coloc_full_signal_',
  #              genes$name[x], '.png', sep = '')
  # png(filename = file, res = 300, units = 'in', height = 6, width = 6)
  # plot(df$ibd, df$gtex, xlim = c(-1,1), ylim = c(-1,1),
  #      xlab = 'UNC mash-adjusted effect sizes', ylab = 'GTEx mash-adjusted effect sizes',
  #      main = paste(genes$name[x], 'eQTL signal'))
  # abline(0,1, col='red')
  # abline(v=0,h=0)
  # dev.off()
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
            file = '/work/users/n/n/nnishi/mashr/mashr_coloc_fit_full_signal_effect_size_distribution_stats.txt')
################################################################################
# filter
df.filtered <- df.overlap[df.overlap$wilcoxon_padj < 0.05,] %>% drop_na(wilcoxon_padj)
table(df.filtered$gene_sharing_bin, df.filtered$larger_distribution)

bin <- 'shared'
# pull out shared eGenes and compare GO cluster mapping?? 
shared.gtex <- subset(df.filtered, gene_sharing_bin == bin & larger_distribution == 'GTEx')
shared.ibd <- subset(df.filtered, gene_sharing_bin == bin & larger_distribution == 'IBD')

shared.gtex <- subset(df.filtered, larger_distribution == 'GTEx')
shared.ibd <- subset(df.filtered, larger_distribution == 'IBD')
# add GO cluster assignments
shared.gtex <- merge(shared.gtex, red.out, by.x = 'gene', by.y = 'ID', all.x = TRUE)
shared.ibd <- merge(shared.ibd, red.out, by.x = 'gene', by.y = 'ID', all.x = TRUE)
# label unassigned genes
shared.gtex <- shared.gtex %>% mutate(secondaryTerm = ifelse(GO_ID == primary, primaryTerm, secondaryTerm)) %>% 
  replace_na(list(primaryCluster = 0, primaryTerm = 'unassigned',secondaryCluster = 0, secondaryTerm = 'unassigned'))
shared.ibd <- shared.ibd %>% mutate(secondaryTerm = ifelse(GO_ID == primary, primaryTerm, secondaryTerm)) %>% 
  replace_na(list(primaryCluster = 0, primaryTerm = 'unassigned',secondaryCluster = 0, secondaryTerm = 'unassigned'))
# create a summary table quantifying the number of genes mapping to each category
g <- shared.gtex %>% group_by(primaryTerm, secondaryTerm) %>% summarize(n = n_distinct(gene)) %>% arrange(primaryTerm, secondaryTerm)
i <- shared.ibd %>% group_by(primaryTerm, secondaryTerm) %>% summarize(n = n_distinct(gene)) %>% arrange(primaryTerm, secondaryTerm)
# rename
colnames(g)[1] <- 'GTEx'
colnames(i)[1] <- 'UNC'
# plot
PieDonut(g, aes(GTEx, secondaryTerm, count = n), pieLabelSize = 3, donutLabelSize = 1.75, 
         showRatioThreshold = F, showRatioDonut = F)
PieDonut(i, aes(UNC, secondaryTerm, count = n), pieLabelSize = 3, donutLabelSize = 1.75, 
         showRatioThreshold = F, showRatioDonut = F)
# run stats
# probably not enough genes to run any stats...

g <- shared.gtex %>% group_by(primaryTerm, secondaryTerm) %>% select(name, primaryTerm, secondaryTerm) %>% distinct(name) %>% arrange(primaryTerm, secondaryTerm,name)
i <- shared.ibd %>% group_by(primaryTerm, secondaryTerm) %>% select(name, primaryTerm, secondaryTerm) %>% distinct(name) %>% arrange(primaryTerm, secondaryTerm,name)


# save workspace
#save.image(file = '/work/users/n/n/nnishi/mashr/mashr_coloc_fit_full_signal.RData')
################################################################################
# plot an example
library(ggplot2)
library(ggsignif)

gg <- data.frame(eqtl = c(df$eqtl, df$eqtl), slope = c(df$gtex, df$ibd), group = c(rep('GTEx', nrow(df)), rep('IBD', nrow(df))))

ggplot(gg, aes(group, abs(slope), fill = group)) + geom_violin(lwd=1) +
  geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) +
  labs(x='', y='mash-adjusted absoulte eQTL Effect Size', title='ABO eQTL signal') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  geom_signif(test = 'wilcox.test', comparisons = list(c('GTEx', 'IBD')))

plot(df$ibd, df$gtex, xlim = c(-1,1), ylim = c(-1,1),
     xlab = 'IBD mash-adjusted effect sizes', ylab = 'GTEx mash-adjusted effect sizes',
     main = 'TNFRSF14 eQTL signal')
abline(0,1, col='red')
abline(v=0,h=0)
################################################################################
# LocusZoom plots
library(locuszoomr)
library(LDlinkR)
library(AnnotationHub)
library(memoise)

# LDlink token
t <- '7f59e5f65738'
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
  # if (length(rslist > 1000)) {
  #   rslist <- rslist[order(l$data$logP, decreasing = TRUE)[seq_len(1000)]]
  # }
  rslist <- unique(c(rslist, r))
  ldm <- mem_LDmatrix(rslist, pop = pop, genome_build = genome_build, token = token, ...)
  ld <- ldm[, index_snp]
  l$data$ld <- ld[match(l$data[, labs], ldm$RS_number)]
  l
}
# load coloc SNPs
df.gwas <- fread('/work/users/n/n/nnishi/locuszoom/Liu_GWAS_locuszoom_all_candidate_genes_all_pop_coloc_snps_20240116.txt')
# reformat
df.gwas <- df.gwas %>% unite('id', chr:pos, sep = ':', remove = FALSE)
df.gwas$chr <- gsub('chr', '', df.gwas$chr)
# select gene of interest
egene <- 'ENSG00000157873' #ABO ENSG00000175164 #TNFRSF14 ENSG00000157873
gene <- 'TNFRSF14'

g <- subset(mc.g, gene == egene) %>% distinct(eqtl, .keep_all = TRUE)
i <- subset(mc.i, gene == egene) %>% distinct(eqtl, .keep_all = TRUE)
snps <- intersect(g$snp, i$snp)
df <- subset(g, snp %in% snps) %>% arrange(match(snp, snps))
df <- df %>% separate(snp, c('chr', 'pos', 'ref', 'alt'), remove = FALSE) %>% unite('id', chr:pos, sep = ':', remove = FALSE)
df$chr <- gsub('chr', '', df$chr)
df$abs.gtex <- abs(df$gtex)
df$abs.ibd <- abs(df$ibd)

gwas.gene <- df.gwas[df.gwas$gene == gene,]
# set locus
loc.gwas <- locus(gwas.gene, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', p = 'P-value', 
                  gene = gene, fix_window = 1e5)
loc.gtex <- locus(df, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', yvar = 'abs.gtex', 
                  gene = gene, fix_window = 1e5)
loc.ibd <- locus(df, ens_db = ensembl.v105, labs = 'rsID', chrom = 'chr', yvar = 'abs.ibd', 
                  gene = gene, fix_window = 1e5)
index_snp <- df[df$snp == subset(genes, gene == egene)$GWAS_index_snp,]$rsID
# add LD info
loc.gwas <- get_LD(loc.gwas, token = t, index_snp = index_snp)
loc.gtex <- get_LD(loc.gtex, token = t, index_snp = index_snp)
loc.ibd <- get_LD(loc.ibd, token = t, index_snp = index_snp)
# plot
oldpar <- set_layers(3)
scatter_plot(loc.gwas, index_snp = index_snp, labels = 'index',
             border = TRUE, pcutoff = NULL, cex = 3)
title(gene, adj = 0, line = 2.55)
title('GWAS', adj = 0.2, line = 1)
scatter_plot(loc.ibd, index_snp = index_snp, labels = 'index',
             border = TRUE, legend_pos = NULL, ylim = c(0, round(max(max(df$gtex, df$ibd)), 1)+0.1), cex = 3)
#abline(h = 0, lty = 2)
title('IBD', adj = 0, line = 1)
scatter_plot(loc.gtex, index_snp = index_snp, labels = 'index', 
             border = TRUE, legend_pos = NULL, ylim = c(0, round(max(max(df$gtex, df$ibd)), 1)+0.1), cex = 3)
#abline(h = 0, lty = 2)
title('GTEx', adj = 0, line = 1)
genetracks(loc.ibd)
# plot effect sizes on same panel
oldpar <- set_layers(2)
scatter_plot(loc.gwas, index_snp = index_snp, labels = 'index',
             border = TRUE, pcutoff = NULL, cex = 3)
title(gene, adj = 0, line = 2.55)
title('GWAS', adj = 0.2, line = 1)
scatter_plot(loc.ibd, index_snp = index_snp, labels = 'index',
             border = TRUE, legend_pos = NULL, ylim = c(0, round(max(max(df$gtex, df$ibd)), 1)+0.1), cex = 3)
#abline(h = 0, lty = 2)
title('mash-adjusted effect sizes', adj = 0, line = 1)
scatter_plot(loc.gtex, index_snp = index_snp, labels = 'index', 
             border = TRUE, legend_pos = NULL, ylim = c(0, round(max(max(df$gtex, df$ibd)), 1)+0.1), cex = 3,
             pch = 24, add = TRUE)
#abline(h = 0, lty = 2)
#title('GTEx', adj = 0, line = 1)
genetracks(loc.ibd)









