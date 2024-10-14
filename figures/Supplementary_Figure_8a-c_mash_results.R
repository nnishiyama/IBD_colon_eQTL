# Supplementary Figure 8a-c. mash fit results
library(data.table)
library(tidyverse)
library(mashr)
library(ggplot2)
library(UpSetR)
library(ComplexHeatmap)

# load mash data
load('data/UNC/RData/mashr_colon_eqtl_random_null_corr.RData')
# total number of effects that mash considers significant
ind <- get_significant_results(m.strong)
# extract posterior matrices & subset to mash significant eQTL
mash.beta <- as.data.frame(m.strong$result$PosteriorMean[ind,])
# reformat
mash.beta$eqtl <- row.names(mash.beta)
mash.beta <- mash.beta %>% separate(eqtl, c('gene', 'snp'), sep = '-', remove = FALSE)
################################################################################
# UpSet plot of eQTL overlap
set.seed(12345)
m.beta <- as.data.frame(m.strong$result$PosteriorMean)
m.beta$eqtl <- row.names(m.beta)
# pull out significant results for each condition
m.b <- m.beta[get_significant_results(m.strong, conditions =1),]
m.g <- m.beta[get_significant_results(m.strong, conditions =2),]
m.i <- m.beta[get_significant_results(m.strong, conditions =3),]
# create a set list of eQTL
set.eqtl <- list(UNC = m.i$eqtl, GTEx = m.g$eqtl, BarcUVa = m.b$eqtl)
# create combination matrix
mat.eqtl <- make_comb_mat(set.eqtl, mode = 'distinct')
# UpSet plot
UpSet(mat.eqtl, set_order = c('UNC', 'GTEx', 'BarcUVa'))
png(filename = 'plots/Supplementary_Figure_8a_mash_UpSet_eQTL_overlap.png',
    res = 300, units = 'in', height = 4, width = 6)
ht <- draw(UpSet(mat.eqtl, set_order = c('UNC', 'GTEx', 'BarcUVa'),
                 pt_size = unit(5, 'mm'), lwd = 3,
                 comb_col =c('#7BAFD4','#13294B','#13294B','#13294B','#13294B', '#13294B','#13294B'),
                 top_annotation = upset_top_annotation(mat.eqtl, extend = 0.2, annotation_name_gp = gpar(fontsize=15),axis_param = list(gp=gpar(fontsize=14))),
                 right_annotation = upset_right_annotation(mat.eqtl, extend = 1, 
                                                           annotation_name_gp = gpar(fontsize=15), 
                                                           axis_param = list(gp=gpar(fontsize=14)),
                                                           gp = gpar(col = c('#13294B','#13294B','#13294B'),
                                                                     fill = c('#13294B','#13294B','#13294B'))),
                 row_names_gp = gpar(fontsize=20)), padding = unit(c(0.1,0.1,1.5,1.5), 'cm'))
od <- column_order(ht)
cs <- comb_size(mat.eqtl)
rs <- set_size(mat.eqtl)
ro <- row_order(ht)
decorate_annotation("intersection_size", {grid.text(cs[od], x = seq_along(cs), rot = 45,
                                                    y = unit(cs[od], "native") + unit(2, "pt"), default.units = "native", 
                                                    just = c('left',"bottom"), gp = gpar(fontsize = 20, col = c('#7BAFD4', '#13294B', '#13294B', '#13294B', '#13294B', '#13294B', '#13294B')))})
decorate_annotation('set_size', {grid.text(rs[ro], x = unit(rs[ro], 'native') + unit(30, 'pt'),
                                           y = rev(seq_len(length(rs))), default.units = 'native', just = 'bottom', 
                                           gp = gpar(fontsize = 20, col = c('#13294B','#13294B','#13294B')))})
dev.off()
################################################################################
# compare correlation coefficients before & after mash
# make a combination matrix
combo <- t(combn(names(mash.beta[,1:3]), 2))
comp <- c()
pr <- c()
low <- c()
up <- c()
# loop through combinations
for (i in 1:nrow(combo)) {
  # run correlations pre- & post-mash
  pre <- cor.test(data.strong$Bhat[ind, combo[i,1]], data.strong$Bhat[ind, combo[i, 2]])
  post <- cor.test(mash.beta[ind, combo[i, 1]], mash.beta[ind, combo[i, 2]])
  pr <- c(pr, pre$estimate, post$estimate)
  low <- c(low, pre$conf.int[1], post$conf.int[1])
  up <- c(up, pre$conf.int[2], post$conf.int[2])
  comp <- c(comp, rep(paste(combo[i,1], combo[i,2], sep = '-'), 2))
}
p <- data.frame(comparison = comp, pearson_r = pr, lower = low, upper = up, 
                beta_estimate = rep(c('before mash', 'after mash'), nrow(combo)))
p$comparison <- gsub('barcuva', 'BarcUVa', p$comparison)
p$comparison <- gsub('gtex', 'GTEx', p$comparison)
p$comparison <- gsub('ibd', 'UNC', p$comparison)
# plot
png(filename = 'plots/Supplementary_Figure_8b_mash_effect_size_correlation_before_after.png',
    res = 300, units = 'in', height = 5, width = 5)
ggplot(p, aes(x = comparison, y = pearson_r, color = beta_estimate)) + geom_point(size = 5) + 
  geom_errorbar(aes(ymin = low, ymax = up), width = 0.1, linewidth = 1) + ylim(0.79,0.95) + 
  xlab('Pairwise study comparison') + ylab('eQTL beta estimate correlation (Pearson\'s r)')
dev.off()
################################################################################
# BarcUVa-GTEx scatter plots before & after mash
# pull out overlapping set 
overlap.eqtl <- extract_comb(mat.eqtl, '111')
mash.overlap <- subset(mash.beta, eqtl %in% overlap.eqtl)
# compare distribution of post-mash effect sizes between studies
kruskal.test(c(abs(mash.overlap$barcuva), abs(mash.overlap$gtex), abs(mash.overlap$ibd)), 
             c(rep('barcuva', nrow(mash.overlap)), rep('gtex', nrow(mash.overlap)), rep('ibd', nrow(mash.overlap))))
pairwise.wilcox.test(c(abs(mash.overlap$barcuva), abs(mash.overlap$gtex), abs(mash.overlap$ibd)), 
             c(rep('barcuva', nrow(mash.overlap)), rep('gtex', nrow(mash.overlap)), rep('ibd', nrow(mash.overlap))),
             p.adjust.method = 'BH')
# scatter plot of effect size estimates after mash
png(filename = 'plots/Supplementary_Figure_8c_scatter_plot_mash_effect_size_BarcUVa_GTEx.png',
    res = 300, units = 'in', height = 4, width = 4)
plot(mash.overlap$barcuva, mash.overlap$gtex, xlim = c(-2,2), ylim = c(-2,2), 
     xlab = 'BarcUVa mash effect size estimates', ylab = 'GTEx mash effect size estimates')
abline(v=0, h=0)
abline(0,1, col = 'red')
abline(lm(mash.overlap$gtex ~ mash.overlap$barcuva)$coefficients, col = 'green')
dev.off()
# pull out effect estimates before mash
before <- as.data.frame(data.strong$Bhat)
before$eqtl <- row.names(before)
before <- subset(before, eqtl %in% overlap.eqtl)
# scatter plot of effect size estimates pre-mash
png(filename = 'plots/Supplementary_Figure_8c_scatter_plot_pre-mash_effect_size_BarcUVa_GTEx.png',
    res = 300, units = 'in', height = 4, width = 4)
plot(before$barcuva, before$gtex, xlim = c(-2,2), ylim = c(-2,2), 
     xlab = 'BarcUVa pre-mash effect size estimates', ylab = 'GTEx pre-mash effect size estimates')
abline(v=0, h=0)
abline(0,1, col = 'red')
abline(lm(before$gtex ~ before$barcuva)$coefficients, col = 'green')
dev.off()
