# Figure 4b. mash effect size distribution violin plot
library(data.table)
library(tidyverse)
library(mashr)
library(ggplot2)
library(UpSetR)
library(ComplexHeatmap)
library(ggsignif)

# load mash data
load('data/UNC/RData/mashr_colon_eqtl_random_null_corr.RData')
# set up combination matrix
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
# pull out overlapping set 
overlap.eqtl <- extract_comb(mat.eqtl, '111')
mash.overlap <- subset(mash.beta, eqtl %in% overlap.eqtl)
# pull out effect estimates before mash
before <- as.data.frame(data.strong$Bhat)
before$eqtl <- row.names(before)
before <- subset(before, eqtl %in% overlap.eqtl)
# compare effect size distributions before & after mash
wb <- wilcox.test(abs(before$barcuva), abs(mash.overlap$barcuva))
wg <- wilcox.test(abs(before$gtex), abs(mash.overlap$gtex))
wi <- wilcox.test(abs(before$ibd), abs(mash.overlap$ibd))
wbg <- wilcox.test(abs(mash.overlap$barcuva), abs(mash.overlap$gtex))
wbi <- wilcox.test(abs(mash.overlap$barcuva), abs(mash.overlap$ibd))
wgi <- wilcox.test(abs(mash.overlap$gtex), abs(mash.overlap$ibd))
p.stat <- c(signif(wb$p.value, digits = 3), signif(wg$p.value, digits = 3), round(wi$p.value, 3))
p.stat <- gsub('^', 'p = ', p.stat)
p.stat <- c(p.stat, rep('p < 2.2e-16', 3))
# plot distrubtions before & after mash
# format for plotting
mo.b <- mash.overlap %>% select(eqtl, barcuva)
mo.g <- mash.overlap %>% select(eqtl, gtex)
mo.i <- mash.overlap %>% select(eqtl, ibd)
colnames(mo.b) <- c('eqtl', 'beta')
colnames(mo.g) <- c('eqtl', 'beta')
colnames(mo.i) <- c('eqtl', 'beta')
mo.b$group <- 'BarcUVa'
mo.g$group <- 'GTEx'
mo.i$group <- 'UNC'
ba <- full_join(mo.b, mo.g)
ba <- full_join(ba, mo.i)
ba$beta_estimate <- 'after mash'
b.b <- before %>% select(eqtl, barcuva)
b.g <- before %>% select(eqtl, gtex)
b.i <- before %>% select(eqtl, ibd)
colnames(b.b) <- c('eqtl', 'beta')
colnames(b.g) <- c('eqtl', 'beta')
colnames(b.i) <- c('eqtl', 'beta')
b.b$group <- 'BarcUVa'
b.g$group <- 'GTEx'
b.i$group <- 'UNC'
b <- full_join(b.b, b.g)
b <- full_join(b, b.i)
b$beta_estimate <- 'before mash'
ba <- full_join(ba, b)
ba$beta_estimate <- factor(ba$beta_estimate, levels = c('before mash', 'after mash'))
# summary stats
ss <- ba %>% group_by(group, beta_estimate) %>% summarize(med=round(median(abs(beta)), 3), ym = -0.02)
# plot
png(filename = 'plots/Figure_4b_violin_plot_mash_effect_size_before_after.png',
    res = 300, units = 'in', height = 8, width = 8)
ggplot(ba, aes(group, abs(beta), fill = beta_estimate)) + geom_violin(lwd=1) + 
  geom_boxplot(outlier.shape = NA, width = 0.1, position = position_dodge(width = 0.9)) + 
  labs(x='', y='Absoulte eQTL Effect Size', title='Shared mash eQTL') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  geom_signif(y_position=c(2.19, 1.9, 1.75, 2.07, 2.16), 
              xmin=c(0.78, 1.78, 2.78, 1.22, 1.22), 
              xmax=c(1.22, 2.22, 3.22, 2.22, 3.22),
              annotation=p.stat[1:5], tip_length=0) +
  geom_signif(y_position=c(1.85), 
              xmin=c(2.22), 
              xmax=c(3.22),
              annotation=p.stat[6], tip_length=0, color = 'red', fontface = 'bold') +
  geom_text(data = ss, aes(group, y = ym, label = paste0('med = ', med)), 
            position = position_dodge(width = 0.9), size = 3)
dev.off()
