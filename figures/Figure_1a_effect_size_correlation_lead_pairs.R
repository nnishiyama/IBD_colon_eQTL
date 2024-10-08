# Figure 1a. effect size correlation for shared lead eSNP-eGene pairs
library(data.table)
library(tidyverse)
library(qvalue)

# load all eQTL results
ibd <- fread('data/UNC/Supplementary_Table_1_IBD_cis-eQTL_lead.significant.txt')
gtex <- fread('data/GTEx/GTEx_Analysis_v8_eQTL/Colon_Transverse.v8.egenes.txt.gz')
barcuva <- fread('data/BarcUVa/BarcUVA-seq_colon_eQTLs.csv')
# reformat ENSEMBL ID
barcuva <- barcuva %>% separate(`Ensembl gene id (GENCODE v19)`, c('gene', 'ver'))
gtex <- gtex %>% separate(gene_id, c('gene','ver'))
ibd <- ibd %>% separate(phe_id, c('gene', 'ver'))
# reformat SNP
barcuva <- barcuva %>% unite('snp', BED_hg38:ALT, sep = ':', remove = FALSE)
gtex <- gtex %>% unite('snp', chr:alt, sep = ':', remove = FALSE)
ibd$snp <- ibd$var_id
# annotate eQTL
barcuva$eqtl <- paste(barcuva$gene, barcuva$snp, sep = '-')
gtex$eqtl <- paste(gtex$gene, gtex$snp, sep = '-')
ibd$eqtl <- paste(ibd$gene, ibd$snp, sep = '-')
# separate signficant eQTL
barcuva_sig <- barcuva[barcuva$qval < 0.05,]
barcuva_non <- barcuva[barcuva$qval >= 0.05,]
gtex_sig <- gtex[gtex$qval < 0.05,]
gtex_non <- gtex[gtex$qval >= 0.05,]
# find shared pairs
bg.pairs <- intersect(barcuva_sig$eqtl, gtex_sig$eqtl)
ib.pairs <- intersect(ibd$eqtl, barcuva_sig$eqtl)
ig.pairs <- intersect(ibd$eqtl, gtex_sig$eqtl)
# find union and intersection
length(unique(c(ib.pairs, ig.pairs)))
length(intersect(ib.pairs, ig.pairs))
# subset data
bg <- subset(barcuva_sig, eqtl %in% bg.pairs) %>% arrange(match(eqtl, bg.pairs))
bi <- subset(barcuva_sig, eqtl %in% ib.pairs) %>% arrange(match(eqtl, ib.pairs))
gb <- subset(gtex_sig, eqtl %in% bg.pairs) %>% arrange(match(eqtl, bg.pairs))
gi <- subset(gtex_sig, eqtl %in% ig.pairs) %>% arrange(match(eqtl, ig.pairs))
ib <- subset(ibd, eqtl %in% ib.pairs) %>% arrange(match(eqtl, ib.pairs))
ig <- subset(ibd, eqtl %in% ig.pairs) %>% arrange(match(eqtl, ig.pairs))
# correlation test
bg.cor <- cor.test(bg$slope, gb$slope)
ib.cor <- cor.test(bi$slope, ib$slope)
ig.cor <- cor.test(gi$slope, ig$slope)
# plot effect size correlation
## IBD-GTEx
png(filename = 'plots/Figure_1a_effect_size_correlation_lead_pairs_IBD-GTEx.png',
    res = 300, units = 'in', height = 6, width = 6)
plot(ig$slope, gi$slope, pch = 19, xlim = c(-2,2), ylim = c(-2,2), xlab = 'UNC eQTL effect size', ylab = 'GTEx eQTL effect size')
text(-2, 2, paste('r =', round(ig.cor$estimate, 3)), adj = 0)
text(-2, 1.7, paste('p =', signif(ig.cor$p.value, digits = 3)), adj = 0)
abline(coef = c(0,1), col = 'red', lty = 2)
abline(v=0, h=0)
dev.off()
## IBD-BarcUVa
png(filename = 'plots/Figure_1a_effect_size_correlation_lead_pairs_IBD-BarcUVa.png',
    res = 300, units = 'in', height = 6, width = 6)
plot(ib$slope, bi$slope, pch = 19, xlim = c(-2,2), ylim = c(-2,2), xlab = 'UNC eQTL effect size', ylab = 'BarcUVa eQTL effect size')
text(-2, 2, paste('r =', round(ib.cor$estimate, 3)), adj = 0)
text(-2, 1.7, paste('p =', signif(ib.cor$p.value, digits = 3)), adj = 0)
abline(coef = c(0,1), col = 'red', lty = 2)
abline(v=0, h=0)
dev.off()
## GTEx-BarcUVa
png(filename = 'plots/Figure_1a_effect_size_correlation_lead_pairs_GTEx-BarcUVa.png',
    res = 300, units = 'in', height = 6, width = 6)
plot(bg$slope, gb$slope, pch = 19, xlim = c(-2,2), ylim = c(-2,2), ylab = 'GTEx eQTL effect size', xlab = 'Barcuva eQTL effect size')
text(-2, 2, paste('r =', round(bg.cor$estimate, 5)), adj = 0)
text(-2, 1.7, paste('p =', signif(bg.cor$p.value, digits = 3)), adj = 0)
abline(coef = c(0,1), col = 'red', lty = 2)
abline(v=0, h=0)
dev.off()
