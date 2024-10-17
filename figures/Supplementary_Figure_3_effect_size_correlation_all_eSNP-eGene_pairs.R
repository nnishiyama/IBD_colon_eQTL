# correlate effect sizes for all shared significant eSNP-eGene pairs

library(data.table)
library(tidyverse)

# load all significant gene-variant pairs
cond.ibd <- fread('data/UNC/Supplementary_Table_2_IBD_cis-eQTL_conditional.txt')
cond.gtex <- fread('data/GTEx/GTEx_Analysis_v8_eQTL/Colon_Transverse.v8.signif_variant_gene_pairs.txt.gz')
# reformat GTEx data
cond.gtex$variant_id <- gsub('_b38', '', cond.gtex$variant_id)
cond.gtex$variant_id <- gsub('_', ':', cond.gtex$variant_id)
# pull out significant BarcUVa variants since no conditional results
eqtl.barcuva <- fread('data/BarcUVa/barcuvaseq.eqtls.allpairs.sumstats.hg38.txt')
barcuva <- fread('data/BarcUVa/BarcUVA-seq_colon_eQTLs.csv')
barcuva <- barcuva[barcuva$qval < 0.05,]
cond.barcuva <- data.frame()
# loop through eQTL
for (i in 1:nrow(barcuva)) {
  print(i)
  b <- eqtl.barcuva[eqtl.barcuva$gene_id == barcuva$`Ensembl gene id (GENCODE v19)`[i] & eqtl.barcuva$pval_nominal <= barcuva$pval_nominal_threshold[i], ]
  cond.barcuva <- rbind(cond.barcuva, b)
}
# reformat ENSEMBL IDs
cond.barcuva <- cond.barcuva %>% separate(gene_id, c('gene', 'ver'))
cond.gtex <- cond.gtex %>% separate(gene_id, c('gene', 'ver'))
cond.ibd <- cond.ibd %>% separate(phe_id, c('gene', 'ver'))
# annotate eqtl
cond.barcuva$eqtl <- paste(cond.barcuva$gene, cond.barcuva$variant_id, sep = '-')
cond.gtex$eqtl <- paste(cond.gtex$gene, cond.gtex$variant_id, sep = '-')
cond.ibd$eqtl <- paste(cond.ibd$gene, cond.ibd$var_id, sep = '-')
# find shared eSNP-eGene pairs
pairs.barcuva <- intersect(cond.ibd$eqtl, cond.barcuva$eqtl)
pairs.gtex <- intersect(cond.ibd$eqtl, cond.gtex$eqtl)
pairs.nibd <- intersect(cond.barcuva$eqtl, cond.gtex$eqtl)
# subset data
cond.barcuva.g <- subset(cond.barcuva, eqtl %in% pairs.nibd) %>% arrange(match(eqtl, pairs.nibd))
cond.gtex.b <- subset(cond.gtex, eqtl %in% pairs.nibd) %>% arrange(match(eqtl, pairs.nibd))
cond.barcuva <- subset(cond.barcuva, eqtl %in% pairs.barcuva) %>% arrange(match(eqtl, pairs.barcuva))
cond.gtex <- subset(cond.gtex, eqtl %in% pairs.gtex) %>% arrange(match(eqtl, pairs.gtex))
cond.ibd.b <- subset(cond.ibd, eqtl %in% pairs.barcuva) %>% arrange(match(eqtl, pairs.barcuva))
cond.ibd.g <- subset(cond.ibd, eqtl %in% pairs.gtex) %>% arrange(match(eqtl, pairs.gtex))
# drop duplicates
cond.ibd.b <- cond.ibd.b %>% distinct(eqtl, .keep_all = TRUE)
cond.ibd.g <- cond.ibd.g %>% distinct(eqtl, .keep_all = TRUE)
all(cond.barcuva$eqtl == cond.ibd.b$eqtl)
all(cond.gtex$eqtl == cond.ibd.g$eqtl)
# run pearson's correlation on effect sizes
bg.cor <- cor.test(cond.barcuva.g$slope, cond.gtex.b$slope)
ib.cor <- cor.test(cond.barcuva$slope, cond.ibd.b$bwd_slope)
ig.cor <- cor.test(cond.gtex$slope, cond.ibd.g$bwd_slope)
# plot
## IBD-GTEx
png(filename = 'plots/Supplementary_Figure_3a_Correlation_all_pairs_IBD-GTEx.png',
    res = 300, units = 'in', height = 6, width = 6)
plot(cond.ibd.g$bwd_slope, cond.gtex$slope, xlab = 'UNC eQTL effect size', 
     ylab = 'GTEx eQTL effect size', xlim = c(-2, 2), ylim = c(-2,2))
text(-2, 2, paste('r =', round(ig.cor$estimate, 3)), adj = 0)
text(-2, 1.7, paste('p =', signif(ig.cor$p.value, digits = 3)), adj = 0)
abline(v = 0, h = 0)
abline(c(0,1), lty = 2, col = 'red')
dev.off()
## IBD-BarcUVa
png(filename = 'plots/Supplementary_Figure_3b_Correlation_all_pairs_IBD-BarcUVa.png',
    res = 300, units = 'in', height = 6, width = 6)
plot(cond.ibd.b$bwd_slope, cond.barcuva$slope, xlab = 'UNC eQTL effect size', 
     ylab = 'BarcUVa eQTL effect size', xlim = c(-2, 2), ylim = c(-2,2))
text(-2, 2, paste('r =', round(ib.cor$estimate, 3)), adj = 0)
text(-2, 1.7, paste('p =', signif(ib.cor$p.value, digits = 3)), adj = 0)
abline(v = 0, h = 0)
abline(c(0,1), lty = 2, col = 'red')
dev.off()
## GTEx-BarcUVa
png(filename = 'plots/Supplementary_Figure_3c_Correlation_all_pairs_GTEx-BarcUVa.png',
    res = 300, units = 'in', height = 6, width = 6)
plot(cond.barcuva.g$slope, cond.gtex.b$slope, ylab = 'GTEx eQTL effect size', 
     xlab = 'BarcUVa eQTL effect size', xlim = c(-2, 2), ylim = c(-2,2))
text(-2, 2, paste('r =', round(bg.cor$estimate, 3)), adj = 0)
text(-2, 1.7, paste('p =', signif(bg.cor$p.value, digits = 3)), adj = 0)
abline(v = 0, h = 0)
abline(c(0,1), lty = 2, col = 'red')
dev.off()
