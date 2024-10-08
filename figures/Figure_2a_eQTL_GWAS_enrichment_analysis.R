# Figure 2a. enrichment heatmap of eVariants in GWAS loci
library(data.table)
library(tidyverse)
library(pheatmap)

load('data/UNC/coloc_all_data.RData')

# load all significant gene-variant pairs
cond.ibd <- fread('data/UNC/Supplementary_Table_2_IBD_cis-eQTL_conditional.txt')
cond.gtex <- fread('data/GTEx/GTEx_Analysis_v8_eQTL/Colon_Transverse.v8.signif_variant_gene_pairs.txt.gz')
# reformat GTEx data
cond.gtex$variant_id <- gsub('_b38', '', cond.gtex$variant_id)
cond.gtex$variant_id <- gsub('_', ':', cond.gtex$variant_id)
# separate out significant GWAS SNPs
sig.cd <- cd[cd$`P-value` < 5e-8,]
sig.ibd <- ibd[ibd$`P-value` < 5e-8,]
sig.uc <- uc[uc$`P-value` < 5e-8,]
################################################################################
# basic enrichment analysis based off of GTEx
# pull variants tested for eQTL
ibd.snps <- unique(eqtl.ibd$var_id)
gtex.snps <- unique(eqtl.gtex$variant_id)
# pull significant eSNPs
ibd.snps.sig <- intersect(ibd.snps, cond.ibd$var_id)
gtex.snps.sig <- intersect(gtex.snps, cond.gtex$variant_id)
# calculate enrichment
es.cd <- (length(intersect(sig.cd$MarkerName, ibd.snps.sig))/length(intersect(sig.cd$MarkerName, ibd.snps)))/(length(ibd.snps.sig)/length(ibd.snps))
es.ibd <- (length(intersect(sig.ibd$MarkerName, ibd.snps.sig))/length(intersect(sig.ibd$MarkerName, ibd.snps)))/(length(ibd.snps.sig)/length(ibd.snps))
es.uc <- (length(intersect(sig.uc$MarkerName, ibd.snps.sig))/length(intersect(sig.uc$MarkerName, ibd.snps)))/(length(ibd.snps.sig)/length(ibd.snps))
es.i <- c('Colon_UNC', es.cd, es.ibd, es.uc)
################################################################################
# pull out significant BarcUVa variants since no conditional results
barcuva <- fread('data/BarcUVa/BarcUVA-seq_colon_eQTLs.csv')
barcuva <- barcuva[barcuva$qval < 0.05,]
b.snps.sig <- c()
# loop through eQTL
for (i in 1:nrow(barcuva)) {
  print(i)
  b <- eqtl.barcuva[eqtl.barcuva$gene_id == barcuva$`Ensembl gene id (GENCODE v19)`[i] & eqtl.barcuva$pval_nominal <= barcuva$pval_nominal_threshold[i], ]$variant_id
  b.snps.sig <- c(b.snps.sig, b)
  remove(b)
}
b.snps.sig <- unique(b.snps.sig)
b.snps <- unique(eqtl.barcuva$variant_id)
# calculate enrichment
es.cd <- (length(intersect(sig.cd$MarkerName, b.snps.sig))/length(intersect(sig.cd$MarkerName, b.snps)))/(length(b.snps.sig)/length(b.snps))
es.ibd <- (length(intersect(sig.ibd$MarkerName, b.snps.sig))/length(intersect(sig.ibd$MarkerName, b.snps)))/(length(b.snps.sig)/length(b.snps))
es.uc <- (length(intersect(sig.uc$MarkerName, b.snps.sig))/length(intersect(sig.uc$MarkerName, b.snps)))/(length(b.snps.sig)/length(b.snps))
es.b <- c('Colon_BarcUVa', es.cd, es.ibd, es.uc)
remove(cd, ibd, uc, eqtl.barcuva, eqtl.gtex, eqtl.ibd)
################################################################################
# calculate chi-squared stats for colon eQTL
# CD
cd.df <- data.frame(analysis = c('ibd', 'ibd', 'gtex', 'gtex', 'barcuva', 'barcuva'),
                    shared = c('yes', 'no', 'yes', 'no', 'yes', 'no'),
                    eqtl = c(length(intersect(ibd.snps.sig, sig.cd$MarkerName)),
                             length(setdiff(ibd.snps.sig, sig.cd$MarkerName)),
                             length(intersect(gtex.snps.sig, sig.cd$MarkerName)),
                             length(setdiff(gtex.snps.sig, sig.cd$MarkerName)),
                             length(intersect(b.snps.sig, sig.cd$MarkerName)),
                             length(setdiff(b.snps.sig, sig.cd$MarkerName))))
chisq.test(xtabs(eqtl~analysis+shared, data = cd.df))
# IBD
ibd.df <- data.frame(analysis = c('ibd', 'ibd', 'gtex', 'gtex', 'barcuva', 'barcuva'),
                     shared = c('yes', 'no', 'yes', 'no', 'yes', 'no'),
                     eqtl = c(length(intersect(ibd.snps.sig, sig.ibd$MarkerName)),
                              length(setdiff(ibd.snps.sig, sig.ibd$MarkerName)),
                              length(intersect(gtex.snps.sig, sig.ibd$MarkerName)),
                              length(setdiff(gtex.snps.sig, sig.ibd$MarkerName)),
                              length(intersect(b.snps.sig, sig.ibd$MarkerName)),
                              length(setdiff(b.snps.sig, sig.ibd$MarkerName))))
chisq.test(xtabs(eqtl~analysis+shared, data = ibd.df))
# UC
uc.df <- data.frame(analysis = c('ibd', 'ibd', 'gtex', 'gtex', 'barcuva', 'barcuva'),
                    shared = c('yes', 'no', 'yes', 'no', 'yes', 'no'),
                    eqtl = c(length(intersect(ibd.snps.sig, sig.uc$MarkerName)),
                             length(setdiff(ibd.snps.sig, sig.uc$MarkerName)),
                             length(intersect(gtex.snps.sig, sig.uc$MarkerName)),
                             length(setdiff(gtex.snps.sig, sig.uc$MarkerName)),
                             length(intersect(b.snps.sig, sig.uc$MarkerName)),
                             length(setdiff(b.snps.sig, sig.uc$MarkerName))))
chisq.test(xtabs(eqtl~analysis+shared, data = uc.df))
################################################################################
# loop through GTEx files
sig.files <- dir(path = 'data/GTEx/GTEx_Analysis_v8_eQTL', pattern = '*.v8.signif_variant_gene_pairs.txt.gz', full.names = TRUE)
full.files <- dir(path = 'data/GTEx/GTEx_Analysis_v8_eQTL_full_summary_statistics', pattern = '*.allpairs.txt.gz', full.names = TRUE)
df <- rbind(es.i, es.b)
for (i in 1:length(full.files)){
  print(i)
  es <- c()
  print(full.files[i])
  full.eqtl <- fread(full.files[i])
  full.eqtl <- full.eqtl %>% na.omit()
  full.eqtl$variant_id <- gsub('_b38', '', full.eqtl$variant_id)
  full.eqtl$variant_id <- gsub('_', ':', full.eqtl$variant_id)
  print(sig.files[i])
  sig.eqtl <- fread(sig.files[i])
  sig.eqtl$variant_id <- gsub('_b38', '', sig.eqtl$variant_id)
  sig.eqtl$variant_id <- gsub('_', ':', sig.eqtl$variant_id)
  tissue <- tail(strsplit(strsplit(sig.files[i], split = '[.]')[[1]][1], split = '[/]')[[1]],1)
  snps <- unique(full.eqtl$variant_id)
  snps.sig <- intersect(snps, sig.eqtl$variant_id)
  # calculate enrichment score
  es.cd <- (length(intersect(sig.cd$MarkerName, snps.sig))/length(intersect(sig.cd$MarkerName, snps)))/(length(snps.sig)/length(snps))
  es.ibd <- (length(intersect(sig.ibd$MarkerName, snps.sig))/length(intersect(sig.ibd$MarkerName, snps)))/(length(snps.sig)/length(snps))
  es.uc <- (length(intersect(sig.uc$MarkerName, snps.sig))/length(intersect(sig.uc$MarkerName, snps)))/(length(snps.sig)/length(snps))
  es <- c(tissue, es.cd, es.ibd, es.uc)
  df <- rbind(df, es)
  remove(full.eqtl, sig.eqtl, snps, snps.sig)
}
df <- as.data.frame(df)
colnames(df) <- c('tissue', 'CD', 'IBD', 'UC')
row.names(df) <- df$tissue
row.names(df)[3:9] <- paste(df$tissue[3:9], '_GTEx', sep ='')
df$CD <- as.numeric(df$CD)
df$IBD <- as.numeric(df$IBD)
df$UC <- as.numeric(df$UC)
# reorder
df <- df %>% mutate(colon = ifelse(tissue %like% 'Colon', 'colon', tissue))
df <- df[order(df$colon != 'colon', df$colon),]
# plot
png(res = 300, units = 'in', height = 4, width = 10,
    filename = 'plots/Figure_2a_GWAS_enrichment_eVariants_heatmap.png')
pheatmap(df[,2:4], cluster_rows = F, cluster_cols = F, border_color = 'white', 
         display_numbers = TRUE, number_color = 'black', fontsize_number = 15, fontsize = 15,
         legend_breaks = c(0, 5, 10, 15), main = 'Enrichment of eVariants in IBD GWAS loci',
         cellwidth = 140, cellheight = 20)
dev.off()
